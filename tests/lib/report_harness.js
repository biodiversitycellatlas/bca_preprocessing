#!/usr/bin/env node
/*
 * Run filtering_report.html's own script outside a browser, apply a set of choices to
 * its state, and print the selection its export button would download.
 *
 * The page does the filtering itself, so the selection apply_metacell_filter.py reads
 * is computed in JavaScript; tests/checks/metacell_filtering.sh feeds this output to
 * the apply step and checks the result against the fixture. The DOM is stubbed just
 * far enough for the page to initialise; Plotly is a no-op.
 *
 * Usage: report_harness.js REPORT.html CHOICES.json [SUMMARY_OUT.json]
 *
 * SUMMARY_OUT, when given, receives the HTML of the summary panel (statistics,
 * settings and running totals) as the page rendered it for the first sample.
 *
 * CHOICES: { "samples": { "<id>": { maxMito, maxDoublet, maxAlpha, blacklist: [],
 *            cells: { ... }, genes: { ... } } },
 *            "genes": { minUmis, pfam: [], pfamKeep: [], manual } }
 *
 * A sample's "cells" are its per-cell rules (minUmis, excludeDoublets, ...). The top-level
 * "genes" go to every sample named, under its own "genes". Every sample named counts as
 * set at every level; the others stay pending.
 *
 * SUMMARY_OUT also receives "flow": what the selection holds with a level pending, and what
 * "Apply settings to all samples" copies.
 */
"use strict";

const fs = require("fs");
const vm = require("vm");

const [htmlPath, choicesPath, summaryPath] = process.argv.slice(2);
if (!htmlPath || !choicesPath) {
    console.error("Usage: report_harness.js REPORT.html CHOICES.json");
    process.exit(2);
}

const html = fs.readFileSync(htmlPath, "utf8");
const choices = JSON.parse(fs.readFileSync(choicesPath, "utf8"));

const dataMatch = html.match(/<script id="report-data" type="application\/json">([\s\S]*?)<\/script>/);
if (!dataMatch) throw new Error("no report-data block");
const scripts = [...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map(m => m[1]);
if (scripts.length !== 1) throw new Error(`expected one inline script, found ${scripts.length}`);

function element(id) {
    return {
        id, value: "", checked: false, disabled: false, textContent: "", innerHTML: "",
        dataset: {}, style: {}, files: null, max: "100", min: "0",
        classList: { toggle() {}, add() {}, remove() {}, contains() { return false; } },
        addEventListener() {}, querySelector() { return element(); }, querySelectorAll() { return []; },
        insertAdjacentHTML() {}, appendChild() {}, remove() {}, click() {}, closest() { return null; },
        setAttribute() {}, getAttribute() { return null; },
        // Plotly's event API, recorded so drags and clicks on the plots can be replayed
        _handlers: {},
        on(event, fn) { this._handlers[event] = fn; },
    };
}

const elements = new Map();
const document = {
    activeElement: null,
    body: element("body"),
    getElementById(id) {
        if (id === "report-data") return { textContent: dataMatch[1] };
        if (!elements.has(id)) elements.set(id, element(id));
        return elements.get(id);
    },
    querySelector() { return element(); },
    querySelectorAll() { return []; },
    createElement() { return element(); },
};

const context = vm.createContext({
    document,
    window: {},
    console,
    Plotly: { react() {}, purge() {} },
    requestAnimationFrame: () => 0,   // renders are driven explicitly below
    setTimeout, URL, Blob: function () {}, navigator: {}, confirm: () => true,
});
context.window = context;

vm.runInContext(scripts[0], context, { filename: "filtering_report.html" });

// The page's top-level bindings live in the context's global lexical scope
context.__choices = choices;
context.__summary = null;
const selection = vm.runInContext(`
    (() => {
        const c = __choices;
        for (const [id, s] of Object.entries(c.samples || {})) {
            const st = state.perSample[id];
            if (!st) throw new Error("no such sample in the report: " + id);
            for (const k of ["maxMito", "maxDoublet", "maxAlpha", "minMedianUmis", "minMaxLfc", "minFcGenes"]) {
                if (k in s) st[k] = s[k];
            }
            if (s.blacklist) st.blacklist = new Set(s.blacklist);
            Object.assign(st.cells, s.cells || {});
            const g = Object.assign({}, c.genes || {}, s.genes || {});
            if ("minUmis" in g) st.genes.minUmis = g.minUmis;
            if (g.pfam) st.genes.pfam = new Set(g.pfam);
            if (g.pfamKeep) st.genes.pfamKeep = new Set(g.pfamKeep);
            if ("manual" in g) st.genes.manual = g.manual;
            st.set = { cells: true, metacells: true, genes: true };
        }

        // The summary panel and every tab render without throwing on these choices, the
        // co-expression views and the metacell map's colourings included
        renderSummary();
        for (const tab of ["cells", "metacells", "genes", "export"]) { state.tab = tab; renderActive(); }
        for (const axes of [["mito_pct", "max_lfc"], ["median_umis", "mito_pct"]]) {
            _qcAxes.splice(0, 2, ...axes);
            for (const by of ["x", "y", "gene", "status"]) { _colourBy = by; state.tab = "metacells"; renderActive(); }
        }
        // A selected domain narrows the most expressed genes to its own, and the co-expression to one of them
        _domain = { u: current().universe, d: D.domains.indexOf("PF00125") };
        state.tab = "genes"; renderActive();
        if (!listedGenes(current()).includes(selectedGene(current()))) throw new Error("the co-expression left the selected domain");
        _domain = null;
        __summary = {
            settings: ["cells", "metacells", "genes"].map(l => document.getElementById("settings-" + l).innerHTML).join(""),
            totals: document.getElementById("totals-table").innerHTML,
            stats: document.getElementById("stats-table").innerHTML,
        };

        const sel = buildSelection();
        // A selection loaded back must restore the same choices
        const reloaded = loadSelection(JSON.parse(JSON.stringify(sel)));
        if (reloaded.problems.length) throw new Error("reload problems: " + reloaded.problems.join("; "));
        const again = buildSelection();
        again.created = sel.created;
        if (JSON.stringify(again) !== JSON.stringify(sel)) throw new Error("selection changed after a reload");

        // The plots are the threshold handles: replay drags of the mito line on the QC plot
        // (one across, one along) and of a size line (log axis, data units), a drag of the rank
        // plot's line (log axis, data units) and a click on a gene-histogram bar, then the two
        // lines of the UMIs-per-cell histogram (log10 axis). Done after the selection is built,
        // so they do not change it; on S1_starsolo, whose mito threshold is on.
        const fire = (id, event, payload) => {
            const h = document.getElementById(id)._handlers[event];
            if (!h) throw new Error("no " + event + " handler on " + id);
            h(payload);
        };
        // The QC plot's shape of a threshold's line, on the axis it is drawn along
        const qcShape = (rule, axis) => {
            const i = document.getElementById("plot-qc")._lines.findIndex(l => l.rule === rule && l.axis === axis);
            if (i < 0) throw new Error("no " + axis + " line of " + rule + " on the QC plot");
            return (end) => "shapes[" + i + "]." + end;
        };
        const firstSample = state.sample;
        state.sample = "S1_starsolo";
        sstate().minCells = 2;
        // Mito across and size along; then mito along
        const qcAxes = _qcAxes.slice();
        _qcAxes.splice(0, 2, "mito_pct", "n_cells");
        state.tab = "metacells"; renderActive();
        let sh = qcShape("maxMito", "x");
        fire("plot-qc", "plotly_relayout", { [sh("x0")]: 12.34, [sh("x1")]: 12.38, [sh("y0")]: 0.2 });
        const mito = sstate().maxMito;
        renderActive();
        sh = qcShape("minCells", "y");
        fire("plot-qc", "plotly_relayout", { [sh("y0")]: 4, [sh("y1")]: 9, [sh("x0")]: 0.5, [sh("x1")]: 0.5 });
        const size = sstate().minCells;
        _qcAxes.splice(0, 2, "median_umis", "mito_pct");
        renderActive();
        sh = qcShape("maxMito", "y");
        fire("plot-qc", "plotly_relayout", { [sh("y0")]: 15.0, [sh("y1")]: 15.2, [sh("x0")]: 0.3 });
        const mitoAlong = sstate().maxMito;
        Object.assign(sstate(), { maxMito: mito, minCells: null });
        _qcAxes.splice(0, 2, ...qcAxes);
        state.tab = "genes"; renderActive();
        fire("plot-rank", "plotly_relayout", { "shapes[0].y0": 50, "shapes[0].y1": 50 });
        const rank = gstate().minUmis;
        fire("plot-genehist", "plotly_click", { points: [{ x: 2 }] });
        const click = gstate().minUmis;
        Object.assign(cstate(), { minUmis: 10, maxUmis: 1000 });
        state.tab = "cells"; renderActive();
        fire("plot-cumi", "plotly_relayout", { "shapes[1].x0": 2, "shapes[1].x1": 2 });
        renderActive();
        fire("plot-cumi", "plotly_relayout", { "shapes[0].x0": 1.5, "shapes[0].x1": 1.5 });
        __summary.drag = { mito, mitoAlong, size, rank, click, cellMin: cstate().minUmis, cellMax: cstate().maxUmis };

        // What a setting excludes is not offered for blacklisting, by any route: M1.2 fails the
        // mito threshold, GENE011 the minimum UMIs. A kept metacell still can be.
        const before = { black: [...sstate().blacklist], manual: gstate().manual };
        markMetacells(["M1.2"], "add");
        markMetacells(["M1.2"], "toggle");
        markGene(D.universes[current().universe].ids.indexOf("GENE011"), true);
        const blocked = { metacell: !sstate().blacklist.has("M1.2"), gene: gstate().manual === before.manual };
        markMetacells(["M2.3"], "toggle");
        blocked.kept_allowed = sstate().blacklist.has("M2.3");
        sstate().blacklist = new Set(before.black);

        // The gene table's search takes a regular expression, field by field, and "Blacklist
        // matching" lists the kept genes it matches: ^gapdh$ is Gapdh (GENE002) alone; ^ge also
        // matches GENE001 and GENE011 and the rest by id, but those two are excluded already
        document.getElementById("gene-filter").value = "^gapdh$";
        state.tab = "genes"; renderActive();
        document.getElementById("gene-black-matching").onclick();
        blocked.matching = splitList(gstate().manual).filter(id => !splitList(before.manual).includes(id));
        gstate().manual = before.manual;
        document.getElementById("gene-filter").value = "^ge";
        renderActive();
        blocked.matching_skips_excluded = !_geneMatches.some(j => ["GENE001", "GENE011"].includes(D.universes[current().universe].ids[j]));
        document.getElementById("gene-filter").value = "";
        __summary.blocked = blocked;

        // Co-expression, on S1_starsolo: its two kept metacells are too few, so every real one
        // counts. The cut-off line drags (clamped to 0..1, to 0.01) on the co-expressed genes and
        // on the curve (where y is r); a curve point, a co-expressed gene's bar, a top gene's bar
        // and a heatmap square (its column) each open their gene; r agrees with a direct
        // Pearson; a top gene's share over one metacell is its CP10k there; and "Blacklist
        // correlated" adds kept genes only.
        const s1 = current();
        _gene = null;
        state.tab = "genes"; renderActive();
        fire("plot-corr", "plotly_relayout", { "shapes[0].x0": 0.456, "shapes[0].x1": 0.456 });
        const coexp = { cut: _corrMin, fallback: corrMetacells(s1, metacellCalls(s1, sstate())).fallback };
        fire("plot-corr", "plotly_relayout", { "shapes[0].x0": 1.7, "shapes[0].x1": 1.7 });
        coexp.clamped = _corrMin;
        fire("plot-curve", "plotly_relayout", { "shapes[0].y0": 0.333, "shapes[0].y1": 0.333 });
        coexp.curve_cut = _corrMin;
        // The band's upper end drags too (index 1); the lower end pushes it up when it passes it
        fire("plot-corr", "plotly_relayout", { "shapes[1].x0": 0.95, "shapes[1].x1": 0.95 });
        const upper = _corrMax;
        fire("plot-corr", "plotly_relayout", { "shapes[0].x0": 0.99, "shapes[0].x1": 0.99 });
        coexp.band = upper === 0.95 && _corrMin === 0.99 && _corrMax === 0.99;
        _corrMax = 1;
        renderActive();
        const curve = document.getElementById("plot-curve")._corr;
        fire("plot-curve", "plotly_click", { points: [{ curveNumber: 0, pointIndex: 1 }] });
        coexp.curve_open = _gene.j === curve[1].j;
        _corrMin = 0;
        renderActive();
        const shown = document.getElementById("plot-corr")._shown;
        fire("plot-corr", "plotly_click", { points: [{ y: 0 }] });
        coexp.bar_open = _gene.j === shown[0];
        renderActive();
        const idx = s1.expression.mc;
        const P = logProfiles(s1);
        const top = correlatedWith(s1, selectedGene(s1), idx)[0];
        const pearson = (x, y) => {
            const mx = x.reduce((a, v) => a + v, 0) / x.length, my = y.reduce((a, v) => a + v, 0) / y.length;
            let sxy = 0, sxx = 0, syy = 0;
            x.forEach((v, t) => { sxy += (v - mx) * (y[t] - my); sxx += (v - mx) ** 2; syy += (y[t] - my) ** 2; });
            return sxy / Math.sqrt(sxx * syy);
        };
        const direct = pearson(idx.map(i => P.rows[P.index.get(selectedGene(s1))][i]), idx.map(i => P.rows[P.index.get(top.j)][i]));
        coexp.pearson = Math.abs(direct - top.r) < 1e-4;
        const hm = document.getElementById("plot-coexp")._genes;
        fire("plot-coexp", "plotly_click", { points: [{ curveNumber: 1, x: 0, y: 1 }] });
        coexp.heatmap = _gene.j === hm[0];
        fire("plot-coexp", "plotly_click", { points: [{ curveNumber: 0, x: 1.5, y: 2 }] });
        coexp.heatmap = coexp.heatmap && _gene.j === hm[2];
        // The strip's finder opens a gene by name: Gapdh is GENE002
        const pick = document.getElementById("gene-pick");
        pick.value = "Gapdh";
        pick.onchange({ target: pick });
        coexp.pick = D.universes[s1.universe].ids[selectedGene(s1)] === "GENE002" && pick.value === "";
        // The gene table's top co-expressed genes are the curve's first ones, in a column of their own
        const tableTop = topCoexpressed(s1, idx);
        coexp.table = [...tableTop.keys()].every(j => {
            const c = correlatedWith(s1, j, idx);
            return tableTop.get(j).map(x => x.j).join() === c.slice(0, 3).map(x => x.j).join();
        }) && document.getElementById("gene-table").innerHTML.includes("Top co-expressed");
        // The selected domain's head over the genes it lists: S1_starsolo blacklists Ribosomal_L3,
        // whose button takes it off the blacklist and puts it back
        _domain = { u: s1.universe, d: D.domains.indexOf("Ribosomal_L3") };
        renderActive();
        const domButton = document.getElementById("dom-toggle");
        const domBefore = domButton.textContent;
        domButton.onclick();
        const off = !gstate().pfam.has("Ribosomal_L3");
        renderActive();
        const domAfter = domButton.textContent;
        domButton.onclick();
        coexp.domain = domBefore === "Remove from blacklist" && off && domAfter === "Blacklist domain" && gstate().pfam.has("Ribosomal_L3");
        _domain = null;
        renderActive();
        const one = metacellShares(s1, [idx[0]]).get(hm[0]), cp10k = geneProfile(s1, hm[0])[idx[0]];
        coexp.shares = Math.abs(one - cp10k / 100) <= 1e-4 * Math.max(1, cp10k / 100);
        renderActive();
        const keptBefore = geneCalls(s1).reasons.map((r, j) => r.length ? null : D.universes[s1.universe].ids[j]).filter(Boolean);
        const manualBefore = gstate().manual;
        document.getElementById("corr-black-module").onclick();
        const added = splitList(gstate().manual).filter(id => !splitList(manualBefore).includes(id));
        coexp.module = added.length > 0 && added.every(id => keptBefore.includes(id));
        gstate().manual = manualBefore;
        // A gene let go from the strip stays let go, a domain selected or not: the card asks for
        // one, nothing is offered for blacklisting, and the metacell map still draws
        _gene = { u: s1.universe, j: null };
        renderActive();
        const none = selectedGene(s1) === null && document.getElementById("nav-gene").innerHTML === "none" &&
            !_corrModule.length && document.getElementById("corr-empty").textContent.startsWith("Click a gene");
        _domain = { u: s1.universe, d: D.domains.indexOf("PF00125") };
        renderActive();
        coexp.deselect = none && selectedGene(s1) === null;
        // A plot cleared loses its listeners (Plotly.purge), so it binds them again when drawn:
        // the curve, cleared with no gene, takes drags and clicks once a gene is open again
        const cleared = document.getElementById("plot-curve")._bound === false;
        _domain = null;
        _gene = { u: s1.universe, j: geneOrder(s1)[0] };
        renderActive();
        coexp.rebind = cleared && document.getElementById("plot-curve")._bound === true;
        _domain = null;
        state.tab = "metacells"; renderActive(); state.tab = "genes";
        _gene = null; _corrMin = 0.8; _corrMax = 1;
        __summary.coexp = coexp;
        state.sample = firstSample;

        // A pending level keeps its sample out of the selection, and the export table says so
        const ids = Object.keys(sel.samples);
        state.perSample[ids[1]].set.genes = false;
        const withPending = Object.keys(buildSelection().samples);
        state.tab = "export"; renderActive();
        const exportTable = document.getElementById("export-table").innerHTML;

        // Opening a level's tab with a sample selected applies that level, and only that one
        state.perSample[ids[1]].set = { cells: false, metacells: false, genes: false };
        state.sample = ids[1];
        showTab("genes");
        const visited = Object.assign({}, state.perSample[ids[1]].set);

        // "Apply settings to all samples" copies the selected sample's settings to every sample
        // and applies them; each keeps its own blacklists
        state.sample = ids[0];
        Object.assign(sstate(), { maxMito: 33 });
        gstate().minUmis = 7;
        // The pseudo-groups copy with the settings; the metacells and genes blacklisted stay the sample's own
        const own = (st) => [...st.blacklist].filter(id => !PSEUDO_GROUPS.includes(id));
        const black = { cells: own(state.perSample[ids[1]]), genes: state.perSample[ids[1]].genes.manual };
        applyToAll();
        const other = state.perSample[ids[1]];
        __summary.flow = {
            with_pending: withPending, export_pending: (exportTable.match(/>pending</g) || []).length, visited,
            applied: { maxMito: other.maxMito, minUmis: other.genes.minUmis, set: other.set,
                       kept_blacklist: JSON.stringify([own(other), other.genes.manual]) === JSON.stringify([black.cells, black.genes]) },
        };
        return sel;
    })()
`, context);

if (summaryPath) fs.writeFileSync(summaryPath, JSON.stringify(context.__summary));
process.stdout.write(JSON.stringify(selection, null, 1) + "\n");
