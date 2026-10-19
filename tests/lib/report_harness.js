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
 * CHOICES: { "samples": { "<id>": { maxMito, maxDoublet, maxAlpha, includeOutliers,
 *            includeExcluded, blacklist: [], genes: { ... } } },
 *            "genes": { minUmis, pfam: [], manual, excludeMito } }
 *
 * The top-level "genes" go to every sample named, under its own "genes". Every sample
 * named counts as set at both levels; the others stay pending.
 *
 * SUMMARY_OUT also receives "flow": what the selection holds with a level pending, what
 * "Apply settings to all samples" copies, and whether a version-1 selection reads back the
 * same.
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
            for (const k of ["maxMito", "maxDoublet", "maxAlpha", "includeOutliers", "includeExcluded"]) {
                if (k in s) st[k] = s[k];
            }
            if (s.blacklist) st.blacklist = new Set(s.blacklist);
            const g = Object.assign({}, c.genes || {}, s.genes || {});
            if ("minUmis" in g) st.genes.minUmis = g.minUmis;
            if (g.pfam) st.genes.pfam = new Set(g.pfam);
            if ("manual" in g) st.genes.manual = g.manual;
            if ("excludeMito" in g) st.genes.excludeMito = g.excludeMito;
            st.set = { cells: true, genes: true };
        }

        // The summary panel and every tab render without throwing on these choices
        renderSummary();
        for (const tab of ["cells", "genes", "export"]) { state.tab = tab; renderActive(); }
        __summary = {
            settings: document.getElementById("settings-cells").innerHTML +
                      document.getElementById("settings-genes").innerHTML,
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

        // A version-1 selection (one set of gene rules for the run, the metacell blacklist beside
        // the cell rules) restores the same choices, as long as every sample shared the gene rules
        const first = Object.values(sel.samples)[0];
        const v1 = { schema: sel.schema, schema_version: 1, samples: {},
                     gene_rules: Object.assign({}, first.gene_rules, { manual: first.gene_rules.blacklist, blacklist: undefined }) };
        for (const [id, e] of Object.entries(sel.samples)) {
            const cr = Object.assign({}, e.cell_rules); delete cr.blacklist;
            v1.samples[id] = { fingerprint: e.fingerprint, cell_rules: cr, blacklist: e.cell_rules.blacklist };
        }
        const fromV1 = loadSelection(JSON.parse(JSON.stringify(v1)));
        const afterV1 = buildSelection();
        afterV1.created = sel.created;

        // The plots are the threshold handles: replay a drag of the mito line, a drag of the
        // rank plot's line (log axis, data units) and a click on a gene-histogram bar. Done
        // after the selection is built, so they do not change it.
        const fire = (id, event, payload) => {
            const h = document.getElementById(id)._handlers[event];
            if (!h) throw new Error("no " + event + " handler on " + id);
            h(payload);
        };
        state.tab = "cells"; renderActive();
        fire("plot-mito", "plotly_relayout", { "shapes[0].x0": 12.34, "shapes[0].x1": 12.38, "shapes[0].y0": 0.2 });
        const mito = sstate().maxMito;
        state.tab = "genes"; renderActive();
        fire("plot-rank", "plotly_relayout", { "shapes[0].y0": 50, "shapes[0].y1": 50 });
        const rank = gstate().minUmis;
        fire("plot-genehist", "plotly_click", { points: [{ x: 2 }] });
        __summary.drag = { mito, rank, click: gstate().minUmis };

        // A pending level keeps its sample out of the selection, and the export table says so
        const ids = Object.keys(sel.samples);
        state.perSample[ids[1]].set.genes = false;
        const withPending = Object.keys(buildSelection().samples);
        state.tab = "export"; renderActive();
        const exportTable = document.getElementById("export-table").innerHTML;

        // Opening a level's tab with a sample selected applies that level, and only that one
        state.perSample[ids[1]].set = { cells: false, genes: false };
        state.sample = ids[1];
        showTab("genes");
        const visited = Object.assign({}, state.perSample[ids[1]].set);

        // "Apply settings to all samples" copies the selected sample's settings to every sample
        // and applies them; each keeps its own blacklists
        state.sample = ids[0];
        Object.assign(sstate(), { maxMito: 33 });
        gstate().minUmis = 7;
        const black = { cells: [...state.perSample[ids[1]].blacklist], genes: state.perSample[ids[1]].genes.manual };
        applyToAll();
        const other = state.perSample[ids[1]];
        __summary.flow = {
            v1_problems: fromV1.problems, v1_same: JSON.stringify(afterV1) === JSON.stringify(sel),
            with_pending: withPending, export_pending: (exportTable.match(/>pending</g) || []).length, visited,
            applied: { maxMito: other.maxMito, minUmis: other.genes.minUmis, set: other.set,
                       kept_blacklist: JSON.stringify([[...other.blacklist], other.genes.manual]) === JSON.stringify([black.cells, black.genes]) },
        };
        return sel;
    })()
`, context);

if (summaryPath) fs.writeFileSync(summaryPath, JSON.stringify(context.__summary));
process.stdout.write(JSON.stringify(selection, null, 1) + "\n");
