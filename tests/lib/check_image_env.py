#!/usr/bin/env python3
"""Check that a container image provides every package its module's environment.yml lists.

A module runs from its conda environment under the conda profiles and from its container
under the container profiles, so the two must provide the same tools. They drift apart
silently: an environment.yml gains a package, or a pin moves, and nobody rebuilds the image.
The task then fails only on the container profiles, or worse, runs with another version.

The installed packages are read from the image beforehand (see containers.sh) as one entry
per line, either a conda-meta file name (``name-version-build.json``) or a Python
distribution directory (``name-version.dist-info``). Every dependency in the environment
file, conda and pip alike, must be installed and satisfy its version constraint. Extra
packages in the image are fine.

Usage: check_image_env.py ENVIRONMENT_YML INSTALLED_LIST

Standard library only, since PyYAML is not guaranteed on the host: environment files use a
small, regular subset of YAML that is parsed line by line here.
"""

import re
import sys


def normalise(name):
    """Package names compare case-insensitively, with '_' and '.' equivalent to '-' (PEP 503)."""
    return re.sub(r"[-_.]+", "-", name).lower()


def read_dependencies(path):
    """Return (kind, spec) for every dependency in an environment file; kind is conda or pip."""
    deps = []
    in_deps = False
    pip_indent = None
    with open(path) as fh:
        for raw in fh:
            line = raw.split("#", 1)[0].rstrip()
            if not line.strip():
                continue
            indent = len(line) - len(line.lstrip())
            text = line.strip()
            if indent == 0:
                in_deps = text == "dependencies:"
                pip_indent = None
                continue
            if not in_deps or not text.startswith("- "):
                continue
            item = text[2:].strip().strip("'\"")
            if pip_indent is not None and indent > pip_indent:
                deps.append(("pip", item))
                continue
            pip_indent = None
            if item == "pip:":
                pip_indent = indent
            else:
                deps.append(("conda", item))
    return deps


def parse_spec(kind, spec):
    """Split a dependency into its name and a list of (operator, version) constraints."""
    spec = spec.split("::", 1)[-1]                 # drop a channel prefix
    match = re.match(r"^\s*([A-Za-z0-9_.\-]+)(\[[^\]]*\])?\s*(.*)$", spec)
    name, rest = match.group(1), match.group(3).replace(" ", "")
    constraints = []
    if kind == "conda" and re.match(r"^=[^=<>!]", rest):
        # name=version or name=version=build: a version prefix, the build is not compared
        version = rest[1:].split("=", 1)[0]
        constraints.append(("=", version))
    elif rest:
        for part in rest.split(","):
            op, version = re.match(r"^(==|>=|<=|!=|~=|>|<|=)?(.+)$", part).groups()
            constraints.append((op or "==", version))
    return name, constraints


def version_key(version):
    """Comparable key for a version string, numeric where the component is numeric."""
    key = []
    for component in re.split(r"[.\-+_]", version):
        number = re.match(r"^(\d+)(.*)$", component)
        key.append((int(number.group(1)), number.group(2)) if number else (-1, component))
    return key


def satisfies(installed, op, wanted):
    wanted = re.sub(r"\.\*$", "", wanted)
    have, want = version_key(installed), version_key(wanted)
    if op == "=":                                  # conda: 3.11 matches 3.11.9
        return have[:len(want)] == want
    if op == "==":
        return have[:len(want)] == want and all(c == (0, "") for c in have[len(want):])
    padded = have + [(0, "")] * max(0, len(want) - len(have))
    if op == ">=":
        return padded >= want
    if op == "<=":
        return padded <= want
    if op == ">":
        return padded > want
    if op == "<":
        return padded < want
    if op == "!=":
        return have[:len(want)] != want
    if op == "~=":                                 # compatible release: >= v, == v[:-1].*
        return padded >= want and have[:len(want) - 1] == want[:-1]
    raise ValueError(f"unsupported operator {op}")


def read_installed(path):
    """Map normalised name -> version for conda packages and for Python distributions."""
    conda, python = {}, {}
    with open(path) as fh:
        for line in fh:
            entry = line.strip().rsplit("/", 1)[-1]
            if entry.endswith(".json"):
                parts = entry[:-len(".json")].rsplit("-", 2)
                if len(parts) == 3:
                    conda[normalise(parts[0])] = parts[1]
            elif entry.endswith(".dist-info"):
                parts = entry[:-len(".dist-info")].rsplit("-", 1)
                if len(parts) == 2:
                    python[normalise(parts[0])] = parts[1]
    return conda, python


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: check_image_env.py ENVIRONMENT_YML INSTALLED_LIST")
    deps = read_dependencies(sys.argv[1])
    conda, python = read_installed(sys.argv[2])
    if not conda:
        sys.exit("No conda packages found in the image listing; is this a conda-based image?")

    problems = []
    for kind, spec in deps:
        name, constraints = parse_spec(kind, spec)
        key = normalise(name)
        # pip names a package as on PyPI and conda as on its channel, so look in both
        sources = [("pip", python), ("conda", conda)]
        if kind == "conda":
            sources.reverse()
        found = next(((label, src[key]) for label, src in sources if key in src), None)
        if found is None:
            problems.append(f"missing     {spec}")
            continue
        label, version = found
        failed = [f"{op}{v}" for op, v in constraints if not satisfies(version, op, v)]
        if failed:
            # A conda dependency found only among the Python packages may be a PyPI namesake
            origin = f" from {label}" if label != kind else ""
            problems.append(f"mismatch    {spec}  (image has {version}{origin})")
        else:
            print(f"ok          {spec}  ({version})")

    for problem in problems:
        print(problem)
    print(f"\n{len(deps) - len(problems)}/{len(deps)} dependencies satisfied")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
