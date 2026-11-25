#!/usr/bin/env python3
"""Keep the outputs later analysis needs from being removed.

The protected set is configs/retained_outputs.tsv (see its header for why).

  check_workflow(workflow)   called from rules/99.retention.smk at every parse:
                             raises if a rule declares a retained output temp()
                             or its shell deletes one it does not produce, or if
                             the run is cut short (--until/--omit-from) without
                             --notemp.
  audit                      ledger check: every retained file that exists is
                             recorded in data/.retention_ledger.tsv; a recorded
                             file that has since disappeared is reported (and
                             makes the CLI exit 1).
  filter-deletions           stdin: NUL-separated paths about to be deleted;
                             stdout: the same minus anything protected (or
                             containing a protected file). Refusals go to stderr.
  retire GLOB --reason TEXT  drop ledger entries on purpose (e.g. results moved
                             to results_archive/), with the reason logged.

Run from the phase2 directory.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import fnmatch
import glob
import os
import re
import sys

LIST = "configs/retained_outputs.tsv"
LEDGER = "data/.retention_ledger.tsv"
RETIRED = "data/.retention_retired.tsv"


def load_patterns(path: str = LIST) -> list[tuple[str, str]]:
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            pat, _, why = line.partition("\t")
            out.append((pat.strip(), why.strip()))
    return out


def _wild_re(pattern: str) -> re.Pattern:
    """{family}/{sample}-style pattern -> regex. A literal '{x}' in a rule's
    output pattern also matches, so rule patterns and real paths both work."""
    parts = re.split(r"(\{[^}]+\})", pattern)
    rx = "".join("[^/]+" if p.startswith("{") else re.escape(p) for p in parts)
    return re.compile(rf"^{rx}$")


def _norm(p: str) -> str:
    p = os.path.normpath(str(p))
    return p[2:] if p.startswith("./") else p


# ---------------------------------------------------------------------------
# parse-time check
# ---------------------------------------------------------------------------
def check_workflow(workflow, list_path: str = LIST) -> None:
    from snakemake.io import is_flagged

    pats = load_patterns(list_path)
    rx = [(_wild_re(p), p) for p, _ in pats]
    # Basenames as regexes, for spotting deletions in shell text.
    base_rx = [(re.compile(_wild_re(os.path.basename(p)).pattern[1:-1]), p) for p, _ in pats]
    problems = []

    for rule in workflow.rules:
        own = []
        for f in rule.output:
            s = str(f)
            own.append(s)
            for r, p in rx:
                if r.match(_norm(s)) and is_flagged(f, "temp"):
                    problems.append(f"rule {rule.name}: retained output {p} is declared temp()")

        retained_inputs = set()
        for key, f in rule.input._allitems():
            if key and not callable(f) and isinstance(f, str):
                if any(r.match(_norm(f)) for r, _ in rx):
                    retained_inputs.add(key)

        shell = getattr(rule, "shellcmd", None) or ""
        for line in shell.splitlines():
            code = line.split("#", 1)[0]
            if not re.search(r"(^|[\s;&|(])(rm|unlink)\s", code) and "-delete" not in code:
                continue
            for key in retained_inputs:
                if re.search(r"\{input\." + re.escape(key) + r"\}(?![.\w])", code):
                    problems.append(f"rule {rule.name}: deletes retained input "
                                    f"{{input.{key}}}: {code.strip()}")
            for brx, p in base_rx:
                if brx.search(code) and not any(_wild_re(p).match(_norm(o)) for o in own):
                    problems.append(f"rule {rule.name}: shell deletes something named like "
                                    f"retained {p}: {code.strip()}")

    # A partial run reclaims temp() files as soon as the consumers INSIDE that
    # run are done, even if a consumer outside it still needs them. That is
    # how conifer's stage 1 (`--until run_virsorter2`) deleted every
    # sorted.bam before contig_coverage_stats could read it (2026-11). So a
    # run cut short with --until/--omit-from must keep temp files (--notemp);
    # the full run afterwards reclaims them once every consumer has finished.
    # Per-job executor processes (--target-jobs) never reclaim, so skip them.
    # Snakemake 8 fills workflow.dag_settings only after the Snakefile has
    # been parsed (it is None here), so --until/--omit-from are read from the
    # command line; storage_settings.notemp is already set at parse time.
    ss = getattr(workflow, "storage_settings", None)
    cut, argv = [], sys.argv[1:]
    flags = ("--until", "--omit-from")
    if "--target-jobs" not in argv:
        for i, a in enumerate(argv):
            name, eq, val = a.partition("=")
            if name in flags:
                vals = [val] if eq else []
                for b in argv[i + 1:] if not eq else []:
                    if b.startswith("-"):
                        break
                    vals.append(b)
                cut += vals
    notemp = bool(getattr(ss, "notemp", False)) or any(a in ("--notemp", "--nt") for a in argv)
    if cut and not notemp:
        problems.append(
            f"partial run (--until/--omit-from {' '.join(cut)}) without --notemp: "
            "temp() files would be reclaimed before consumers outside this run "
            "(e.g. contig_coverage_stats) read them. Add --notemp.")

    if problems:
        msg = ("Retention guard (configs/retained_outputs.tsv) -- refusing to run:\n  "
               + "\n  ".join(problems)
               + "\nThese files are read by later analysis. Make the consumer a declared "
                 "input instead of deleting, or change the list with a reason.")
        raise RuntimeError(msg)


# ---------------------------------------------------------------------------
# ledger
# ---------------------------------------------------------------------------
def _existing(pats) -> dict[str, int]:
    found = {}
    for p, _ in pats:
        g = re.sub(r"\{[^}]+\}", "*", p)
        for path in glob.glob(g):
            if _wild_re(p).match(path):
                try:
                    found[path] = os.path.getsize(path) if os.path.isfile(path) else -1
                except OSError:
                    pass
    return found


def _read_ledger() -> dict[str, tuple[str, str]]:
    led = {}
    if os.path.exists(LEDGER):
        with open(LEDGER) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                led[row["path"]] = (row["first_seen"], row["size"])
    return led


def audit(quiet: bool = False) -> int:
    pats = load_patterns()
    now = dt.datetime.now().isoformat(timespec="seconds")
    led = _read_ledger()
    found = _existing(pats)
    missing = sorted(p for p in led if p not in found)
    new = [p for p in found if p not in led]
    for p in new:
        led[p] = (now, str(found[p]))
    for p in missing:
        led.pop(p)  # reported below; keep a record of the loss instead
    os.makedirs(os.path.dirname(LEDGER), exist_ok=True)
    tmp = LEDGER + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["path", "first_seen", "size"])
        for p in sorted(led):
            w.writerow([p, *led[p]])
    os.replace(tmp, LEDGER)
    if missing:
        with open("data/.retention_losses.tsv", "a") as fh:
            for p in missing:
                fh.write(f"{now}\t{p}\n")
        print(f"[retention] ERROR: {len(missing)} retained output(s) recorded earlier are GONE "
              f"(logged in data/.retention_losses.tsv):", file=sys.stderr)
        for p in missing[:50]:
            print(f"  {p}", file=sys.stderr)
        if len(missing) > 50:
            print(f"  ... and {len(missing) - 50} more", file=sys.stderr)
        return 1
    if not quiet:
        print(f"[retention] ok: {len(led)} retained outputs tracked ({len(new)} new)", file=sys.stderr)
    return 0


def retire(pattern: str, reason: str) -> int:
    led = _read_ledger()
    gone = [p for p in led if fnmatch.fnmatch(p, pattern)]
    now = dt.datetime.now().isoformat(timespec="seconds")
    with open(RETIRED, "a") as fh:
        for p in gone:
            fh.write(f"{now}\t{p}\t{reason}\n")
            led.pop(p)
    with open(LEDGER, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["path", "first_seen", "size"])
        for p in sorted(led):
            w.writerow([p, *led[p]])
    print(f"[retention] retired {len(gone)} ledger entr(y/ies) matching {pattern}: {reason}", file=sys.stderr)
    return 0


def filter_deletions() -> int:
    pats = load_patterns()
    rx = [(_wild_re(p), p) for p, _ in pats]
    protected = None  # existing protected paths, globbed once, only if a dir is seen
    data = sys.stdin.buffer.read().split(b"\0")
    refused = 0
    out = sys.stdout.buffer
    for raw in data:
        if not raw:
            continue
        path = _norm(raw.decode())
        hit = next((p for r, p in rx if r.match(path)), None)
        if hit is None and os.path.isdir(path):
            # a directory that contains a protected file
            if protected is None:
                protected = sorted(_existing(pats))
            pre = path.rstrip("/") + "/"
            inside = next((x for x in protected if x.startswith(pre)), None)
            if inside:
                hit = inside
        if hit:
            refused += 1
            print(f"[retention] refusing to delete {path} (retained: {hit})", file=sys.stderr)
            continue
        out.write(raw + b"\0")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("audit"); a.add_argument("--quiet", action="store_true")
    sub.add_parser("filter-deletions")
    r = sub.add_parser("retire"); r.add_argument("glob"); r.add_argument("--reason", required=True)
    args = ap.parse_args()
    if args.cmd == "audit":
        return audit(args.quiet)
    if args.cmd == "filter-deletions":
        return filter_deletions()
    return retire(args.glob, args.reason)


if __name__ == "__main__":
    sys.exit(main())
