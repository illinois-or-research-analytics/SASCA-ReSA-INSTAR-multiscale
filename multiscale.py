#!/usr/bin/env python3
"""Combine and compare the micro-, meso- and macroscale summaries of runs.

  multiscale.py merge RUN_OUT_DIR                       # -> RUN_OUT_DIR/summary_all.csv, RUN_OUT_DIR/stages.txt
  multiscale.py compare LABEL=DIR,DIR,... [...] [--scales micro,meso,macro] [--out FILE.csv]

RUN_OUT_DIR is one run's output of the pipeline (slurm/run_pipeline.sh), with the stage folders
  micro/summary.csv   meso/ikc_summary.csv   meso/leiden_summary.csv   macro/macro_summary.csv
`merge` puts every measure into summary_all.csv as "<scale>.<measure>" (measure,value) and records which stages were found
in stages.txt, so a run whose Leiden stage was still running or failed still merges (its meso.leiden_* measures are absent).
`compare` prints a markdown table of measures by group; replicates given as LABEL=DIR,DIR,... show as "min to max".
"""
import argparse, os
import numpy as np
import pandas as pd

STAGES = {"micro": ["micro/summary.csv"], "meso": ["meso/ikc_summary.csv", "meso/leiden_summary.csv"], "macro": ["macro/macro_summary.csv"]}


def cmd_merge(a):
    rows, status = [], []
    for scale, files in STAGES.items():
        for f in files:
            p = os.path.join(a.run_out, f)
            ok = os.path.exists(p)
            status.append(f"{'present' if ok else 'MISSING'}  {f}")
            if ok:
                d = pd.read_csv(p)
                rows += [(f"{scale}.{m}", v) for m, v in zip(d.measure, d.value)]
    pd.DataFrame(rows, columns=["measure", "value"]).to_csv(os.path.join(a.run_out, "summary_all.csv"), index=False)
    open(os.path.join(a.run_out, "stages.txt"), "w").write("\n".join(status) + "\n")
    print("\n".join(status)); print(f"wrote {len(rows)} measures to {os.path.join(a.run_out, 'summary_all.csv')}")


def fmt(x):
    if pd.isna(x): return ""
    if x == int(x) and abs(x) >= 100: return f"{int(x):,}"
    return f"{x:.4g}"


def cmd_compare(a):
    scales = tuple(s + "." for s in a.scales.split(","))
    groups = {}
    for spec in a.dirs:
        label, paths = spec.split("=", 1) if "=" in spec else (os.path.basename(os.path.normpath(spec)), spec)
        groups[label] = pd.concat([pd.read_csv(os.path.join(d, "summary_all.csv")).set_index("measure").value.rename(i)
                                   for i, d in enumerate(paths.split(","))], axis=1)
    measures = []
    for S in groups.values():
        measures += [m for m in S.index if m.startswith(scales) and m not in measures]
    first_seen = {m: i for i, m in enumerate(measures)}                      # keep each scale's measures in their original order
    order = a.scales.split(",")
    measures.sort(key=lambda m: (order.index(m.split(".", 1)[0]), first_seen[m]))
    T = {}
    for label, S in groups.items():
        cells = []
        for m in measures:
            vals = pd.to_numeric(S.loc[m], errors="coerce").dropna().astype(float) if m in S.index else pd.Series(dtype=float)
            lo, hi = (vals.min(), vals.max()) if len(vals) else (np.nan, np.nan)
            cells.append(fmt(lo) if fmt(lo) == fmt(hi) else f"{fmt(lo)} to {fmt(hi)}")
        T[label] = cells
    T = pd.DataFrame(T, index=measures)
    if a.out: T.to_csv(a.out)
    print("| measure | " + " | ".join(T.columns) + " |"); print("|" + "---|" * (len(T.columns) + 1))
    for m, row in T.iterrows(): print(f"| {m} | " + " | ".join(row) + " |")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("merge"); m.add_argument("run_out")
    c = sub.add_parser("compare"); c.add_argument("dirs", nargs="+"); c.add_argument("--scales", default="micro,meso,macro"); c.add_argument("--out")
    a = ap.parse_args()
    {"merge": cmd_merge, "compare": cmd_compare}[a.cmd](a)
