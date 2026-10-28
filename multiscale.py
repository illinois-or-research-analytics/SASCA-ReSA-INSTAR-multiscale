#!/usr/bin/env python3
"""Combine and compare the micro-, meso- and macroscale summaries of runs.

  multiscale.py merge RUN_OUT_DIR                       # -> RUN_OUT_DIR/summary_all.csv, RUN_OUT_DIR/stages.txt
  multiscale.py compare LABEL=DIR,DIR,... [...] [--scales micro,meso,macro] [--out FILE.csv]

RUN_OUT_DIR is one run's output of the pipeline (slurm/run_pipeline.sh), with the stage folders
  micro/summary.csv   meso/ikc_summary.csv   meso/leiden_summary.csv   macro/macro_summary.csv
`merge` puts every measure into summary_all.csv as "<scale>.<measure>" (measure,value) and records which stages were found
in stages.txt, so a run whose Leiden stage was still running or failed still merges (its meso.leiden_* measures are absent).
`compare` prints a markdown table of measures by group; replicates given as LABEL=DIR,DIR,... show as "min to max".
`report` prints a readable markdown summary of a curated set of measures, grouped by scale, with plain-English names,
percentages and a "what it means" column (multiscale.py report LABEL=DIR,DIR,... [...] [--out FILE.md]). A star marks the
measures that answer each scale's central question, a tick marks integrity checks that must hold, and when every run has
two or more replicates a column flags the measures whose runs differ by more than their replicate spread.
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


# (measure, plain-English name, kind, role, what it means). kind: n = count, f = number, p = share shown as a percent.
# role: "key" = answers the scale's central question, "check" = must hold in every valid run, "" = context for reading the rest.
REPORT = [
    ("Run size and integrity", "", [
        ("micro.n_nodes", "Papers", "n", "", "Seed papers plus simulated ones."),
        ("micro.n_agents", "Simulated papers", "n", "", "Papers added by the simulation (agents)."),
        ("micro.n_edges", "Citations", "n", "", "Reference links between papers."),
        ("micro.checks_failed", "Integrity checks failed", "n", "check", "Structural checks on the output; must be 0."),
    ]),
    ("Micro scale: individual papers", "Central question: how are citations spread over papers, and what predicts who gets them?", [
        ("micro.indeg_all_median", "Citations received, median paper", "n", "", "The typical paper."),
        ("micro.indeg_all_p99", "Citations received, top 1% paper", "n", "", "The threshold for the top 1% of papers."),
        ("micro.indeg_all_max", "Citations received, most-cited paper", "n", "", "The extreme of the tail."),
        ("micro.never_cited_all", "Never cited", "p", "key", "Share of papers nobody cites."),
        ("micro.uncited_at_2y", "Uncited after 2 years", "p", "", "Share of papers with no citation two years after publication."),
        ("micro.uncited_at_10y", "Uncited after 10 years", "p", "", "Same, after ten years."),
        ("micro.still_cited_final_year", "Old papers still being cited", "p", "", "Among papers aged 10 or more, the share that receive a citation in the final simulated year."),
        ("micro.tail_gamma", "Tail exponent (gamma)", "f", "key", "Lower means a heavier tail of very highly cited papers."),
        ("micro.beta_eff_age_6-10", "Attachment exponent, papers aged 6-10", "f", "key", "1 = chance of being cited grows in proportion to citations already held; above 1 = rich get richer faster."),
        ("micro.pair_prob_fitness", "Higher fitness wins a head-to-head", "f", "", "0.5 = fitness does not matter; the permuted control sits at 0.5."),
        ("micro.pair_prob_author_reputation", "Higher author reputation wins", "f", "", "Same test for author reputation."),
        ("micro.r2_all_traits", "Citation variation explained by traits", "p", "key", "Fitness, author count and author reputation together."),
        ("micro.variation_unexplained_by_traits", "Citation variation left to luck", "p", "", "What the traits do not explain."),
    ]),
    ("Meso scale: communities", "Central question: do papers fall into groups that cite mostly within themselves?", [
        ("meso.leiden_ge10_clusters", "Leiden clusters of 10+ papers", "n", "", "Communities found by Leiden at the stated resolution."),
        ("meso.leiden_ge10_paper_share", "Papers in those clusters", "p", "", "Share of all papers sitting in a cluster of 10 or more."),
        ("meso.leiden_ge10_median_conductance", "Median conductance of those clusters", "f", "key", "0 = cites only inside the cluster, 1 = fully mixed with the rest. Small clusters cannot score low."),
        ("meso.leiden_ge10_share_lt05_by_count", "Clusters that are well separated", "p", "", "Share of clusters with conductance below 0.5."),
        ("meso.ikc_clusters", "Dense cores found by IKC", "n", "", "Clusters surviving the modularity test at k=20."),
        ("meso.ikc_clustered_share", "Papers in IKC cores", "p", "key", "Share of all papers inside such a core; independent of the Leiden labels."),
        ("meso.ikc_cluster_size_largest", "Largest IKC core", "n", "", "Papers in the biggest core."),
    ]),
    ("Macro scale: the whole network", "Central question: how concentrated is attention overall, and is the structure sound?", [
        ("macro.top10_share", "Top 10 papers' share of citations", "p", "key", "How much of all citation mass the ten most-cited papers hold."),
        ("macro.top1pct_of_cited_share", "Top 1% of cited papers' share", "p", "", "Share of citations going to the top 1% of papers that were cited."),
        ("macro.gini_indegree", "Inequality of citations (Gini)", "f", "key", "0 = equal, 1 = one paper gets everything."),
        ("macro.global_transitivity", "Global transitivity", "f", "", "Chance that two papers citing the same paper cite each other; hubs push it down."),
        ("macro.mean_local_clustering", "Mean local clustering", "f", "", "Average closure around each paper."),
        ("macro.triangles_per_paper", "Citation triangles per paper", "f", "", "Sets of three papers all linked to one another by citations, per paper."),
        ("macro.triangles_cyclic", "Cyclic triangles", "n", "check", "Must be 0: citations only point backward in time."),
    ]),
]
MARK = {"key": "★", "check": "✔", "": ""}


def cmd_report(a):
    groups = {}
    for spec in a.dirs:
        label, paths = spec.split("=", 1) if "=" in spec else (os.path.basename(os.path.normpath(spec)), spec)
        groups[label] = pd.concat([pd.read_csv(os.path.join(d, "summary_all.csv")).set_index("measure").value.rename(i)
                                   for i, d in enumerate(paths.split(","))], axis=1)
    spread = all(S.shape[1] >= 2 for S in groups.values()) and len(groups) > 1

    def values(S, m):
        return pd.to_numeric(S.loc[m], errors="coerce").dropna().astype(float) if m in S.index else pd.Series(dtype=float)

    def cell(S, m, kind):
        v = values(S, m)
        if not len(v): return "-"
        f = {"n": lambda x: f"{x:,.0f}", "p": lambda x: f"{100 * x:.1f}%", "f": lambda x: f"{x:.3g}"}[kind]
        lo, hi = f(v.min()), f(v.max())
        return lo if lo == hi else f"{lo} to {hi}"

    def separated(m):                 # every pair of runs has non-overlapping replicate ranges
        r = [(v.min(), v.max()) for v in (values(S, m) for S in groups.values()) if len(v)]
        return len(r) == len(groups) and all(a[1] < b[0] or b[1] < a[0] for i, a in enumerate(r) for b in r[i + 1:])

    reps = ", ".join(f"{l} ({S.shape[1]} replicate{'s' if S.shape[1] > 1 else ''})" for l, S in groups.items())
    out = ["# Multiscale readout", "", f"Runs: {reps}. Replicates are shown as a range. Shares are percentages of papers unless stated.", "",
           "**★** answers the scale's central question. **✔** is an integrity check that must hold."
           + (" **Differs** marks runs separated by more than their replicate spread." if spread else
              " With one replicate per run there is no spread to judge differences against, so no row is flagged as a real difference."), ""]
    for title, intro, rows in REPORT:
        head = ["", "Measure"] + list(groups) + (["Differs"] if spread else []) + ["What it means"]
        out += [f"## {title}", "", intro, ""] if intro else [f"## {title}", ""]
        out += ["| " + " | ".join(head) + " |", "|:-:|---|" + "---:|" * len(groups) + (":-:|" if spread else "") + "---|"]
        for m, name, k, role, gloss in rows:
            out.append(f"| {MARK[role]} | {name} | " + " | ".join(cell(S, m, k) for S in groups.values())
                       + (f" | {'yes' if separated(m) else ''}" if spread else "") + f" | {gloss} |")
        out.append("")
    text = "\n".join(out)
    if a.out: open(a.out, "w").write(text)
    print(text)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("merge"); m.add_argument("run_out")
    c = sub.add_parser("compare"); c.add_argument("dirs", nargs="+"); c.add_argument("--scales", default="micro,meso,macro"); c.add_argument("--out")
    r = sub.add_parser("report"); r.add_argument("dirs", nargs="+"); r.add_argument("--out")
    a = ap.parse_args()
    {"merge": cmd_merge, "compare": cmd_compare, "report": cmd_report}[a.cmd](a)
