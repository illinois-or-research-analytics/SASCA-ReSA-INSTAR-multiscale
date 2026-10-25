#!/usr/bin/env python3
"""Macroscale readout: concentration of citations and triangle counts.

  macroscale_readout.py run RUN_OUTPUT_DIR OUT_DIR [--triangles-bin PATH/triangles | --skip-triangles]

RUN_OUTPUT_DIR holds nodelist.csv and edgelist.csv (source,target; source cites target). Writes OUT_DIR/macro_summary.csv.

Concentration (over all papers, uncited ones included): share of all citations received by the 10 most-cited papers
(top10_share) and by the most-cited 1% of cited papers (top1pct_of_cited_share), and the Gini coefficient of in-degree
(0 = equal, 1 = one paper has everything).
Triangles come from the C++ program in this folder (build: `make -C macroscale`), which counts undirected triangles
(directions ignored, duplicate edges merged) and reports transitivity, mean local clustering (over papers with at least
two connections) and directed out-/in-closure. It was validated against an independent sparse-matrix count on a real
50k-paper slice (triangles, connected triples and transitivity match exactly). Transitivity is dominated by open triples
around hubs, so triangles per paper and mean local clustering are the steadier comparisons.
"""
import argparse, os, subprocess
import numpy as np
import pandas as pd


def concentration(run_dir):
    el = pd.read_csv(f"{run_dir}/edgelist.csv", usecols=["target"], engine="pyarrow")
    n = len(pd.read_csv(f"{run_dir}/nodelist.csv", usecols=["node_id"], engine="pyarrow"))
    cnt = np.sort(el.target.value_counts().to_numpy())[::-1]                 # in-degree of every cited paper, largest first
    E = int(cnt.sum())
    s = np.concatenate([np.zeros(n - len(cnt), dtype=np.int64), cnt[::-1]])  # ascending, with the never-cited papers as zeros
    gini = float((2 * np.arange(1, n + 1) - n - 1).dot(s) / (n * E))
    return {"macro_papers": n, "macro_citations": E, "top10_share": float(cnt[:10].sum() / E),
            "top1pct_of_cited_share": float(cnt[:len(cnt) // 100].sum() / E), "gini_indegree": gini}


def triangles(run_dir, binary, out):
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("SLURM_CPUS_PER_TASK", os.environ.get("OMP_NUM_THREADS", "4")))
    txt = f"{out}/triangles_output.txt"
    with open(txt, "w") as f:
        r = subprocess.run([binary, f"{run_dir}/edgelist.csv"], stdout=f, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0:
        raise SystemExit(f"triangles failed with exit status {r.returncode}; see {txt}")
    t = dict(l.strip().split("=", 1) for l in open(txt) if "=" in l)
    return {"triangles": float(t["triangles"]), "triangles_transitive": float(t["triangles_transitive"]),
            "triangles_cyclic": float(t["triangles_cyclic"]), "triangles_with_reciprocal_arc": float(t["triangles_with_reciprocal_arc"]),
            "global_transitivity": float(t["global_transitivity"]), "mean_local_clustering": float(t["mean_local_clustering"]),
            "out_closure_global": float(t["out_closure_global"]), "in_closure_global": float(t["in_closure_global"])}


def cmd_run(a):
    os.makedirs(a.out, exist_ok=True)
    S = concentration(a.run_dir)
    if not a.skip_triangles:
        if not a.triangles_bin or not os.path.exists(a.triangles_bin):
            raise SystemExit("triangles program not found: build it with `make -C macroscale` and pass --triangles-bin, or use --skip-triangles")
        S.update(triangles(a.run_dir, a.triangles_bin, a.out))
        S["triangles_per_paper"] = S["triangles"] / S["macro_papers"]
    pd.DataFrame(list(S.items()), columns=["measure", "value"]).to_csv(f"{a.out}/macro_summary.csv", index=False)
    print(pd.DataFrame(list(S.items()), columns=["measure", "value"]).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run"); r.add_argument("run_dir"); r.add_argument("out")
    r.add_argument("--triangles-bin"); r.add_argument("--skip-triangles", action="store_true")
    a = ap.parse_args()
    cmd_run(a)
