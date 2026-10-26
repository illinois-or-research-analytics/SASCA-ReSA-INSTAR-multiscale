#!/usr/bin/env python3
"""Mesoscale readout for SASCA-ReSA(-INSTAR) and CantataS runs.

  mesoscale_readout.py ikc    RUN_OUTPUT_DIR OUT_DIR --ikc-py PATH/IKC.py [--k 20] [--python PY] [--keep-edges]
  mesoscale_readout.py leiden RUN_OUTPUT_DIR OUT_DIR [--resolution 0.1] [--seed 1] [--keep-membership]

RUN_OUTPUT_DIR holds nodelist.csv (node_id, type, year, ...) and edgelist.csv (source,target; source cites target).

ikc     Iterative k-core clustering (IKC, Wedell et al.): cores are extracted recursively from the densest down, and a
        component is kept as a cluster only if it passes a positive-modularity test, so extracted nodes leave the graph
        and clusters are not k-cores. The IKC program itself is not part of this repository: pass its path (the version
        used here is chackoge/ERNIE_Plus, Illinois/clustering/eleanor/code/IKC.py, commit d3325da946; it needs networkit).
        Writes ikc_summary.csv, ikc_k<K>.csv (node, cluster, k_value, modularity) and ikc_stdout.txt.
leiden  Leiden CPM clustering of the undirected citation graph at one resolution, then per-cluster conductance
        (edges leaving the cluster / the smaller of the cluster's and the rest of the graph's edge volume; lower = better
        separated). Writes leiden_summary.csv and leiden_clusters.csv (cluster, size, cut, volume, conductance).

Both write OUT_DIR/*_summary.csv as measure,value rows. Conductance depends on cluster size (a perfect 8-paper clique in a
graph of mean degree 37 cannot go below about 0.8), so compare conductance between runs only alongside the share of papers
that sit in the clusters it describes (leiden_ge10_paper_share).
"""
import argparse, os, re, subprocess, sys, time
import numpy as np
import pandas as pd


def write_summary(path, S):
    pd.DataFrame(list(S.items()), columns=["measure", "value"]).to_csv(path, index=False)


def load(run_dir, nodelist_cols=("node_id", "type", "year")):
    nl = pd.read_csv(f"{run_dir}/nodelist.csv", usecols=list(nodelist_cols), engine="pyarrow")
    el = pd.read_csv(f"{run_dir}/edgelist.csv", engine="pyarrow")
    return nl, el


# ---------------------------------------------------------------- IKC
def cmd_ikc(a):
    os.makedirs(a.out, exist_ok=True)
    nl, el = load(a.run_dir)
    n = len(nl)
    tsv, out_csv, log = f"{a.out}/ikc_edges.tsv", f"{a.out}/ikc_k{a.k}.csv", f"{a.out}/ikc_stdout.txt"
    el.iloc[:, :2].to_csv(tsv, sep="\t", header=False, index=False)           # IKC wants tab-separated, no header
    t0 = time.time()
    with open(log, "w") as lf:
        r = subprocess.run([a.python or sys.executable, a.ikc_py, "-e", tsv, "-o", out_csv, "-k", str(a.k)],
                           stdout=lf, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        raise SystemExit(f"IKC failed with exit status {r.returncode}; see {log}")
    runtime = time.time() - t0
    if not a.keep_edges:
        os.remove(tsv)
    text = open(log).read()
    rej = lambda what: int((re.search(rf"not {what} :\s+(\d+)", text) or [0, 0])[1])
    c = pd.read_csv(out_csv, header=None, names=["node", "cluster", "k_value", "modularity"])
    cl = c[c.k_value >= a.k]                                                   # clusters with at least k connections per paper
    indeg = el.target.value_counts()
    S = {"ikc_k": a.k, "ikc_clusters": cl.cluster.nunique(), "ikc_clustered_papers": len(cl), "ikc_clustered_share": len(cl) / n,
         "ikc_rejected_modularity": rej("modular"), "ikc_rejected_kvalid": rej("k-valid"), "ikc_runtime_s": round(runtime, 1)}
    nan = float("nan")
    if len(cl):
        sizes = cl.groupby("cluster").size().sort_values(ascending=False)
        giant = cl[cl.cluster == sizes.index[0]]
        types = nl.set_index("node_id").loc[cl.node, "type"]
        S.update({"ikc_cluster_size_median": float(sizes.median()), "ikc_cluster_size_largest": int(sizes.iloc[0]),
                  "ikc_depth_max": int(cl.k_value.max()), "ikc_seed_papers_clustered": int((types == "seed").sum()),
                  "ikc_giant_modularity": float(giant.modularity.iloc[0]),
                  "ikc_giant_citation_share": float(indeg.reindex(giant.node).fillna(0).sum() / len(el)),
                  "ikc_giant_top100_inside": int(np.isin(indeg.nlargest(100).index, giant.node).sum()),
                  "ikc_giant_top1000_inside": int(np.isin(indeg.nlargest(1000).index, giant.node).sum())})
    else:
        S.update({k: nan for k in ("ikc_cluster_size_median", "ikc_cluster_size_largest", "ikc_depth_max", "ikc_seed_papers_clustered",
                                   "ikc_giant_modularity", "ikc_giant_citation_share", "ikc_giant_top100_inside", "ikc_giant_top1000_inside")})
    write_summary(f"{a.out}/ikc_summary.csv", S)
    print(pd.DataFrame(list(S.items()), columns=["measure", "value"]).to_string(index=False))


# ---------------------------------------------------------------- Leiden CPM + conductance
def cmd_leiden(a):
    import igraph as ig, leidenalg
    os.makedirs(a.out, exist_ok=True)
    nl, el = load(a.run_dir, nodelist_cols=("node_id",))
    el = el.iloc[:, :2]; el.columns = ["source", "target"]
    n, E = len(nl), len(el)
    g = ig.Graph.DataFrame(el, directed=False, use_vids=False)
    t0 = time.time()
    part = leidenalg.find_partition(g, leidenalg.CPMVertexPartition, resolution_parameter=a.resolution, seed=a.seed)
    runtime = time.time() - t0
    names = np.array(g.vs["name"]); mem = np.array(part.membership); k = mem.max() + 1; sizes = np.bincount(mem, minlength=k)
    if a.keep_membership:
        pd.DataFrame({"node_id": names, "cluster": mem}).to_csv(f"{a.out}/leiden_membership_{a.resolution}.csv", index=False)
    idx = pd.Index(names); u, v = idx.get_indexer(el.source), idx.get_indexer(el.target)
    cu, cv = mem[u], mem[v]; cross = cu != cv
    cut = np.bincount(cu[cross], minlength=k) + np.bincount(cv[cross], minlength=k)
    deg = np.bincount(u, minlength=len(names)) + np.bincount(v, minlength=len(names)); vol = np.bincount(mem, weights=deg, minlength=k)
    denom = np.minimum(vol, 2 * E - vol); cond = np.where(denom > 0, cut / np.maximum(denom, 1), np.nan)
    pd.DataFrame({"cluster": np.arange(k), "size": sizes, "cut": cut, "volume": vol, "conductance": cond}).to_csv(f"{a.out}/leiden_clusters.csv", index=False)
    S = {"leiden_resolution": a.resolution, "leiden_clusters": int(k), "leiden_singletons": int((sizes == 1).sum()), "leiden_runtime_s": round(runtime, 1)}
    pct = lambda x, p: float(np.percentile(x, p)) if len(x) else float("nan")
    for lo in (2, 10):
        sel = sizes >= lo; s, cd = sizes[sel], cond[sel]
        S.update({f"leiden_ge{lo}_clusters": int(sel.sum()), f"leiden_ge{lo}_paper_share": float(s.sum() / n),
                  f"leiden_ge{lo}_size_median": pct(s, 50), f"leiden_ge{lo}_size_p90": pct(s, 90), f"leiden_ge{lo}_size_p99": pct(s, 99),
                  f"leiden_ge{lo}_size_max": int(s.max()) if len(s) else 0, f"leiden_ge{lo}_median_conductance": pct(cd, 50),
                  f"leiden_ge{lo}_share_lt05_by_count": float((cd < 0.5).mean()) if len(cd) else float("nan")})
    big = sizes >= 10; low = big & (cond < 0.5)
    S["leiden_ge10_share_lt05_by_papers"] = float(sizes[low].sum() / sizes[big].sum()) if big.any() else float("nan")
    S["leiden_ge10_papers_in_lt05_clusters_share"] = float(sizes[low].sum() / n)
    write_summary(f"{a.out}/leiden_summary.csv", S)
    print(pd.DataFrame(list(S.items()), columns=["measure", "value"]).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    i = sub.add_parser("ikc"); i.add_argument("run_dir"); i.add_argument("out"); i.add_argument("--ikc-py", required=True)
    i.add_argument("--k", type=int, default=20); i.add_argument("--python"); i.add_argument("--keep-edges", action="store_true")
    l = sub.add_parser("leiden"); l.add_argument("run_dir"); l.add_argument("out"); l.add_argument("--resolution", type=float, default=0.1)
    l.add_argument("--seed", type=int, default=1); l.add_argument("--keep-membership", action="store_true")
    a = ap.parse_args()
    {"ikc": cmd_ikc, "leiden": cmd_leiden}[a.cmd](a)
