"""Smoke test: grow a small toy citation network, run the node-level readout on it and check the outputs.

    python tests/smoke_test.py

Builds ~8k papers (seeds plus 15 years of agents) in a temporary directory in the same file format the simulator
writes (output/nodelist.csv, output/edgelist.csv), runs `run` and `compare`, and asserts that every fidelity check
passes, the summary has the expected measures, and the permuted-trait control sits near 0.5.
"""
import os, subprocess, sys, tempfile
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); SCRIPT = os.path.join(HERE, "..", "microscale", "microscale_readout.py")


def make_toy_run(path, seed=0):
    rng = np.random.default_rng(seed)
    years, types = [], []
    for _ in range(300): years.append(int(rng.integers(2000, 2006))); types.append("seed")
    years = list(np.sort(years))
    src, tgt = [], []
    def cite(u, pool_mask):
        pool = np.flatnonzero(pool_mask)
        if len(pool) < 3: return
        k = min(int(rng.integers(3, 12)), len(pool)); w = indeg[pool] + 1.0; w /= w.sum()
        for t in rng.choice(pool, size=k, replace=False, p=w): src.append(u); tgt.append(int(t)); indeg[int(t)] += 1
    indeg = np.zeros(300 + 9000, dtype=np.int64)
    for u in range(300):                                    # seed-to-seed edges, older targets only
        cite(u, (np.array(years + [9999] * (len(indeg) - 300)) < years[u]))
    n = 300
    for y in range(2006, 2021):
        n_new = int(0.25 * n) + 50; ys = np.array(years)
        for _ in range(n_new):
            if n >= len(indeg): break
            years.append(y); types.append("agent"); u = n; n += 1
            cite(u, np.concatenate([ys[:u] < y, np.zeros(len(indeg) - u, bool)]))
        ys = np.array(years)
    N = len(years); outdeg = np.bincount(src, minlength=N); indeg = np.bincount(tgt, minlength=N)
    agent = np.array(types) == "agent"
    nl = pd.DataFrame(dict(
        node_id=np.arange(N), type=types, year=years, alpha=np.where(agent, rng.random(N), -1.0),
        pa_weight=np.where(agent, rng.random(N), -1.0), fit_weight=np.where(agent, rng.random(N), -1.0),
        num_authors_weight=np.where(agent, rng.random(N), -1.0), author_reputation_weight=np.where(agent, rng.random(N), -1.0),
        fit_peak_value=np.where(agent, np.minimum(1000, (1 / rng.random(N) ** 0.6).astype(int)), 1),
        in_degree=indeg, out_degree=outdeg, assigned_out_degree=np.where(agent, outdeg, -1),
        num_authors=np.maximum(1, (rng.pareto(2.0, N) + 1).astype(int)), initial_author_reputation=rng.integers(0, 20, N),
        final_author_reputation=rng.integers(0, 40, N), author_id=rng.integers(0, N // 3, N)))
    os.makedirs(os.path.join(path, "output"), exist_ok=True)
    nl.to_csv(os.path.join(path, "output", "nodelist.csv"), index=False)
    pd.DataFrame(dict(source=src, target=tgt)).to_csv(os.path.join(path, "output", "edgelist.csv"), index=False)
    return N, len(src)


def main():
    with tempfile.TemporaryDirectory() as tmp:
        outs = []
        for i in range(2):
            run = os.path.join(tmp, f"toy{i}"); N, E = make_toy_run(run, seed=i)
            out = os.path.join(tmp, f"res{i}"); outs.append(out)
            r = subprocess.run([sys.executable, SCRIPT, "run", os.path.join(run, "output"), out, "--keep-vectors"], capture_output=True, text=True)
            assert r.returncode == 0, r.stderr[-2000:]
            checks = open(os.path.join(out, "checks.txt")).read()
            assert "FAIL" not in checks, checks
        S = pd.read_csv(os.path.join(outs[0], "summary.csv")).set_index("measure").value
        for m in ["indeg_all_median", "indeg_all_p99", "indeg_all_max", "never_cited_all", "uncited_at_5y", "peak_age_median", "tail_gamma",
                  "beta_eff_age_6-10", "pair_prob_fitness", "variation_unexplained_by_traits", "top1pct_of_cited_share"]:
            assert m in S.index, f"missing {m}"
        assert abs(S["pair_prob_permuted_control"] - 0.5) < 0.02, S["pair_prob_permuted_control"]
        assert S["checks_failed"] == 0
        assert os.path.exists(os.path.join(outs[0], "node_vector.parquet"))
        r = subprocess.run([sys.executable, SCRIPT, "compare", f"toy=" + ",".join(outs)], capture_output=True, text=True)
        assert r.returncode == 0 and "| indeg_all_max |" in r.stdout, r.stderr[-2000:]
        print(f"toy runs: {N} nodes, {E} edges"); print("smoke test passed")


if __name__ == "__main__":
    main()
