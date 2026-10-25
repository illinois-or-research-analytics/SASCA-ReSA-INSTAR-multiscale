"""Smoke test of the mesoscale and macroscale stages and of merge/compare, on the toy network from smoke_test.py.

    python tests/smoke_test_pipeline.py

Runs the micro stage, the Leiden stage (needs igraph and leidenalg), the macro stage (concentration always; triangles only if
the program has been built with `make -C macroscale`), then `multiscale.py merge` and `compare`. The IKC stage is not tested
here because it needs the external IKC.py and networkit; test it on a real run with slurm/run_pipeline.sh.
"""
import importlib.util, os, subprocess, sys, tempfile
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
sys.path.insert(0, HERE)
from smoke_test import make_toy_run


def run(*args):
    r = subprocess.run([sys.executable, *args], capture_output=True, text=True)
    assert r.returncode == 0, (args, r.stderr[-2000:])
    return r.stdout


def main():
    have = lambda m: importlib.util.find_spec(m) is not None
    tri = os.path.join(ROOT, "macroscale", "triangles")
    with tempfile.TemporaryDirectory() as tmp:
        outs = []
        for i in range(2):
            runs_dir = os.path.join(tmp, "runs"); name = f"toy{i}"
            make_toy_run(os.path.join(runs_dir, name), seed=i)
            inp, out = os.path.join(runs_dir, name, "output"), os.path.join(tmp, "out", name); outs.append(out)
            run(os.path.join(ROOT, "microscale", "microscale_readout.py"), "run", inp, os.path.join(out, "micro"))
            if have("igraph") and have("leidenalg"):
                run(os.path.join(ROOT, "mesoscale", "mesoscale_readout.py"), "leiden", inp, os.path.join(out, "meso"), "--resolution", "0.1")
            use_tri = os.path.exists(tri)
            run(os.path.join(ROOT, "macroscale", "macroscale_readout.py"), "run", inp, os.path.join(out, "macro"),
                *(["--triangles-bin", tri] if use_tri else ["--skip-triangles"]))
            run(os.path.join(ROOT, "multiscale.py"), "merge", out)
        S = pd.read_csv(os.path.join(outs[0], "summary_all.csv")).set_index("measure").value
        need = ["micro.indeg_all_median", "macro.top10_share", "macro.gini_indegree"]
        if have("igraph") and have("leidenalg"): need += ["meso.leiden_clusters", "meso.leiden_ge10_median_conductance"]
        if os.path.exists(tri): need += ["macro.triangles", "macro.global_transitivity"]
        for m in need: assert m in S.index, f"missing {m}"
        assert 0 <= S["macro.gini_indegree"] <= 1 and 0 < S["macro.top10_share"] <= 1
        table = run(os.path.join(ROOT, "multiscale.py"), "compare", "toy=" + ",".join(outs))
        assert "| micro.indeg_all_median |" in table and "| macro.top10_share |" in table
        print(f"merged {len(S)} measures; stages: igraph/leidenalg={'yes' if have('igraph') and have('leidenalg') else 'skipped'}, "
              f"triangles={'yes' if os.path.exists(tri) else 'skipped (not built)'}; IKC not tested here")
        print("pipeline smoke test passed")


if __name__ == "__main__":
    main()
