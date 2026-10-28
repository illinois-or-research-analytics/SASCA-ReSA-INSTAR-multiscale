#!/usr/bin/env python3
"""Microscale (paper-level) readout for SASCA-ReSA(-INSTAR) runs.

  microscale_readout.py run RUN_OUTPUT_DIR OUT_DIR [--config FILE] [--seed-nodelist F --seed-edgelist F]
  microscale_readout.py compare DIR_OR_LABEL=DIR,DIR,... [...] [--out FILE.csv]   # replicates show as min to max

`run` reads one run's output/{nodelist,edgelist}.csv and writes to OUT_DIR:
  checks.txt              fidelity checks (PASS/FAIL with detail)
  node_vector.parquet     one row per paper (size, size at age, survival, lifespan, diet, assigned traits)   [--keep-vectors]
  author_vector.parquet   one row per author                                                                [--keep-vectors]
  patterns.csv            population patterns (agents; seeds carry pre-simulation history)
  kernel_fits.csv, kernel_pooled.csv, tail_exponent.csv
  luck_vs_phenotype.csv   how much of the variation in citations lines up with the assigned traits
  summary.csv             the fixed list of measures reported for every run (measure,value)
`compare` prints the summary.csv files side by side (markdown) and optionally saves a CSV.

Conventions: citation age = year(citing) - year(cited); "by age X" counts ages 0..X; a value is NA (never 0)
when the paper has not reached the age; percentiles are nearest-rank over the stated population.
"""
import argparse, os, re, sys
import numpy as np, pandas as pd
from scipy import stats

AGE_CHECKPOINTS = (2, 5, 10)
BANDS = [(1, 1), (2, 2), (3, 5), (6, 10), (11, 20), (21, 10**6)]
BAND_LABELS = ["1", "2", "3-5", "6-10", "11-20", ">20"]
TRAITS = {"fitness": "fit_peak_value", "num_authors": "num_authors", "author_reputation": "initial_author_reputation"}


# ------------------------------------------------------------------ loading
def read_config(path):
    cfg = {}
    if path and os.path.exists(path):
        for line in open(path):
            m = re.match(r"^\s*([A-Za-z_]+)\s*=\s*([^;\s]*)", line)
            if m:
                cfg[m.group(1)] = m.group(2)
    return cfg


def load(run_dir):
    nl = pd.read_csv(f"{run_dir}/nodelist.csv", engine="pyarrow")
    el = pd.read_csv(f"{run_dir}/edgelist.csv", engine="pyarrow").astype(np.int64)
    n = int(nl.node_id.max()) + 1
    nl = nl.set_index("node_id").reindex(np.arange(n))          # ids absent from the nodelist become NaN rows
    nl["type"] = nl["type"].fillna("absent")
    return nl, el, n


# ------------------------------------------------------------------ fidelity checks
def run_checks(nl, el, n, cfg, seed_nl_path, seed_el_path):
    out = []
    def chk(name, ok, detail=""):
        out.append(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}".rstrip())
    types = nl["type"].to_numpy(); present = types != "absent"
    is_seed, is_agent = types == "seed", types == "agent"
    year = nl["year"].fillna(-1).to_numpy().astype(np.int64)
    src, tgt = el.source.to_numpy(), el.target.to_numpy()
    chk("node types are seed/agent", set(types) <= {"seed", "agent", "absent"}, str(pd.Series(types).value_counts().to_dict()))
    chk("edge endpoints are in the nodelist", present[src].all() and present[tgt].all())
    indeg = np.bincount(tgt, minlength=n); outdeg = np.bincount(src, minlength=n)
    chk("nodelist in_degree matches edgelist", np.array_equal(indeg[present], nl.in_degree.to_numpy()[present]))
    chk("nodelist out_degree matches edgelist", np.array_equal(outdeg[present], nl.out_degree.to_numpy()[present]))
    src_agent = is_agent[src]
    chk("agents never cite a later paper", not ((year[src] < year[tgt]) & src_agent).any(), f"n={int(((year[src] < year[tgt]) & src_agent).sum())}")
    chk("no self-loops", not (src == tgt).any())
    chk("no duplicate edges", not el.duplicated().any(), f"n={int(el.duplicated().sum())}")
    assigned = nl.assigned_out_degree.fillna(-1).to_numpy().astype(np.int64)
    short = (assigned[is_agent] != outdeg[is_agent]).sum()
    chk("realized out-degree equals assigned (agents)", short == 0, f"mismatches={int(short)}")
    # same-year citations: each of floor(frac * new agents) agents makes exactly one citation to a same-year agent
    if "same_year_citations" in cfg:
        frac = float(cfg["same_year_citations"])
        sy = src_agent & (year[src] == year[tgt])
        per_src = pd.Series(src[sy]).value_counts()
        born = pd.Series(year[is_agent]).value_counts().sort_index()
        have = pd.Series(year[per_src.index.to_numpy()]).value_counts().reindex(born.index).fillna(0)
        want = np.floor(born * frac)
        chk("same-year citations: floor(frac x new agents) agents per year", (have == want).all(), f"frac={frac}, years off={int((have != want).sum())}")
        chk("same-year citations: one each, agent to agent", (per_src == 1).all() and is_agent[tgt[sy]].all(), f"multi={int((per_src > 1).sum())}")
    if seed_nl_path and seed_el_path and os.path.exists(seed_nl_path) and os.path.exists(seed_el_path):
        s_nl = pd.read_csv(seed_nl_path, engine="pyarrow"); s_el = pd.read_csv(seed_el_path, engine="pyarrow").astype(np.int64)
        idcol = "node_id" if "node_id" in s_nl.columns else s_nl.columns[0]
        s_ids = s_nl[idcol].to_numpy(); s_year = s_nl["year"].to_numpy()
        s_deg = np.bincount(np.concatenate([s_el.iloc[:, 0], s_el.iloc[:, 1]]), minlength=s_ids.max() + 1)
        isolates = s_ids[s_deg[s_ids] == 0]
        dropped = np.setdiff1d(s_ids, np.flatnonzero(is_seed))
        chk("seeds dropped from the output are exactly the input's isolated nodes", np.array_equal(np.sort(dropped), np.sort(isolates)),
            f"dropped={dropped[:5].tolist()} isolates={isolates[:5].tolist()}")
        chk("no seeds beyond the input", len(np.setdiff1d(np.flatnonzero(is_seed), s_ids)) == 0)
        keep = np.isin(s_ids, np.flatnonzero(is_seed))
        chk("seed years equal the input", np.array_equal(year[s_ids[keep]], s_year[keep]))
        ss = el[~src_agent]
        a = s_el.iloc[:, :2].to_numpy(); b = ss.to_numpy()
        chk("seed-to-seed edges equal the input edgelist", len(a) == len(b) and np.array_equal(a[np.lexsort((a[:, 1], a[:, 0]))], b[np.lexsort((b[:, 1], b[:, 0]))]),
            f"output={len(b)} input={len(a)}")
    return out


# ------------------------------------------------------------------ per-node and per-author vectors
def build_vectors(nl, el, n):
    year = nl["year"].fillna(-1).to_numpy().astype(np.int64)
    src, tgt = el.source.to_numpy(), el.target.to_numpy()
    cage = year[src] - year[tgt]
    agents = nl[nl.type == "agent"]
    end_year = int(agents.year.max())
    v = nl.loc[nl.type != "absent", ["type", "year", "author_id"]].rename(columns={"year": "birth_year"}).copy()
    v["age_at_end"] = end_year - v.birth_year
    rdf = pd.DataFrame(dict(u=src, age=cage, young=cage < 5)); g = rdf.groupby("u")
    v["refs"] = g.size().reindex(v.index).fillna(0).astype(int)
    v["ref_age_median"] = g.age.median().reindex(v.index)
    v["price_index"] = g.young.mean().reindex(v.index)
    per_age = pd.DataFrame(dict(v_=tgt, age=cage)).groupby(["v_", "age"]).size().rename("c").reset_index()
    v["cites_final"] = per_age.groupby("v_").c.sum().reindex(v.index).fillna(0).astype(int)
    for X in AGE_CHECKPOINTS:
        cx = per_age[per_age.age <= X].groupby("v_").c.sum().reindex(v.index).fillna(0)
        v[f"cites_{X}y"] = cx.where(v.age_at_end >= X).astype("Int64")
    ga = per_age.groupby("v_").age
    v["first_cited_age"] = ga.min().reindex(v.index).astype("Int64")
    pk = per_age.sort_values(["v_", "c", "age"], ascending=[True, False, True]).drop_duplicates("v_").set_index("v_")
    v["peak_age"] = pk.age.reindex(v.index).astype("Int64")
    v["peak_cites"] = pk.c.reindex(v.index).astype("Int64")
    v["last_cited_age"] = ga.max().reindex(v.index).astype("Int64")
    traits = ["alpha", "pa_weight", "fit_weight", "num_authors_weight", "author_reputation_weight",
              "fit_peak_value", "num_authors", "initial_author_reputation", "final_author_reputation"]
    v = v.join(nl[traits])
    a = v.assign(is_seed=v.type == "seed").groupby("author_id").agg(
        papers=("birth_year", "size"), first_year=("birth_year", "min"), last_year=("birth_year", "max"),
        citations=("cites_final", "sum"), seed_papers=("is_seed", "sum"))
    a["career_length"] = a.last_year - a.first_year
    return v, a.reset_index()


def patterns(v, a):
    P = []
    def add(pattern, measure, s):
        s = pd.Series(s).dropna().astype(float)
        P.append(dict(pattern=pattern, measure=measure, n=len(s), mean=s.mean(), median=s.median(), p90=s.quantile(.9), p99=s.quantile(.99), max=s.max()))
    ag = v[v.type == "agent"]
    add("size", "cites_final", ag.cites_final)
    for X in AGE_CHECKPOINTS:
        add("size_at_age", f"cites_{X}y", ag[f"cites_{X}y"])
        add("survival", f"uncited_at_{X}y (share)", (ag[f"cites_{X}y"].dropna() == 0).astype(float))
    old = ag[ag.age_at_end >= 10]
    add("lifespan", "peak_age (age>=10)", old.peak_age)
    add("lifespan", "last_cited_age (age>=10)", old.last_cited_age)
    add("lifespan", "still_cited_in_final_year (age>=10, share)", (old.last_cited_age == old.age_at_end).astype(float))
    add("diet", "refs", ag.refs); add("diet", "ref_age_median", ag.ref_age_median); add("diet", "price_index", ag.price_index)
    aa = a[a.seed_papers == 0]
    add("productivity", "papers_per_author", aa.papers); add("productivity", "career_length", aa.career_length)
    return pd.DataFrame(P)


# ------------------------------------------------------------------ emergent measures
def nearest_rank(x, qs):
    return [int(q) for q in np.percentile(x, qs, method="nearest")]


def indegree_summary(v, el):
    d = v.cites_final.to_numpy(); da = v.loc[v.type == "agent", "cites_final"].to_numpy()
    S = {}
    for tag, x in (("all", d), ("agents", da)):
        mn, md, p90, p99, mx = nearest_rank(x, [0, 50, 90, 99, 100])
        S.update({f"indeg_{tag}_min": mn, f"indeg_{tag}_median": md, f"indeg_{tag}_p90": p90, f"indeg_{tag}_p99": p99, f"indeg_{tag}_max": mx,
                  f"never_cited_{tag}": float((x == 0).mean())})
    cited = np.sort(d[d > 0])[::-1]
    S["top10_share"] = float(cited[:10].sum() / cited.sum())
    S["top1pct_of_cited_share"] = float(cited[:len(cited) // 100].sum() / cited.sum())
    # the single maximum is one paper and very noisy, so also summarize the top 1% and top 0.1% of all papers by in-degree
    ranked = np.sort(d)[::-1]
    for tag, frac in (("1pct", 0.01), ("01pct", 0.001)):
        top = ranked[:max(1, int(np.ceil(frac * len(ranked))))]
        mn, md, p90, p99, mx = nearest_rank(top, [0, 50, 90, 99, 100])
        S.update({f"tail{tag}_n": len(top), f"tail{tag}_min": mn, f"tail{tag}_median": md, f"tail{tag}_q90": p90, f"tail{tag}_q99": p99, f"tail{tag}_max": mx})
    return S


def tail_exponent(indeg_agents):
    """Discrete power-law MLE (Clauset-Shalizi-Newman), k_min by minimum KS distance, >=100 tail points."""
    k = np.sort(indeg_agents[indeg_agents > 0]); best = None
    for kmin in np.unique(k):
        tail = k[k >= kmin]; m = len(tail)
        if m < 100: break
        g = 1 + m / np.log(tail / (kmin - 0.5)).sum()
        vals, cnt = np.unique(tail, return_counts=True)
        emp = 1 - np.concatenate([[0], np.cumsum(cnt)[:-1]]) / m
        D = np.abs(emp - ((vals - 0.5) / (kmin - 0.5)) ** (1 - g)).max()
        if best is None or D < best["ks"]:
            best = dict(gamma=g, gamma_se=(g - 1) / np.sqrt(m), kmin=int(kmin), n_tail=m, ks=D)
    return best


def kernel_fits(nl, el, n):
    """Effective attachment kernel A(k) by node-age band: weighted log-log slope over k ranges (bins with >=100 node-years)."""
    types = nl["type"].to_numpy(); year = np.full(n, 10**6, np.int64); ok = types != "absent"
    year[ok] = nl["year"].to_numpy()[ok].astype(np.int64)
    tgt = el.target.to_numpy(); eyear = year[el.source.to_numpy()]
    order = np.argsort(eyear, kind="stable"); tgt_s, ey_s = tgt[order], eyear[order]
    years_sim = np.unique(year[types == "agent"])
    cum = np.bincount(tgt_s[ey_s < years_sim[0]], minlength=n); parts = []
    for t in years_sim:
        lo, hi = np.searchsorted(ey_s, [t, t + 1]); new = np.bincount(tgt_s[lo:hi], minlength=n)
        elig = year < t
        if elig.any():
            mean_t = new[elig].sum() / elig.sum(); age = t - year[elig]
            band = np.full(age.shape, -1)
            for i, (a, b) in enumerate(BANDS): band[(age >= a) & (age <= b)] = i
            df = pd.DataFrame(dict(band=band, k=cum[elig], c=new[elig] / mean_t))
            parts.append(df.groupby(["band", "k"]).agg(node_years=("c", "size"), c_norm=("c", "sum")).reset_index())
        cum += new
    K = pd.concat(parts); rows, pooled = [], []
    for i, lab in enumerate(BAND_LABELS):
        p = K[K.band == i].groupby("k").agg(node_years=("node_years", "sum"), c_norm=("c_norm", "sum")).reset_index()
        p["A"] = p.c_norm / p.node_years; p.insert(0, "band", lab); pooled.append(p)
        row = dict(band=lab, node_years=int(p.node_years.sum()))
        for lo_k, hi_k in [(1, 5), (5, 20), (20, 10**9)]:
            q = p[(p.k >= lo_k) & (p.k < hi_k) & (p.node_years >= 100) & (p.A > 0)]
            tag = f"k{lo_k}_{hi_k if hi_k < 10**9 else 'inf'}"
            if len(q) < 4: row[f"slope_{tag}"] = np.nan; continue
            row[f"slope_{tag}"] = np.polyfit(np.log(q.k), np.log(q.A), 1, w=np.sqrt(q.node_years))[0]
        rows.append(row)
    return pd.DataFrame(rows), pd.concat(pooled)


def luck_vs_phenotype(v, rng, npairs=2_000_000):
    """How much of the variation in final citations lines up with the papers' assigned (cited-side) traits, holding
    birth year fixed? A description of how strongly outcomes depend on the system's own parts, not a forecasting model.
    pair test: of two same-year agents whose trait values and citation counts both differ, the share where the
    higher-trait paper has more citations (0.5 = a coin flip). A permuted trait should give 0.5.
    variation_unexplained_by_traits = 1 - R^2 of a rank regression of within-year citation rank on all three traits
    (an upper bound on 'luck': it only captures monotone, additive effects)."""
    ag = v[(v.type == "agent") & (v.age_at_end >= 5)].reset_index(drop=True)
    cohort_idx = [np.flatnonzero(ag.birth_year.to_numpy() == y) for y in np.unique(ag.birth_year)]
    out = ag.cites_final.to_numpy().astype(float); R = {}
    def pair_prob(x):
        conc = tot = 0
        for idx in cohort_idx:
            m = int(npairs * len(idx) / len(ag)); i = idx[rng.integers(len(idx), size=m)]; j = idx[rng.integers(len(idx), size=m)]
            dx, do = x[i] - x[j], out[i] - out[j]; ok = (dx != 0) & (do != 0)
            conc += int((np.sign(dx[ok]) == np.sign(do[ok])).sum()); tot += int(ok.sum())
        return conc / tot, tot
    for name, col in TRAITS.items():
        R[f"pair_prob_{name}"], R[f"pair_n_{name}"] = pair_prob(ag[col].to_numpy().astype(float))
    null = ag[TRAITS["fitness"]].to_numpy().astype(float).copy()
    for idx in cohort_idx: null[idx] = rng.permutation(null[idx])
    R["pair_prob_permuted_control"], _ = pair_prob(null)
    y = ag.groupby("birth_year").cites_final.rank(pct=True).to_numpy()
    X = np.column_stack([ag.groupby("birth_year")[c].rank(pct=True).to_numpy() for c in TRAITS.values()])
    def r2(Xm):
        A = np.column_stack([np.ones(len(y)), Xm]); b = np.linalg.lstsq(A, y, rcond=None)[0]
        return 1 - ((y - A @ b) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    for k, name in enumerate(TRAITS): R[f"r2_{name}"] = r2(X[:, [k]])
    R["r2_all_traits"] = r2(X); R["variation_unexplained_by_traits"] = 1 - R["r2_all_traits"]; R["n_agents"] = len(ag)
    return R


# ------------------------------------------------------------------ run / compare
def cmd_run(a):
    os.makedirs(a.out, exist_ok=True)
    cfg = read_config(a.config or os.path.join(os.path.dirname(os.path.abspath(a.run_dir)), "config.ini"))
    seed_nl = a.seed_nodelist or cfg.get("nodelist"); seed_el = a.seed_edgelist or cfg.get("edgelist")
    rng = np.random.default_rng(a.rng)
    nl, el, n = load(a.run_dir)
    lines = run_checks(nl, el, n, cfg, seed_nl, seed_el); open(f"{a.out}/checks.txt", "w").write("\n".join(lines) + "\n")
    v, au = build_vectors(nl, el, n)
    if a.keep_vectors:   # ~65 MB per run; regenerate any time
        v.reset_index().to_parquet(f"{a.out}/node_vector.parquet", index=False); au.to_parquet(f"{a.out}/author_vector.parquet", index=False)
    P = patterns(v, au); P.to_csv(f"{a.out}/patterns.csv", index=False)
    S = {"n_nodes": len(v), "n_seeds": int((v.type == "seed").sum()), "n_agents": int((v.type == "agent").sum()), "n_edges": len(el),
         "checks_failed": sum(l.startswith("FAIL") for l in lines)}
    S.update(indegree_summary(v, el))
    pm = P.set_index("measure")
    for X in AGE_CHECKPOINTS: S[f"uncited_at_{X}y"] = pm.loc[f"uncited_at_{X}y (share)", "mean"]
    S["cites_5y_median"] = pm.loc["cites_5y", "median"]; S["peak_age_median"] = pm.loc["peak_age (age>=10)", "median"]
    S["still_cited_final_year"] = pm.loc["still_cited_in_final_year (age>=10, share)", "mean"]
    S["ref_age_median"] = pm.loc["ref_age_median", "median"]; S["price_index_mean"] = pm.loc["price_index", "mean"]
    S["papers_per_author_median"] = pm.loc["papers_per_author", "median"]; S["papers_per_author_max"] = pm.loc["papers_per_author", "max"]
    t = tail_exponent(v.loc[v.type == "agent", "cites_final"].to_numpy())
    pd.DataFrame([t]).to_csv(f"{a.out}/tail_exponent.csv", index=False); S["tail_gamma"] = t["gamma"]; S["tail_kmin"] = t["kmin"]
    kf, kp = kernel_fits(nl, el, n); kf.to_csv(f"{a.out}/kernel_fits.csv", index=False); kp.to_csv(f"{a.out}/kernel_pooled.csv", index=False)
    for lab in ("3-5", "6-10", "11-20"): S[f"beta_eff_age_{lab}"] = kf.set_index("band").loc[lab, "slope_k20_inf"]
    L = luck_vs_phenotype(v, rng); pd.DataFrame([L]).to_csv(f"{a.out}/luck_vs_phenotype.csv", index=False)
    for k, val in L.items():
        if not k.startswith("pair_n_") and k != "n_agents": S[k] = val
    pd.DataFrame(list(S.items()), columns=["measure", "value"]).to_csv(f"{a.out}/summary.csv", index=False)
    print("\n".join(lines)); print(); print(pd.DataFrame(list(S.items()), columns=["measure", "value"]).to_string(index=False))


def _fmt(x):
    if pd.isna(x): return ""
    if x == int(x) and abs(x) >= 100: return f"{int(x):,}"
    return f"{x:.4g}"


def cmd_compare(a):
    """Each argument is a run dir, or LABEL=DIR,DIR,... for replicates; replicates show as min-max."""
    groups = {}
    for spec in a.dirs:
        label, paths = spec.split("=", 1) if "=" in spec else (os.path.basename(os.path.normpath(spec)), spec)
        groups[label] = pd.concat([pd.read_csv(f"{d}/summary.csv").set_index("measure").value.rename(i) for i, d in enumerate(paths.split(","))], axis=1)
    measures = list(next(iter(groups.values())).index); T = {}
    for label, S in groups.items():
        cells = []
        for m in measures:
            vals = S.loc[m].dropna().astype(float) if m in S.index else pd.Series(dtype=float)
            lo, hi = (vals.min(), vals.max()) if len(vals) else (np.nan, np.nan)
            cells.append(_fmt(lo) if _fmt(lo) == _fmt(hi) else f"{_fmt(lo)} to {_fmt(hi)}")
        T[label] = cells
    T = pd.DataFrame(T, index=measures)
    if a.out: T.to_csv(a.out)
    header = "| measure | " + " | ".join(T.columns) + " |"
    print(header); print("|" + "---|" * (len(T.columns) + 1))
    for m, row in T.iterrows(): print(f"| {m} | " + " | ".join(row) + " |")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run"); r.add_argument("run_dir"); r.add_argument("out"); r.add_argument("--config")
    r.add_argument("--seed-nodelist"); r.add_argument("--seed-edgelist"); r.add_argument("--rng", type=int, default=1); r.add_argument("--keep-vectors", action="store_true")
    c = sub.add_parser("compare"); c.add_argument("dirs", nargs="+"); c.add_argument("--out")
    a = ap.parse_args()
    {"run": cmd_run, "compare": cmd_compare}[a.cmd](a)
