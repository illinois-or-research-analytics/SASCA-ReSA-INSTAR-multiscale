# SASCA-ReSA-INSTAR-multiscale

Multiscale readout for simulated citation networks from [SASCA-ReSA](https://github.com/illinois-or-research-analytics/SASCA-ReSA)
and [SASCA-ReSA-INSTAR](https://github.com/illinois-or-research-analytics/SASCA-ReSA-INSTAR).

**Principle.** These are abstract generative models of a complex system. We do not try to match real citation data or to forecast
anything. We measure the grown network at several scales, report which patterns are present and how robust they are across
replicates, and rank conditions against a baseline. Spread across replicates is a result (it shows how path-dependent the system is),
not noise to explain away.

| scale | what is measured | status |
|---|---|---|
| **microscale** (a paper) | size, size at age, survival, lifespan, diet, productivity; attachment kernel and tail exponent; luck vs phenotype | done: `microscale/` |
| **mesoscale** (clusters of papers) | IKC dense-core clusters (coverage, depth, modularity, hubs inside); Leiden CPM clusters with conductance | done: `mesoscale/` |
| **macroscale** (the whole network) | concentration of citations (top-10 and top-1% shares, Gini) and triangle counts (transitivity, local clustering) | done: `macroscale/` |

## Quick start
```bash
pip install -r requirements.txt            # numpy, pandas, scipy, pyarrow
python tests/smoke_test.py                 # builds a toy simulation and runs everything on it (about a minute)

# one run: RUN_OUTPUT_DIR holds nodelist.csv and edgelist.csv written by the simulator
python microscale/microscale_readout.py run RUN_OUTPUT_DIR results/myrun
# add --keep-vectors to also save the per-paper and per-author tables (about 65 MB for 1.2M papers)

# compare runs; LABEL=DIR,DIR,... groups replicates and shows each measure as "min to max"
python microscale/microscale_readout.py compare \
    baseline=results/r0,results/r1,results/r2  treatment=results/t0,results/t1 --out results/compare.csv
```
On a Slurm cluster: `slurm/microscale_readout.sbatch` (one array task per run; usage in its header).

The script reads `config.ini` next to `RUN_OUTPUT_DIR` (override with `--config`) to find the seed network and
`same_year_citations`; without it, the seed and same-year checks are skipped.

## All three scales on a Slurm cluster
```bash
cp slurm/pipeline.env.example my.env          # set RUNS_DIR, OUT_DIR, PYTHON, IKC_PY, MODULES, SBATCH_EXTRA
slurm/run_pipeline.sh -e my.env RUN1 RUN2     # RUN = a folder under RUNS_DIR with output/{nodelist,edgelist}.csv
python multiscale.py compare baseline=OUT_DIR/RUN1,OUT_DIR/RUN2 treatment=OUT_DIR/RUN3   # when the jobs have finished
```
For each run the driver submits the stages `micro`, `ikc`, `leiden` and `macro` in parallel and a `merge` job that waits for them
(it runs even if a stage failed, so a late or failed Leiden job does not block the rest). Results go to `OUT_DIR/RUN/{micro,meso,macro}/`
and are merged into `OUT_DIR/RUN/summary_all.csv` (every measure as `<scale>.<measure>`; `stages.txt` says which outputs were found).
Typical cost for a 1.2M-paper, 22M-citation run: micro about 2 minutes, IKC about 3, triangles about 15 seconds, and Leiden 15-25
minutes, so the whole readout takes about half an hour. The meso stage needs `igraph`, `leidenalg` and `networkit`; the macro stage
needs the triangle program (`make -C macroscale`, built automatically by the driver if missing). The IKC program is not in this
repository: set `IKC_PY` to `Illinois/clustering/eleanor/code/IKC.py` from `chackoge/ERNIE_Plus` (commit `d3325da946`).
Individual stages can be run directly: `mesoscale/mesoscale_readout.py ikc|leiden`, `macroscale/macroscale_readout.py run`.
`python tests/smoke_test_pipeline.py` tests everything except IKC (and the triangle program, unless built) on a toy network.

## What `run` writes
| file | contents |
|---|---|
| `checks.txt` | fidelity checks, PASS/FAIL with detail |
| `summary.csv` | the fixed list of measures reported for every run (`measure,value`) |
| `patterns.csv` | population patterns (agents only; mean, median, 90th, 99th, max) |
| `tail_exponent.csv`, `kernel_fits.csv`, `kernel_pooled.csv`, `luck_vs_phenotype.csv` | detail behind the summary |
| `node_vector.parquet`, `author_vector.parquet` | one row per paper / author (only with `--keep-vectors`) |

## The measures
Conventions: citation age is the citing paper's year minus the cited paper's year. "By age X" counts ages 0 to X. A value that needs
a paper to have reached an age is **NA, never 0**, for papers that have not reached it, so young papers cannot masquerade as poorly
cited ones. Percentiles are nearest-rank. Population patterns use agents only, because seed papers carry history from before the simulation.

**Checks** (inputs the simulator should reproduce): node degrees match the edgelist; no citations to later papers, self-loops or duplicate
edges; every agent's realized out-degree equals its assigned quota; same-year citations are exactly `floor(same_year_citations x new agents)`
per year, one each, agent to agent; seeds and seed-to-seed edges match the input network, except nodes the input leaves isolated (the
simulator builds its node set from the edgelist).

**Patterns, in the spirit of pattern-oriented ecology** (each paper is an "individual"):
| pattern | measures |
|---|---|
| size | citations received (final); in-degree min, median, 90th, 99th, max; share never cited; share of all citations held by the top 10 papers and by the top 1% of cited papers |
| size at age | citations by age 2, 5 and 10 |
| survival | share still uncited at age 2, 5 and 10 |
| lifespan | age of the citation peak; age when last cited; share of papers aged 10 or more still cited in the final year |
| diet | median age of the papers a paper cites; share of references under 5 years old (largely set by the simulator's recency table, so mostly an input check) |
| productivity | papers per author; career length (authors' output is designed to follow Lotka's law, so also an input check) |

**Emergent structure**
- **Tail exponent** gamma of the in-degree distribution, N(k) ~ k^(-gamma): discrete power-law fit with k_min chosen by minimum KS distance
  (Clauset, Shalizi & Newman 2009), agents only. Not the same thing as the attachment exponent beta that sets how agents score candidates.
- **Effective attachment kernel** (Jeong, Neda & Barabasi 2003): for each simulated year, how the citations a paper receives grow with
  the in-degree it already has, within bands of paper age so that recency does not mask the trend. Reported as the log-log slope for
  in-degree 20 or more, by age band (`beta_eff_age_*`; about 1.1 means roughly linear).
- **Luck vs phenotype**: how much of the variation in citations lines up with a paper's assigned traits (fitness, number of authors,
  initial author reputation), holding birth year fixed. This describes how strongly outcomes depend on the system's own parts; it is not a forecasting model.
  - *Pair test* (`pair_prob_*`): of two same-year papers whose trait and citations both differ, the share where the higher-trait paper has
    more citations. 0.5 is a coin flip; a permuted-trait control (`pair_prob_permuted_control`) should give 0.5.
  - `variation_unexplained_by_traits`: 1 minus R^2 of a regression of within-year citation rank on the three traits' ranks. An upper bound
    on luck, since it only captures monotone, additive effects. Agents aged 5 or more.

**Mesoscale.** *IKC* (iterative k-core clustering) extracts dense cores recursively from the densest down; a core is kept only if it
passes a positive-modularity test and every paper in it has at least `k` connections inside it, and extracted papers leave the graph,
so its clusters are not plain k-cores. Reported at `k`=20: number of clusters, share of papers in them, depth, modularity of the largest
cluster, the share of all citations it receives, and how many of the 100 / 1,000 most-cited papers it contains. *Leiden CPM* at one
resolution: number of clusters, size distribution, and **conductance** (the share of a cluster's connections that lead outside it; lower
means better separated), reported as the median and as the share of clusters below 0.5, by cluster count and weighted by papers.
Conductance has a floor set by cluster size (a perfect 8-paper clique in a graph of mean degree 37 cannot go below about 0.8), so read it
together with the share of papers in the clusters it describes (`leiden_ge10_paper_share`).

**Macroscale.** *Concentration*: share of all citations received by the 10 most-cited papers and by the most-cited 1% of cited papers, and
the Gini coefficient of in-degree (0 = equal, 1 = one paper has all). *Triangles*: undirected closed triples (directions ignored,
duplicate edges merged), per paper, plus transitivity, mean local clustering and directed out-/in-closure. Transitivity is dominated by
open triples around hubs, so triangles per paper and mean local clustering are the steadier comparisons. The triangle program was
validated against an independent sparse-matrix count on a real 50k-paper slice.

## Example
[docs/example_results.md](docs/example_results.md): a 30-year baseline against four treatments, replicate ranges.

## Notes
- Percentiles over "all nodes" include papers never cited, so the minimum is 0.
- The largest in-degree is a single paper's value and varies a lot between replicates; the 99th percentile and the top-10 and top-1% shares are steadier.
- Developed against SASCA-ReSA 6.2.8 and SASCA-ReSA-INSTAR 6.3 and 6.4 output (`nodelist.csv` columns: `node_id, type, year, alpha, pa_weight,
  fit_weight, num_authors_weight, author_reputation_weight, fit_peak_value, in_degree, out_degree, assigned_out_degree, num_authors, initial_author_reputation,
  final_author_reputation, author_id`; `edgelist.csv`: `source, target` meaning source cites target).
