# Example: a baseline and four treatments

Node-level readout (`microscale_readout.py run` on each run, then `compare` with replicate groups). 30 simulated years from the multifield
seed network (26,686 papers), about 1.17M papers and 21.9M citations at the end. Cells show the range across replicates; a single value
means the replicates agree to the digits shown.

| group | runs | simulator | what differs |
|---|---|---|---|
| old baseline | 6 | SASCA-ReSA 6.2.8 (5 runs), INSTAR 6.3.0 (1 run) | alpha and agent weights random, no planted nodes or cartels, `fully_random_citations` 0.05, score exponents 3 |
| fresh baseline | 3 | INSTAR 6.4.0 | same, but `fully_random_citations` 0.01 |
| pa_exponent 1 | 3 | INSTAR 6.3.0 | old baseline with `pa_exponent=1` (the other three scores still 3) |
| alpha 0.9 | 2 | SASCA-ReSA 6.2.8 | old baseline with alpha fixed at 0.9 |
| alpha 0.1 | 2 | SASCA-ReSA 6.2.8 | old baseline with alpha fixed at 0.1 |

**Compare each treatment with the old baseline**, which shares its random-citation setting. The fresh baseline differs from the old one in
one parameter (0.05 to 0.01); it leaves slightly more papers uncited and is the reference for later runs.

| measure | old baseline | fresh baseline | pa_exponent 1 | alpha 0.9 | alpha 0.1 |
|---|---|---|---|---|---|
| in-degree median | 6 | 5 | 6 | 8 | 4 |
| in-degree 90th percentile | 27 | 26 | 31 | 37 | 20 |
| in-degree 99th percentile | 145 to 154 | 154 to 158 | 168 | 139 | 176 to 178 |
| in-degree maximum | 126,378 to 229,528 | 114,788 to 126,390 | 97,035 to 106,755 | 44,057 to 52,505 | 125,803 to 134,044 |
| never cited | 0.1519 to 0.1525 | 0.1568 to 0.1575 | 0.1476 to 0.1478 | 0.1362 to 0.1363 | 0.1804 to 0.1807 |
| share of citations held by top 10 papers | 0.042 to 0.049 | 0.035 to 0.042 | 0.030 to 0.033 | 0.016 to 0.018 | 0.041 to 0.043 |
| share held by top 1% of cited papers | 0.419 to 0.429 | 0.438 to 0.439 | 0.346 to 0.347 | 0.244 | 0.544 to 0.545 |
| uncited at age 2 | 0.1112 to 0.1126 | 0.1238 to 0.1246 | 0.0951 to 0.0957 | 0.0676 | 0.1788 to 0.1805 |
| uncited at age 5 | 0.0166 to 0.0169 | 0.0215 to 0.0216 | 0.0121 to 0.0124 | 0.0054 to 0.0055 | 0.0463 to 0.0470 |
| uncited at age 10 | 0.0023 to 0.0026 | 0.0039 to 0.0041 | 0.0017 to 0.0018 | 0.0004 | 0.0113 |
| median age of citation peak | 4 | 4 | 4 | 5 | 3 |
| aged 10+ and still cited in final year | 0.550 to 0.555 | 0.536 to 0.539 | 0.578 to 0.581 | 0.728 to 0.733 | 0.363 to 0.369 |
| tail exponent gamma | 2.33 to 2.36 | 2.30 to 2.31 | 2.36 | 2.75 to 2.76 | 2.03 |
| effective kernel slope, age 6-10 | 1.07 to 1.08 | 1.08 to 1.09 | 1.08 | 0.94 | 1.11 |
| pair test: fitness | 0.599 to 0.602 | 0.602 to 0.604 | 0.600 to 0.602 | 0.608 to 0.610 | 0.621 to 0.627 |
| pair test: number of authors | 0.595 to 0.598 | 0.598 to 0.601 | 0.588 to 0.590 | 0.579 to 0.580 | 0.643 |
| pair test: author reputation | 0.640 to 0.642 | 0.637 to 0.642 | 0.630 to 0.631 | 0.609 to 0.611 | 0.654 to 0.655 |
| pair test: permuted control | 0.498 to 0.502 | 0.500 | 0.499 to 0.501 | 0.500 to 0.502 | 0.500 to 0.502 |
| variation not lined up with traits | 0.789 to 0.792 | 0.787 to 0.792 | 0.812 to 0.816 | 0.846 to 0.849 | 0.698 to 0.703 |

## Reading it
- Replicates agree closely within a condition on almost every measure. The exception is the single largest in-degree: which paper becomes
  the biggest hub is a matter of early chance, so it varies by up to a factor of 2 between replicates of the same condition.
- Locality (alpha) moves most things. High alpha (citing close to the generator) spreads citations out (smaller hubs, fewer uncited papers,
  a steeper tail, papers cited longer); low alpha concentrates them (the top 1% of cited papers hold more than half of all citations).
- Linear attachment for the in-degree score (`pa_exponent=1`) mostly trims the very top (smaller largest hubs, lower top-1% share) while the
  effective kernel and tail exponent stay close to the baseline's. The fitness and author-count scores were still cubed in these runs.
- In every condition the traits tilt outcomes modestly. Of two same-year papers that differ, the higher-trait one has more citations about
  60 to 65% of the time, and roughly 70 to 85% of the within-year variation does not line up with the traits. High alpha weakens the link and
  low alpha strengthens it (the second reading is an interpretation, not a tested mechanism).
