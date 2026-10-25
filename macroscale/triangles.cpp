// triangles.cpp: local-structure metrics for a directed citation edgelist.
// Input: CSV with header "source,target" (source cites target).
// Build: g++ -O3 -fopenmp -o triangles triangles.cpp
// Usage: ./triangles edgelist.csv
//
// Reports undirected transitivity / clustering, a directed breakdown of
// triangles (transitive, cyclic, with reciprocal arcs), and the directed
// "out-closure" (do a paper's references cite each other) and "in-closure"
// (are the papers citing the same paper linked) coefficients.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>
#include <omp.h>

struct Edge {
    uint32_t a, b;  // a < b
    uint8_t f;      // bit0: a->b, bit1: b->a
};

static inline uint64_t choose2(uint64_t k) { return k * (k - 1) / 2; }

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s edgelist.csv\n", argv[0]);
        return 1;
    }
    FILE* fp = fopen(argv[1], "r");
    if (!fp) {
        perror("fopen");
        return 1;
    }

    // 1. read arcs, skipping the header and self loops
    std::vector<Edge> edges;
    uint32_t max_id = 0;
    uint64_t n_selfloops = 0;
    char line[256];
    while (fgets(line, sizeof line, fp)) {
        char* end;
        unsigned long long s = strtoull(line, &end, 10);
        if (end == line || *end != ',') continue;  // header or malformed
        char* p = end + 1;
        unsigned long long t = strtoull(p, &end, 10);
        if (end == p) continue;
        if (s == t) {
            n_selfloops++;
            continue;
        }
        uint32_t a = (uint32_t)std::min(s, t), b = (uint32_t)std::max(s, t);
        edges.push_back({a, b, (uint8_t)(s < t ? 1 : 2)});
        max_id = std::max(max_id, b);
    }
    fclose(fp);
    uint64_t n_rows = edges.size();

    // 2. merge duplicate pairs, OR-ing the direction flags
    std::sort(edges.begin(), edges.end(), [](const Edge& x, const Edge& y) {
        return x.a != y.a ? x.a < y.a : x.b < y.b;
    });
    size_t m = 0;
    for (size_t i = 0; i < edges.size(); i++) {
        if (m > 0 && edges[m - 1].a == edges[i].a && edges[m - 1].b == edges[i].b)
            edges[m - 1].f |= edges[i].f;
        else
            edges[m++] = edges[i];
    }
    edges.resize(m);
    uint64_t n_arcs = 0, n_recip_pairs = 0;
    for (const Edge& e : edges) {
        n_arcs += (e.f & 1) + ((e.f >> 1) & 1);
        if (e.f == 3) n_recip_pairs++;
    }

    // 3. degrees on raw ids
    std::vector<uint32_t> deg(max_id + 1, 0), outd(max_id + 1, 0), ind(max_id + 1, 0);
    for (const Edge& e : edges) {
        deg[e.a]++;
        deg[e.b]++;
        if (e.f & 1) { outd[e.a]++; ind[e.b]++; }
        if (e.f & 2) { outd[e.b]++; ind[e.a]++; }
    }

    // 4. rank nodes by (degree, id); relabel so edges point low rank -> high rank
    std::vector<uint32_t> ord;
    for (uint32_t i = 0; i <= max_id; i++)
        if (deg[i] > 0) ord.push_back(i);
    std::sort(ord.begin(), ord.end(), [&](uint32_t x, uint32_t y) {
        return deg[x] != deg[y] ? deg[x] < deg[y] : x < y;
    });
    size_t n = ord.size();
    std::vector<uint32_t> rank(max_id + 1, 0);
    for (size_t r = 0; r < n; r++) rank[ord[r]] = (uint32_t)r;

    // 5. forward adjacency (CSR). entry = (neighbor rank << 8) | flags, where
    // flags bit0 means "this node -> neighbor", bit1 means "neighbor -> this node"
    std::vector<uint64_t> off(n + 1, 0);
    for (const Edge& e : edges) off[std::min(rank[e.a], rank[e.b]) + 1]++;
    for (size_t i = 0; i < n; i++) off[i + 1] += off[i];
    std::vector<uint64_t> adj(m);
    {
        std::vector<uint64_t> pos(off.begin(), off.end() - 1);
        for (const Edge& e : edges) {
            uint32_t ra = rank[e.a], rb = rank[e.b];
            uint8_t f = e.f;
            uint32_t u = ra, v = rb;
            if (ra > rb) {
                u = rb;
                v = ra;
                f = (uint8_t)(((f & 1) << 1) | ((f >> 1) & 1));
            }
            adj[pos[u]++] = ((uint64_t)v << 8) | f;
        }
    }
    edges.clear();
    edges.shrink_to_fit();
#pragma omp parallel for schedule(dynamic, 1024)
    for (size_t u = 0; u < n; u++) std::sort(adj.begin() + off[u], adj.begin() + off[u + 1]);

    // 6. enumerate triangles u < v < w (ranks)
    std::vector<uint32_t> tri(n, 0), oc(n, 0), ic(n, 0);
    uint64_t T = 0, n_trans = 0, n_cyc = 0, n_recip_tri = 0;
#pragma omp parallel for schedule(dynamic, 256) \
    reduction(+ : T, n_trans, n_cyc, n_recip_tri)
    for (size_t u = 0; u < n; u++) {
        for (uint64_t ei = off[u]; ei < off[u + 1]; ei++) {
            uint32_t v = (uint32_t)(adj[ei] >> 8);
            uint32_t fuv = (uint32_t)(adj[ei] & 0xff);
            uint64_t i = off[u], iend = off[u + 1];
            uint64_t j = off[v], jend = off[v + 1];
            while (i < iend && j < jend) {
                uint32_t wi = (uint32_t)(adj[i] >> 8), wj = (uint32_t)(adj[j] >> 8);
                if (wi < wj) {
                    i++;
                } else if (wi > wj) {
                    j++;
                } else {
                    uint32_t w = wi;
                    uint32_t fuw = (uint32_t)(adj[i] & 0xff);
                    uint32_t fvw = (uint32_t)(adj[j] & 0xff);
                    int ou = (fuv & 1) + (fuw & 1);
                    int ov = ((fuv >> 1) & 1) + (fvw & 1);
                    int ow = ((fuw >> 1) & 1) + ((fvw >> 1) & 1);
                    int iu = ((fuv >> 1) & 1) + ((fuw >> 1) & 1);
                    int iv = (fuv & 1) + ((fvw >> 1) & 1);
                    int iw = (fuw & 1) + (fvw & 1);
                    T++;
                    if (fuv == 3 || fuw == 3 || fvw == 3)
                        n_recip_tri++;
                    else if (ou == 1 && ov == 1 && ow == 1)
                        n_cyc++;
                    else
                        n_trans++;
#pragma omp atomic
                    tri[u]++;
#pragma omp atomic
                    tri[v]++;
#pragma omp atomic
                    tri[w]++;
                    if (ou == 2) {
#pragma omp atomic
                        oc[u]++;
                    }
                    if (ov == 2) {
#pragma omp atomic
                        oc[v]++;
                    }
                    if (ow == 2) {
#pragma omp atomic
                        oc[w]++;
                    }
                    if (iu == 2) {
#pragma omp atomic
                        ic[u]++;
                    }
                    if (iv == 2) {
#pragma omp atomic
                        ic[v]++;
                    }
                    if (iw == 2) {
#pragma omp atomic
                        ic[w]++;
                    }
                    i++;
                    j++;
                }
            }
        }
    }

    // 7. summarize
    uint64_t triplets = 0, out_pairs = 0, in_pairs = 0, oc_sum = 0, ic_sum = 0;
    double lc_sum = 0, oc_mean_sum = 0, ic_mean_sum = 0;
    uint64_t lc_n = 0, oc_n = 0, ic_n = 0;
    for (size_t r = 0; r < n; r++) {
        uint32_t x = ord[r];
        uint64_t d = deg[x], o = outd[x], in = ind[x];
        triplets += choose2(d);
        out_pairs += choose2(o);
        in_pairs += choose2(in);
        oc_sum += oc[r];
        ic_sum += ic[r];
        if (d >= 2) { lc_sum += (double)tri[r] / choose2(d); lc_n++; }
        if (o >= 2) { oc_mean_sum += (double)oc[r] / choose2(o); oc_n++; }
        if (in >= 2) { ic_mean_sum += (double)ic[r] / choose2(in); ic_n++; }
    }

    printf("file=%s\n", argv[1]);
    printf("rows_read=%llu\n", (unsigned long long)n_rows);
    printf("self_loops_skipped=%llu\n", (unsigned long long)n_selfloops);
    printf("nodes_with_edges=%zu\n", n);
    printf("distinct_arcs=%llu\n", (unsigned long long)n_arcs);
    printf("undirected_edges=%zu\n", m);
    printf("reciprocal_pairs=%llu\n", (unsigned long long)n_recip_pairs);
    printf("triangles=%llu\n", (unsigned long long)T);
    printf("triangles_transitive=%llu\n", (unsigned long long)n_trans);
    printf("triangles_cyclic=%llu\n", (unsigned long long)n_cyc);
    printf("triangles_with_reciprocal_arc=%llu\n", (unsigned long long)n_recip_tri);
    printf("connected_triples=%llu\n", (unsigned long long)triplets);
    printf("global_transitivity=%.6f\n", triplets ? 3.0 * T / triplets : 0.0);
    printf("mean_local_clustering=%.6f\n", lc_n ? lc_sum / lc_n : 0.0);
    printf("out_closure_global=%.6f\n", out_pairs ? (double)oc_sum / out_pairs : 0.0);
    printf("out_closure_mean=%.6f\n", oc_n ? oc_mean_sum / oc_n : 0.0);
    printf("in_closure_global=%.6f\n", in_pairs ? (double)ic_sum / in_pairs : 0.0);
    printf("in_closure_mean=%.6f\n", ic_n ? ic_mean_sum / ic_n : 0.0);
    return 0;
}
