struct DynamicBridges {
    int n;
    std::vector<int> dsu_2ecc;      // rep of 2-edge-connected component
    std::vector<int> dsu_cc;        // rep of connected component
    std::vector<int> dsu_cc_size;   // size of connected component
    std::vector<int> parent;        // parent in spanning forest (on CC level)
    std::vector<int> last_visit;    // for LCA-like merge_path
    int lca_iteration = 0;
    int bridges = 0;

    DynamicBridges(int n = 0) { init(n); }

    void init(int n_) {
        n = n_;
        dsu_2ecc.assign(n, 0);
        dsu_cc.assign(n, 0);
        dsu_cc_size.assign(n, 1);
        parent.assign(n, -1);
        last_visit.assign(n, 0);
        lca_iteration = 0;
        bridges = 0;
        for (int i = 0; i < n; ++i) {
            dsu_2ecc[i] = i;
            dsu_cc[i] = i;
            dsu_cc_size[i] = 1;
            parent[i] = -1;
        }
    }

    int find_2ecc(int v) {
        if (v == -1) return -1;
        if (dsu_2ecc[v] == v) return v;
        return dsu_2ecc[v] = find_2ecc(dsu_2ecc[v]);
    }

    int find_cc(int v) {
        v = find_2ecc(v);
        if (dsu_cc[v] == v) return v;
        return dsu_cc[v] = find_cc(dsu_cc[v]);
    }

    void make_root(int v) {
        v = find_2ecc(v);
        int root = v;
        int child = -1;
        while (v != -1) {
            int p = find_2ecc(parent[v]);
            parent[v] = child;
            int cc = dsu_cc[v];
            dsu_cc[v] = root;
            child = v;
            v = p;
        }
        dsu_cc_size[root] = dsu_cc_size[child];
    }

    void merge_path(int a, int b) {
        ++lca_iteration;
        std::vector<int> path_a, path_b;
        int lca = -1;

        while (lca == -1) {
            if (a != -1) {
                a = find_2ecc(a);
                path_a.push_back(a);
                if (last_visit[a] == lca_iteration) {
                    lca = a;
                    break;
                }
                last_visit[a] = lca_iteration;
                a = parent[a];
            }
            if (b != -1) {
                b = find_2ecc(b);
                path_b.push_back(b);
                if (last_visit[b] == lca_iteration) {
                    lca = b;
                    break;
                }
                last_visit[b] = lca_iteration;
                b = parent[b];
            }
        }

        // Compress 2ECCs along paths and decrease bridge count
        for (int v : path_a) {
            dsu_2ecc[v] = lca;
            if (v == lca) break;
            --bridges;
        }
        for (int v : path_b) {
            dsu_2ecc[v] = lca;
            if (v == lca) break;
            --bridges;
        }
    }

    // Add undirected edge (a,b)
    void add_edge(int a, int b) {
        a = find_2ecc(a);
        b = find_2ecc(b);
        if (a == b) return; // already in same 2ECC, no effect

        int ca = find_cc(a);
        int cb = find_cc(b);

        if (ca != cb) {
            // Edge connects two different CCs -> it's a new bridge
            ++bridges;

            // union by size on CC level
            if (dsu_cc_size[ca] > dsu_cc_size[cb]) {
                std::swap(a, b);
                std::swap(ca, cb);
            }

            make_root(a);
            parent[a] = b;
            dsu_cc[a] = b;
            dsu_cc_size[cb] += dsu_cc_size[a];
        } else {
            // same CC but different 2ECC: we introduce a cycle -> some bridges vanish
            merge_path(a, b);
        }
    }

    // Check if u and v are in the same 2-edge-connected component
    bool same_2ecc(int u, int v) {
        return find_2ecc(u) == find_2ecc(v);
    }

    // Check if u and v are in the same connected component
    bool same_cc(int u, int v) {
        return find_cc(u) == find_cc(v);
    }
};
