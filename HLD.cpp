#include <bits/stdc++.h>
//src=https://github.com/ceciver/cp/blob/master/template/hld.cpp
using namespace std;

// ========= CONFIGURE T, unit, f =========

using T = long long;          // value type (can be int / long long / struct / etc.)
T unit = 0;                   // neutral element for sum
T f(T a, T b) { return a + b; } // merge function (here: sum)

// ========= Lazy Segment Tree (range add, range set, range sum) =========

struct Tree {
    int n;

    struct Node {
        T sum;
        T add_lazy;
        T set_lazy;
        bool has_set;
        Node(T s = 0) : sum(s), add_lazy(0), set_lazy(0), has_set(false) {}
    };

    vector<Node> st;

    Tree(int n_ = 0) {
        if (n_ > 0) init(n_);
    }

    Tree(const vector<T> &base) {
        init((int)base.size());
        build(1, 0, n, base);
    }

    void init(int n_) {
        n = n_;
        st.assign(4 * n, Node());
    }

    void build(int p, int l, int r, const vector<T> &base) {
        if (r - l == 1) {
            st[p].sum = (l < (int)base.size() ? base[l] : unit);
            return;
        }
        int m = (l + r) / 2;
        build(p * 2, l, m, base);
        build(p * 2 + 1, m, r, base);
        st[p].sum = f(st[p * 2].sum, st[p * 2 + 1].sum);
    }

    void apply_set(int p, int l, int r, T val) {
        st[p].sum = (r - l) * val;
        st[p].set_lazy = val;
        st[p].add_lazy = 0;
        st[p].has_set = true;
    }

    void apply_add(int p, int l, int r, T val) {
        st[p].sum += (r - l) * val;
        if (st[p].has_set) {
            st[p].set_lazy += val;
        } else {
            st[p].add_lazy += val;
        }
    }

    void push(int p, int l, int r) {
        if (r - l == 1) return;
        int m = (l + r) / 2;
        if (st[p].has_set) {
            apply_set(p * 2, l, m, st[p].set_lazy);
            apply_set(p * 2 + 1, m, r, st[p].set_lazy);
            st[p].has_set = false;
        }
        if (st[p].add_lazy != 0) {
            apply_add(p * 2, l, m, st[p].add_lazy);
            apply_add(p * 2 + 1, m, r, st[p].add_lazy);
            st[p].add_lazy = 0;
        }
    }

    // range set [ql, qr)
    void range_set(int ql, int qr, T val) { range_set(1, 0, n, ql, qr, val); }

    void range_set(int p, int l, int r, int ql, int qr, T val) {
        if (qr <= l || r <= ql) return;
        if (ql <= l && r <= qr) {
            apply_set(p, l, r, val);
            return;
        }
        push(p, l, r);
        int m = (l + r) / 2;
        range_set(p * 2, l, m, ql, qr, val);
        range_set(p * 2 + 1, m, r, ql, qr, val);
        st[p].sum = f(st[p * 2].sum, st[p * 2 + 1].sum);
    }

    // range add [ql, qr)
    void range_add(int ql, int qr, T val) { range_add(1, 0, n, ql, qr, val); }

    void range_add(int p, int l, int r, int ql, int qr, T val) {
        if (qr <= l || r <= ql) return;
        if (ql <= l && r <= qr) {
            apply_add(p, l, r, val);
            return;
        }
        push(p, l, r);
        int m = (l + r) / 2;
        range_add(p * 2, l, m, ql, qr, val);
        range_add(p * 2 + 1, m, r, ql, qr, val);
        st[p].sum = f(st[p * 2].sum, st[p * 2 + 1].sum);
    }

    // range query [ql, qr)
    T query(int ql, int qr) { return query(1, 0, n, ql, qr); }

    T query(int p, int l, int r, int ql, int qr) {
        if (qr <= l || r <= ql) return unit;
        if (ql <= l && r <= qr) return st[p].sum;
        push(p, l, r);
        int m = (l + r) / 2;
        return f(query(p * 2, l, m, ql, qr), query(p * 2 + 1, m, r, ql, qr));
    }
};

// ========= HLD =========

struct Edge {
    int u, v;
    T w;
};

template <bool VALS_EDGES>
struct HLD {
    int n;
    int currentTime = 0;

    vector<vector<int>> adjacency;
    vector<int> parent, subtreeSize, depth, head, position, inversePosition;
    Tree tree;

    HLD(const vector<vector<int>>& adjacency_)
        : n((int)adjacency_.size()),
          adjacency(adjacency_),
          parent(n, -1),
          subtreeSize(n, 1),
          depth(n, 0),
          head(n, 0),
          position(n, 0),
          inversePosition(n, 0),
          tree(n)
    {
        dfsSize(0);
        dfsHeavyLight(0);
    }

    void dfsSize(int vertex) {
        if (parent[vertex] != -1) {
            vector<int> &neighbors = adjacency[vertex];
            neighbors.erase(find(neighbors.begin(), neighbors.end(), parent[vertex]));
        }
        for (int &child : adjacency[vertex]) {
            parent[child] = vertex;
            depth[child] = depth[vertex] + 1;
            dfsSize(child);
            subtreeSize[vertex] += subtreeSize[child];
            if (!adjacency[vertex].empty() &&
                subtreeSize[child] > subtreeSize[adjacency[vertex][0]]) {
                swap(child, adjacency[vertex][0]);
            }
        }
    }

    void dfsHeavyLight(int vertex) {
        position[vertex] = currentTime;
        inversePosition[currentTime] = vertex;
        ++currentTime;

        if (adjacency[vertex].empty()) return;

        int heavyChild = adjacency[vertex][0];
        head[heavyChild] = head[vertex];
        dfsHeavyLight(heavyChild);

        for (size_t i = 1; i < adjacency[vertex].size(); ++i) {
            int child = adjacency[vertex][i];
            head[child] = child;
            dfsHeavyLight(child);
        }
    }

    void processPath(int u, int v, const function<void(int,int)>& op) {
        while (head[u] != head[v]) {
            if (depth[head[u]] > depth[head[v]]) swap(u, v);
            int start = head[v];
            op(position[start], position[v] + 1);
            v = parent[start];
        }
        if (depth[u] > depth[v]) swap(u, v);
        op(position[u] + (int)VALS_EDGES, position[v] + 1);
    }

    // ----- path ops -----

    void addPath(int u, int v, T delta) {
        processPath(u, v, [this, &delta](int l, int r) {
            tree.range_add(l, r, delta);
        });
    }

    void setPath(int u, int v, T x) {
        processPath(u, v, [this, &x](int l, int r) {
            tree.range_set(l, r, x);
        });
    }

    T queryPath(int u, int v) {
        T result = unit;
        processPath(u, v, [this, &result](int l, int r) {
            result = f(result, tree.query(l, r));
        });
        return result;
    }

    // ----- subtree ops -----

    void addSubtree(int v, T delta) {
        int L = position[v] + (int)VALS_EDGES;
        int R = position[v] + subtreeSize[v];
        tree.range_add(L, R, delta);
    }

    void setSubtree(int v, T x) {
        int L = position[v] + (int)VALS_EDGES;
        int R = position[v] + subtreeSize[v];
        tree.range_set(L, R, x);
    }

    T querySubtree(int v) {
        int L = position[v] + (int)VALS_EDGES;
        int R = position[v] + subtreeSize[v];
        return tree.query(L, R);
    }

    // ----- build helpers -----

    void build_from_nodes(const vector<T> &a) {
        vector<T> base(n);
        for (int v = 0; v < n; ++v) base[position[v]] = a[v];
        tree = Tree(base);
    }

    void build_from_edges(const vector<Edge> &edges) {
        vector<T> base(n, unit);
        for (const auto &e : edges) {
            int child = (parent[e.u] == e.v ? e.u : e.v);
            base[position[child]] = e.w;
        }
        tree = Tree(base);
    }

    // (extras: lca / distance / jump / kthOnPath could go here if you want)
};

// ========= Example usage in main =========

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, q;
    cin >> n >> q;          // number of nodes, queries
    vector<vector<int>> adj(n);
    for (int i = 0; i < n - 1; ++i) {
        int u, v;
        cin >> u >> v;
        --u; --v;          // convert to 0-based
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    vector<T> val(n);
    for (int i = 0; i < n; ++i) cin >> val[i]; // initial value on each vertex

    HLD<false> hld(adj);        // false -> values on vertices, not edges
    hld.build_from_nodes(val);  // build segtree over initial values

    // Example query types:
    // 1 u v     -> query sum on path u-v
    // 2 u v x   -> add x to path u-v
    // 3 v       -> query sum on subtree of v
    // 4 v x     -> add x to subtree of v
    // (all u, v are 1-based in input)

    while (q--) {
        int type;
        cin >> type;
        if (type == 1) {
            int u, v;
            cin >> u >> v;
            --u; --v;
            cout << hld.queryPath(u, v) << "\n";
        } else if (type == 2) {
            int u, v;
            T x;
            cin >> u >> v >> x;
            --u; --v;
            hld.addPath(u, v, x);
        } else if (type == 3) {
            int v;
            cin >> v;
            --v;
            cout << hld.querySubtree(v) << "\n";
        } else if (type == 4) {
            int v;
            T x;
            cin >> v >> x;
            --v;
            hld.addSubtree(v, x);
        }
    }
    return 0;
}
