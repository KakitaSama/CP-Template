#include <bits/stdc++.h>
using namespace std;

// ---------- DSU with rollback ----------
struct RollbackDSU {
    int n;
    vector<int> parent, sz;
    vector<pair<int,int>> st; // (b, old_size_of_root_a)

    RollbackDSU(int n_ = 0) {
        init(n_);
    }

    void init(int n_) {
        n = n_;
        parent.resize(n + 1);
        sz.assign(n + 1, 1);
        for (int i = 1; i <= n; ++i) parent[i] = i;
        st.clear();
    }

    int find(int x) {
        while (x != parent[x]) x = parent[x];
        return x;
    }

    bool unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return false;
        if (sz[a] < sz[b]) swap(a, b); // attach b -> a
        st.push_back(make_pair(b, sz[a])); // store b and old sz[a]
        parent[b] = a;
        sz[a] += sz[b];
        return true;
    }

    int snapshot() {
        return (int)st.size();
    }

    void rollback(int snap) {
        while ((int)st.size() > snap) {
            pair<int,int> last = st.back();
            st.pop_back();
            int b = last.first;
            int old_sz_a = last.second;
            int a = parent[b];
            sz[a] = old_sz_a;
            parent[b] = b;
        }
    }
};

// ---------- Queries & globals ----------
struct Query {
    char type; // 'a' = add, 'r' = rem, 'c' = conn
    int u, v;
};

int N, Q;
vector<Query> queries;
vector<vector<pair<int,int>>> segtree; // segment tree nodes store edges (u, v)
vector<string> ans;
RollbackDSU dsu;

// ---------- Segment tree: add edge to interval ----------
void add_edge(int idx, int L, int R, int ql, int qr, pair<int,int> e) {
    if (qr < L || R < ql) return;
    if (ql <= L && R <= qr) {
        segtree[idx].push_back(e);
        return;
    }
    int mid = (L + R) / 2;
    add_edge(idx * 2, L, mid, ql, qr, e);
    add_edge(idx * 2 + 1, mid + 1, R, ql, qr, e);
}

// ---------- DFS over segment tree with rollback ----------
void dfs(int idx, int L, int R) {
    int snap = dsu.snapshot();

    // apply all edges active on this whole segment
    for (size_t i = 0; i < segtree[idx].size(); ++i) {
        int u = segtree[idx][i].first;
        int v = segtree[idx][i].second;
        dsu.unite(u, v);
    }

    if (L == R) {
        if (queries[L].type == 'c') {
            int u = queries[L].u;
            int v = queries[L].v;
            ans[L] = (dsu.find(u) == dsu.find(v)) ? "YES" : "NO";
        }
    } else {
        int mid = (L + R) / 2;
        dfs(idx * 2, L, mid);
        dfs(idx * 2 + 1, mid + 1, R);
    }

    dsu.rollback(snap);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    cin >> N >> Q;
    queries.resize(Q + 1);

    // Read queries (normalize edges so u < v)
    for (int i = 1; i <= Q; ++i) {
        string op;
        int u, v;
        cin >> op >> u >> v;
        if (u > v) swap(u, v); // normalize (u, v)

        char t;
        if (op[0] == 'a') t = 'a';      // add
        else if (op[0] == 'r') t = 'r'; // rem
        else t = 'c';                   // conn

        Query q;
        q.type = t;
        q.u = u;
        q.v = v;
        queries[i] = q;
    }

    // Build edge alive intervals
    map<pair<int,int>, int> last_add;
    vector<tuple<int,int,int,int>> intervals; // (l, r, u, v)

    for (int i = 1; i <= Q; ++i) {
        char t = queries[i].type;
        int u = queries[i].u;
        int v = queries[i].v;

        if (t == 'a') {
            last_add[make_pair(u, v)] = i;
        } else if (t == 'r') {
            pair<int,int> e = make_pair(u, v);
            int start = last_add[e];
            last_add.erase(e);
            if (start <= i - 1) {
                intervals.push_back(make_tuple(start, i - 1, u, v));
            }
        }
        // 'c' -> nothing for intervals
    }

    // Remaining edges: alive until end Q
    for (map<pair<int,int>, int>::iterator it = last_add.begin(); it != last_add.end(); ++it) {
        pair<int,int> e = it->first;
        int start = it->second;
        if (start <= Q) {
            intervals.push_back(make_tuple(start, Q, e.first, e.second));
        }
    }

    // Build segment tree with these intervals
    segtree.assign(4 * Q + 4, vector<pair<int,int>>());
    for (size_t i = 0; i < intervals.size(); ++i) {
        int l = std::get<0>(intervals[i]);
        int r = std::get<1>(intervals[i]);
        int u = std::get<2>(intervals[i]);
        int v = std::get<3>(intervals[i]);
        add_edge(1, 1, Q, l, r, make_pair(u, v));
    }

    ans.assign(Q + 1, "");
    dsu.init(N);

    if (Q > 0) dfs(1, 1, Q);

    // Output in order
    for (int i = 1; i <= Q; ++i) {
        if (queries[i].type == 'c') {
            cout << ans[i] << '\n';
        }
    }

    return 0;
}
