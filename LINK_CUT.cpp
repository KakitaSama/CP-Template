// ============================================================
// Link-Cut Tree (Splay-based) - Mega Template with Extras
// ------------------------------------------------------------
// CORE: dynamic forest (make_root, link, cut, connected, find_root)
// PATH FEATURES: node values, path add, path sum, path max
// LCA(root, u, v)
// EXTRAS:
//   - cut_parent(x): cut edge between x and its parent
//   - kth_on_path(u, v, k): return k-th node on path u->v (1-indexed)
//
// Node indices: 1..N, 0 = null
// Complexity: each operation O(log N) amortized
// ============================================================

struct LinkCutTree {
    // ===================== CUSTOMIZABLE DATA =====================
    struct Data {
        long long sum; // path sum
        long long mx;  // path max
    };

    static Data make_data(long long v = 0) {
        return {v, v};
    }

    static Data merge(const Data &a, const Data &b) {
        return {a.sum + b.sum, std::max(a.mx, b.mx)};
    }

    static const long long INF_NEG;

    // ===================== INTERNAL NODE STRUCT ==================
    struct Node {
        int ch[2]{0, 0}; // children in splay tree
        int p = 0;       // parent in splay/LCT
        bool rev = false;        // lazy path reversal
        bool has_lazy_add = false;
        long long lazy_add = 0;  // lazy for path-add

        int sz = 1;      // size of splay subtree (number of nodes)

        Data val = make_data(0);  // value at node
        Data path = make_data(0); // aggregate over splay subtree

        Node() = default;
    };

    std::vector<Node> t;

    LinkCutTree(int n = 0) {
        if (n) init(n);
    }

    void init(int n) {
        t.assign(n + 1, Node());
        for (int i = 1; i <= n; ++i) {
            t[i].val  = make_data(0);
            t[i].path = make_data(0);
            t[i].sz   = 1;
        }
    }

    // ===================== CORE: SPLAY / LCT BASICS ==============

    bool is_root(int x) {
        int p = t[x].p;
        return p == 0 || (t[p].ch[0] != x && t[p].ch[1] != x);
    }

    void push_up(int x) {
        t[x].sz = 1;
        Data res = t[x].val;

        int L = t[x].ch[0], R = t[x].ch[1];
        if (L) {
            t[x].sz += t[L].sz;
            res = merge(t[L].path, res);
        }
        if (R) {
            t[x].sz += t[R].sz;
            res = merge(res, t[R].path);
        }
        t[x].path = res;
    }

    void apply_rev(int x) {
        if (!x) return;
        std::swap(t[x].ch[0], t[x].ch[1]);
        t[x].rev ^= 1;
    }

    void apply_add_node(int x, long long delta) {
        if (!x) return;
        t[x].val.sum += delta;
        t[x].val.mx  += delta;
        t[x].path.sum += delta * t[x].sz;
        t[x].path.mx  += delta;
        t[x].lazy_add += delta;
        t[x].has_lazy_add = true;
    }

    void push_down(int x) {
        if (t[x].rev) {
            apply_rev(t[x].ch[0]);
            apply_rev(t[x].ch[1]);
            t[x].rev = false;
        }
        if (t[x].has_lazy_add) {
            long long d = t[x].lazy_add;
            if (t[x].ch[0]) apply_add_node(t[x].ch[0], d);
            if (t[x].ch[1]) apply_add_node(t[x].ch[1], d);
            t[x].lazy_add = 0;
            t[x].has_lazy_add = false;
        }
    }

    void rotate(int x) {
        int p = t[x].p;
        int g = t[p].p;
        bool is_right = (x == t[p].ch[1]);
        int b = t[x].ch[!is_right];

        if (!is_root(p)) {
            if (t[g].ch[0] == p) t[g].ch[0] = x;
            else if (t[g].ch[1] == p) t[g].ch[1] = x;
        }
        t[x].p = g;

        t[x].ch[!is_right] = p;
        t[p].p = x;

        t[p].ch[is_right] = b;
        if (b) t[b].p = p;

        push_up(p);
        push_up(x);
    }

    void push_all(int x) {
        static std::vector<int> st;
        st.clear();
        int y = x;
        while (!is_root(y)) {
            st.push_back(y);
            y = t[y].p;
        }
        st.push_back(y);
        for (int i = (int)st.size() - 1; i >= 0; --i) {
            push_down(st[i]);
        }
    }

    void splay(int x) {
        push_all(x);
        while (!is_root(x)) {
            int p = t[x].p;
            int g = t[p].p;
            if (!is_root(p)) {
                bool zigzig = ((t[p].ch[0] == x) == (t[g].ch[0] == p));
                if (zigzig) rotate(p);
                else        rotate(x);
            }
            rotate(x);
        }
        push_up(x);
    }

    // access(x): expose path root->x, return last accessed node
    int access(int x) {
        int last = 0;
        for (int y = x; y; y = t[y].p) {
            splay(y);
            t[y].ch[1] = last;
            push_up(y);
            last = y;
        }
        splay(x);
        return last;
    }

    // make_root(x): make x the root of its tree
    void make_root(int x) {
        access(x);
        apply_rev(x);
        push_up(x);
    }

    int find_root(int x) {
        access(x);
        while (t[x].ch[0]) {
            push_down(x);
            x = t[x].ch[0];
        }
        splay(x);
        return x;
    }

    bool connected(int u, int v) {
        if (u == v) return true;
        make_root(u);
        return find_root(v) == u;
    }

    // link(u,v): add edge u-v (must be different trees)
    bool link(int u, int v) {
        make_root(u);
        if (find_root(v) == u) return false; // already connected -> would form cycle
        t[u].p = v;
        return true;
    }

    // cut(u,v): remove edge u-v if it exists
    bool cut(int u, int v) {
        make_root(u);
        access(v);
        // now path u-v is in v's splay; edge u-v must be as below
        if (t[v].ch[0] != u || t[u].ch[1] != 0) return false;
        t[v].ch[0] = 0;
        t[u].p = 0;
        push_up(v);
        return true;
    }

    // ===================== PATH FEATURES ==========================

    // Set value of node x
    void set_val(int x, long long v) {
        access(x);
        t[x].val = make_data(v);
        push_up(x);
    }

    // Add delta to all nodes on path u-v
    void add_path(int u, int v, long long delta) {
        make_root(u);
        access(v);
        apply_add_node(v, delta);
        push_up(v);
    }

    // Query aggregate Data on path u-v
    Data query_path(int u, int v) {
        make_root(u);
        access(v);
        return t[v].path;
    }

    long long query_path_sum(int u, int v) {
        return query_path(u, v).sum;
    }

    long long query_path_max(int u, int v) {
        return query_path(u, v).mx;
    }

    // ===================== LCA (ROOTED) ===========================
    // LCA with respect to root r.
    // Assumes u and v are connected; returns 0 if not.
    int lca(int r, int u, int v) {
        if (!connected(u, v)) return 0;
        make_root(r);
        access(u);
        int l = access(v); // standard LCT trick
        return l;
    }

    // ===================== EXTRA STUFF ============================

    // 1) cut_parent(x): cut edge between x and its parent (if any)
    //
    // Meaning: detach x from the part of the tree "above" it.
    // Use when you know x has exactly one parent in the current forest.
    void cut_parent(int x) {
        access(x);
        int p = t[x].ch[0];
        if (!p) return;      // x is already a root
        t[x].ch[0] = 0;
        t[p].p = 0;
        push_up(x);
    }

    // 2) kth_on_path(u, v, k): k-th node from u to v (1-indexed)
    //
    // After make_root(u); access(v); the in-order traversal of v's splay tree
    // corresponds to nodes on the path u -> v.
    int kth_on_path(int u, int v, int k) {
        make_root(u);
        access(v); // v is root of splay representing path u->v
        int x = v;
        while (true) {
            push_down(x);
            int L = t[x].ch[0];
            int left_sz = (L ? t[L].sz : 0);
            if (k <= left_sz) {
                x = L;
            } else if (k == left_sz + 1) {
                splay(x);
                return x;
            } else {
                k -= left_sz + 1;
                x = t[x].ch[1];
            }
        }
    }
};

const long long LinkCutTree::INF_NEG = (long long)-4e18;

// ===================== EXAMPLE USAGE ==============================
/*
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n = 5;
    LinkCutTree lct(n);

    // set initial values
    for (int i = 1; i <= n; ++i) lct.set_val(i, i); // val[i] = i

    // build tree 1-2-3-4-5
    lct.link(1, 2);
    lct.link(2, 3);
    lct.link(3, 4);
    lct.link(4, 5);

    // path sum 2-5: 2+3+4+5 = 14
    cout << lct.query_path_sum(2, 5) << "\n";

    // add +10 on path 3-5 (3,4,5)
    lct.add_path(3, 5, 10);

    // new sum 2-5: 2 + 13 + 14 + 15 = 44
    cout << lct.query_path_sum(2, 5) << "\n";

    // max on path 1-5: max(1,2,23,24,25) = 25
    cout << lct.query_path_max(1, 5) << "\n";

    // LCA with root 1: LCA(4,5) = 4
    cout << "LCA(4,5) = " << lct.lca(1, 4, 5) << "\n";

    // k-th node on path 1->5:
    // path = [1,2,3,4,5], k=3 -> 3
    cout << "kth_on_path(1,5,3) = " << lct.kth_on_path(1, 5, 3) << "\n";

    // cut node 4 from its parent (which is 3)
    lct.cut_parent(4);
    cout << (lct.connected(3, 4) ? "YES" : "NO") << "\n"; // NO

    return 0;
}
*/
