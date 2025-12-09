#include <bits/stdc++.h>
using namespace std;

/*
 * Persistent Segment Tree (Point updates, Range sum queries)
 *
 * - Values are long long, combine = sum.
 * - Each update creates a new version (root index).
 * - You can query any version independently.
 *
 * Indexing: array positions are [0..n-1].
 */

struct PersistentSegTree {
    struct Node {
        long long val;
        int left, right; // children indices in 'tree'
        Node(long long v = 0, int l = 0, int r = 0) : val(v), left(l), right(r) {}
    };

    int n;
    vector<Node> tree;  // all nodes for all versions
    vector<int> root;   // root[i] = root node index of version i

    PersistentSegTree() {}

    PersistentSegTree(int n_, const vector<long long>* arr = nullptr) {
        init(n_, arr);
    }

    void init(int n_, const vector<long long>* arr = nullptr) {
        n = n_;
        tree.clear();
        tree.reserve(40 * n); // heuristic, adjust if needed
        tree.push_back(Node()); // dummy 0-node (unused)

        root.clear();
        int r = build(0, n - 1, arr);
        root.push_back(r); // version 0
    }

    // Build initial tree (version 0) from arr or as all zeros
    int build(int L, int R, const vector<long long>* arr) {
        int id = newNode();
        if (L == R) {
            tree[id].val = (arr ? (*arr)[L] : 0LL);
            return id;
        }
        int mid = (L + R) >> 1;
        tree[id].left  = build(L, mid, arr);
        tree[id].right = build(mid + 1, R, arr);
        pull(id);
        return id;
    }

    int newNode() {
        tree.push_back(Node());
        return (int)tree.size() - 1;
    }

    void pull(int id) {
        tree[id].val = tree[tree[id].left].val + tree[tree[id].right].val;
    }

    // Create a new version by setting position 'pos' to 'newVal'
    // prevRoot: root of previous version
    int update_set(int prevRoot, int L, int R, int pos, long long newVal) {
        int id = newNode();
        tree[id] = tree[prevRoot]; // copy previous node (structure)
        if (L == R) {
            tree[id].val = newVal;
            return id;
        }
        int mid = (L + R) >> 1;
        if (pos <= mid) {
            tree[id].left = update_set(tree[prevRoot].left, L, mid, pos, newVal);
        } else {
            tree[id].right = update_set(tree[prevRoot].right, mid + 1, R, pos, newVal);
        }
        pull(id);
        return id;
    }

    // Create a new version by adding 'delta' to position 'pos'
    int update_add(int prevRoot, int L, int R, int pos, long long delta) {
        int id = newNode();
        tree[id] = tree[prevRoot];
        if (L == R) {
            tree[id].val = tree[prevRoot].val + delta;
            return id;
        }
        int mid = (L + R) >> 1;
        if (pos <= mid) {
            tree[id].left = update_add(tree[prevRoot].left, L, mid, pos, delta);
        } else {
            tree[id].right = update_add(tree[prevRoot].right, mid + 1, R, pos, delta);
        }
        pull(id);
        return id;
    }

    // Range sum query on version with root 'rt'
    long long query(int rt, int L, int R, int ql, int qr) const {
        if (qr < L || R < ql || rt == 0) return 0;
        if (ql <= L && R <= qr) return tree[rt].val;
        int mid = (L + R) >> 1;
        return query(tree[rt].left, L, mid, ql, qr)
             + query(tree[rt].right, mid + 1, R, ql, qr);
    }

    // Convenience wrappers that also store new version index in root[]
    // returns index of new version
    int set_point(int version, int pos, long long newVal) {
        int newRoot = update_set(root[version], 0, n - 1, pos, newVal);
        root.push_back(newRoot);
        return (int)root.size() - 1;
    }

    int add_point(int version, int pos, long long delta) {
        int newRoot = update_add(root[version], 0, n - 1, pos, delta);
        root.push_back(newRoot);
        return (int)root.size() - 1;
    }

    long long range_sum(int version, int l, int r) const {
        return query(root[version], 0, n - 1, l, r);
    }
};


int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n = 5;
    vector<long long> a = {1, 2, 3, 4, 5};

    PersistentSegTree pst(n, &a);

    // version 0: [1,2,3,4,5]
    cout << pst.range_sum(0, 0, 4) << "\n"; // 15

    // create version 1: set a[2] = 10
    int v1 = pst.set_point(0, 2, 10); // version 1: [1,2,10,4,5]
    cout << pst.range_sum(v1, 0, 4) << "\n"; // 22

    // version 2: from version 1, add +3 at position 0 → [4,2,10,4,5]
    int v2 = pst.add_point(v1, 0, 3);
    cout << pst.range_sum(v2, 0, 2) << "\n"; // 4+2+10 = 16

    // old version is unchanged:
    cout << pst.range_sum(0, 0, 2) << "\n"; // still 1+2+3 = 6
}
