#include <bits/stdc++.h>
using namespace std;

// Li Chao Tree (min) for lines y = m*x + b, x is long long
// To convert to **max Li Chao**:
//   1) Option A (easiest): store lines as usual, but negate values:
//        - When adding line y = m*x + b, insert line (-m, -b).
//        - When querying, return -query(x).
//   2) Option B (structural): change all `<` to `>` in add_line logic,
//      and change the neutral +INF to -INF in query.

struct LiChao {
    struct Line {
        long long m, b;
        // For MIN Li Chao: b = +INF is "identity" (always worse).
        // For MAX Li Chao (Option B): use b = -INF instead.
        Line(long long _m = 0, long long _b = (long long)4e18) : m(_m), b(_b) {}
        long long get(long long x) const { return m * x + b; }
    };

    struct Node {
        Line line;
        Node *left, *right;
        Node(const Line &l) : line(l), left(nullptr), right(nullptr) {}
    };

    long long X_MIN, X_MAX;
    Node *root;

    LiChao() : X_MIN(0), X_MAX(0), root(nullptr) {}

    LiChao(long long x_min, long long x_max) { init(x_min, x_max); }

    void init(long long x_min, long long x_max) {
        X_MIN = x_min;
        X_MAX = x_max;
        root = nullptr;
    }

    // Add line on full domain [X_MIN, X_MAX]
    void add_line(long long m, long long b) {
        Line newLine(m, b);
        // For MAX tree, Option A:
        //   Line newLine(-m, -b);
        add_line(root, X_MIN, X_MAX, newLine);
    }

    // Add line only on segment [l, r]
    void add_segment(long long m, long long b, long long l, long long r) {
        Line newLine(m, b);
        // For MAX tree, Option A:
        //   Line newLine(-m, -b);
        add_segment(root, X_MIN, X_MAX, newLine, l, r);
    }

    // Query minimum at x (for MIN Li Chao)
    // For MAX Li Chao, Option A: return -query(x)
    long long query(long long x) const {
        return query(root, X_MIN, X_MAX, x);
    }

private:
    // Core insertion for MIN Li Chao:
    //   - We keep at node->line the line that is best at mid.
    //   - We recurse with the worse line into the half where it might still win.
    // For MAX Li Chao (Option B):
    //   - replace all `<` with `>` in comparisons.
    void add_line(Node *&node, long long l, long long r, Line newLine) {
        if (!node) {
            node = new Node(newLine);
            return;
        }
        long long mid = (l + r) >> 1;

        Line cur = node->line;

        // MIN tree comparisons
        bool leftBetter  = newLine.get(l) < cur.get(l);
        bool midBetter   = newLine.get(mid) < cur.get(mid);
        bool rightBetter = newLine.get(r) < cur.get(r);

        // For MAX tree (Option B) use:
        //   bool leftBetter  = newLine.get(l) > cur.get(l);
        //   bool midBetter   = newLine.get(mid) > cur.get(mid);
        //   bool rightBetter = newLine.get(r) > cur.get(r);

        if (midBetter) {
            node->line = newLine;
            newLine = cur;
        }

        if (l == r) return;

        if (leftBetter != midBetter) {
            add_line(node->left, l, mid, newLine);
        } else if (rightBetter != midBetter) {
            add_line(node->right, mid + 1, r, newLine);
        }
    }

    void add_segment(Node *&node, long long l, long long r, Line newLine, long long ql, long long qr) {
        if (qr < l || r < ql) return;
        if (ql <= l && r <= qr) {
            add_line(node, l, r, newLine);
            return;
        }
        if (!node) node = new Node(Line()); // identity line
        long long mid = (l + r) >> 1;
        add_segment(node->left,  l,     mid, newLine, ql, qr);
        add_segment(node->right, mid+1, r,   newLine, ql, qr);
    }

    long long query(Node *node, long long l, long long r, long long x) const {
        if (!node) return (long long)4e18; // +INF for MIN tree
        // For MAX tree (Option B) use -4e18 as neutral and max() instead of min()
        long long res = node->line.get(x);
        if (l == r) return res;
        long long mid = (l + r) >> 1;
        if (x <= mid) return min(res, query(node->left, l, mid, x));
        else          return min(res, query(node->right, mid+1, r, x));
    }
};
