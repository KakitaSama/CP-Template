#include <bits/stdc++.h>
using namespace std;

struct PersistentLiChaoMax {
    struct Line {
        long long m, b;  // y = m*x + b
        Line(long long _m = 0, long long _b = (long long)-4e18) : m(_m), b(_b) {}
        inline long long get(long long x) const { return m * x + b; }
    };

    struct Node {
        Line ln;
        Node *l, *r;
        Node() : ln(), l(nullptr), r(nullptr) {}
        Node(const Line& _ln) : ln(_ln), l(nullptr), r(nullptr) {}
    };

    static const int MAX_NODES = 4000000; // should be enough for Q<=1e5
    Node pool[MAX_NODES];
    int poolPtr = 0;

    long long X_MIN, X_MAX;

    PersistentLiChaoMax(long long x_min, long long x_max) {
        X_MIN = x_min;
        X_MAX = x_max;
    }

    inline Node* new_node() {
        return &pool[poolPtr++];
    }

    inline Node* new_node(const Node& other) {
        pool[poolPtr] = other;
        return &pool[poolPtr++];
    }

    Node* addLine(Node* root, long long m, long long b) {
        Line nl(m, b);
        return addLineInternal(root, X_MIN, X_MAX, nl);
    }

    long long query(Node* root, long long x) const {
        return queryInternal(root, X_MIN, X_MAX, x);
    }

private:
    Node* addLineInternal(Node* node, long long l, long long r, Line newLine) {
        if (!node) {
            node = new_node();
            node->ln = newLine;
            node->l = node->r = nullptr;
            return node;
        }

        Node* cur = new_node(*node);

        long long mid = (l + r) >> 1;

        Line low = cur->ln;
        Line high = newLine;

        bool leftBetter  = high.get(l) > low.get(l);
        bool midBetter   = high.get(mid) > low.get(mid);
        bool rightBetter = high.get(r) > low.get(r);

        if (midBetter) {
            cur->ln = high;
            newLine = low;
        } else {
            newLine = high;
        }

        if (l == r) return cur;

        if (leftBetter != midBetter) {
            cur->l = addLineInternal(cur->l, l, mid, newLine);
        } else if (rightBetter != midBetter) {
            cur->r = addLineInternal(cur->r, mid + 1, r, newLine);
        }
        return cur;
    }

    long long queryInternal(Node* node, long long l, long long r, long long x) const {
        if (!node) return (long long)-4e18;
        long long res = node->ln.get(x);
        if (l == r) return res;
        long long mid = (l + r) >> 1;
        if (x <= mid) return max(res, queryInternal(node->l, l, mid, x));
        else          return max(res, queryInternal(node->r, mid + 1, r, x));
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int Q;
    cin >> Q;

    const long long X_MIN = -1000000;
    const long long X_MAX =  1000000;

    PersistentLiChaoMax lichao(X_MIN, X_MAX);

    // index range: 0..1e6 according to statement
    static PersistentLiChaoMax::Node* version[1000005];
    for (int i = 0; i <= 1000000; ++i) version[i] = nullptr;

    int current = 0;

    // initial scroll: index=0, x=0, y=0
    version[0] = lichao.addLine(nullptr, 0, 0);

    while (Q--) {
        int type;
        cin >> type;
        if (type == 1) {
            int idx;
            cin >> idx;
            current = idx;
        } else if (type == 2) {
            long long X, Y;
            int idx;
            cin >> X >> Y >> idx;
            // parent = current
            version[idx] = lichao.addLine(version[current], X, Y);
            current = idx;
        } else if (type == 3) {
            long long V;
            cin >> V;
            long long ans = lichao.query(version[current], V);
            cout << ans << "\n";
        }
    }

    return 0;
}
