#include <bits/stdc++.h>
using namespace std;

struct Aho {
    static const int ALPHA = 26;

    struct Node {
        int next[ALPHA];
        int fail;
        int dict;          // dictionary link: nearest pattern node on fail-chain
        long long leafCnt; // patterns ending exactly here (with multiplicity)
        long long val;     // aggregated patterns on suffix chain

        Node() {
            memset(next, -1, sizeof(next));
            fail = 0;
            dict = 0;
            leafCnt = 0;
            val = 0;
        }
    };

    vector<Node> trie;

    Aho() {
        trie.clear();
        trie.push_back(Node()); // root = 0
    }

    void add_string(const string &s) {
        int node = 0;
        for (char ch : s) {
            int c = ch - 'a';
            if (trie[node].next[c] == -1) {
                trie[node].next[c] = (int)trie.size();
                trie.push_back(Node());
            }
            node = trie[node].next[c];
        }
        trie[node].leafCnt++; // multiplicity
    }

    void build() {
        queue<int> q;

        // Initialize root transitions
        for (int c = 0; c < ALPHA; ++c) {
            int nxt = trie[0].next[c];
            if (nxt != -1) {
                trie[nxt].fail = 0;
                q.push(nxt);
            } else {
                trie[0].next[c] = 0; // go back to root
            }
        }

        vector<int> order;
        order.reserve(trie.size());

        // BFS to build fail links and transitions
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            order.push_back(v);

            int f = trie[v].fail;
            for (int c = 0; c < ALPHA; ++c) {
                int nxt = trie[v].next[c];
                if (nxt != -1) {
                    trie[nxt].fail = trie[f].next[c];
                    q.push(nxt);
                } else {
                    trie[v].next[c] = trie[f].next[c];
                }
            }
        }

        // Build dictionary links and aggregated values
        // Root:
        trie[0].dict = 0;
        trie[0].val = trie[0].leafCnt; // normally 0

        for (int v : order) {
            int f = trie[v].fail;

            // dictionary link: nearest fail ancestor that has some pattern
            if (trie[f].leafCnt > 0) {
                trie[v].dict = f;
            } else {
                trie[v].dict = trie[f].dict;
            }

            // aggregate val: number of patterns on suffix chain of v
            // val[v] = leafCnt[v] + val[fail[v]]
            trie[v].val = trie[v].leafCnt + trie[f].val;
        }
    }

    // For each position i in text, return number of pattern endings at i
    vector<long long> process(const string &text) {
        vector<long long> res(text.size());
        int node = 0;
        for (int i = 0; i < (int)text.size(); ++i) {
            int c = text[i] - 'a';
            node = trie[node].next[c];
            // val[node] = total number of patterns that end here as suffixes
            res[i] = trie[node].val;
        }
        return res;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string t;
    int n;
    cin >> t >> n;

    vector<string> s(n);
    for (int i = 0; i < n; ++i) {
        cin >> s[i];
    }

    // Forward automaton: count pattern endings at each position
    Aho aho_fwd;
    for (int i = 0; i < n; ++i) {
        aho_fwd.add_string(s[i]);
    }
    aho_fwd.build();
    vector<long long> end_cnt = aho_fwd.process(t);

    // Reverse automaton: count pattern starts via reversed strings
    string t_rev = t;
    reverse(t_rev.begin(), t_rev.end());

    Aho aho_rev;
    for (int i = 0; i < n; ++i) {
        string r = s[i];
        reverse(r.begin(), r.end());
        aho_rev.add_string(r);
    }
    aho_rev.build();
    vector<long long> end_rev = aho_rev.process(t_rev);

    int m = (int)t.size();
    vector<long long> start_cnt(m);
    // Occurrence of reversed pattern ending at k in t_rev
    // <-> original pattern starting at pos = m - 1 - k in t
    for (int pos = 0; pos < m; ++pos) {
        start_cnt[pos] = end_rev[m - 1 - pos];
    }

    long long ans = 0;
    for (int p = 0; p + 1 < m; ++p) {
        ans += end_cnt[p] * start_cnt[p + 1];
    }

    cout << ans << "\n";
    return 0;
}
