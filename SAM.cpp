#include <bits/stdc++.h>
using namespace std;

struct SuffixAutomaton {
    struct State {
        int link;              // suffix link
        int len;               // max length for this state
        int next[26];          // transitions
        long long occ;         // endpos count (occurrences of ALL substrings in this state)

        State() {
            link = -1;
            len  = 0;
            occ  = 0;
            memset(next, -1, sizeof(next));
        }
    };

    vector<State> st;
    int last;                  // state of whole string
    vector<long long> dpDistinct;
    vector<long long> dpAll;   // substrings with multiplicity

    SuffixAutomaton(int maxLen = 0) {
        st.reserve(2 * maxLen + 5);
        st.push_back(State()); // state 0 = root
        last = 0;
    }

    void extend(char ch) {
        int c = ch - 'a';
        int cur = (int)st.size();
        st.push_back(State());
        st[cur].len = st[last].len + 1;
        st[cur].occ = 1;       // new end position

        int p = last;
        while (p != -1 && st[p].next[c] == -1) {
            st[p].next[c] = cur;
            p = st[p].link;
        }

        if (p == -1) {
            st[cur].link = 0;
        } else {
            int q = st[p].next[c];
            if (st[p].len + 1 == st[q].len) {
                st[cur].link = q;
            } else {
                int clone = (int)st.size();
                st.push_back(st[q]);        // copy q
                st[clone].len = st[p].len + 1;
                st[clone].occ = 0;          // clone is not a real end

                while (p != -1 && st[p].next[c] == q) {
                    st[p].next[c] = clone;
                    p = st[p].link;
                }
                st[q].link = st[cur].link = clone;
            }
        }
        last = cur;
    }

    void build(const string &s) {
        for (char ch : s) extend(ch);
    }

    // Topological order by len (increasing)
    vector<int> topo_order() {
        int maxLen = 0;
        for (auto &x : st) maxLen = max(maxLen, x.len);

        vector<int> cnt(maxLen + 1);
        for (auto &x : st) cnt[x.len]++;

        for (int i = 1; i <= maxLen; ++i) cnt[i] += cnt[i - 1];

        vector<int> order(st.size());
        for (int i = (int)st.size() - 1; i >= 0; --i) {
            order[--cnt[st[i].len]] = i;
        }
        return order;
    }

    // Propagate occ over suffix links: endpos-count
    void build_occurrences() {
        auto ord = topo_order();        // increasing len
        for (int i = (int)ord.size() - 1; i > 0; --i) {
            int v = ord[i];
            int p = st[v].link;
            if (p >= 0) st[p].occ += st[v].occ;
        }
    }

    // ---------- DP for distinct substrings from each state ----------
    long long dfsDistinct(int v) {
        long long &res = dpDistinct[v];
        if (res != -1) return res;
        res = 0;
        for (int c = 0; c < 26; ++c) {
            int to = st[v].next[c];
            if (to == -1) continue;
            // one substring: the path adding this char
            // plus all longer substrings continuing from 'to'
            res += 1 + dfsDistinct(to);
        }
        return res;
    }

    void build_dpDistinct() {
        dpDistinct.assign(st.size(), -1);
        dfsDistinct(0); // root
    }

    // ---------- DP for all substrings (with multiplicity) ----------
    // from each state: sum over transitions (occ[to] + dpAll[to])
    long long dfsAllSub(int v) {
        long long &res = dpAll[v];
        if (res != -1) return res;
        res = 0;
        for (int c = 0; c < 26; ++c) {
            int to = st[v].next[c];
            if (to == -1) continue;
            res += st[to].occ;         // substring = path to 'to'
            res += dfsAllSub(to);      // longer substrings
        }
        return res;
    }

    void build_dpAll() {
        dpAll.assign(st.size(), -1);
        dfsAllSub(0);
    }

    // ---------- k-th distinct substring (1-indexed) ----------
    // assumes build_dpDistinct() has been called
    string kthDistinct(long long k) {
        string res;
        int v = 0;
        while (k > 0) {
            for (int c = 0; c < 26; ++c) {
                int to = st[v].next[c];
                if (to == -1) continue;
                long long block = 1 + dpDistinct[to];
                if (k > block) {
                    k -= block;
                } else {
                    // choose this char
                    res.push_back('a' + c);
                    k--;          // use the substring exactly equal to 'res'
                    if (k == 0) return res;
                    v = to;
                    break;
                }
            }
        }
        return res; // shouldn't reach if k is valid
    }

    // ---------- k-th substring (with multiplicity) ----------
    // assumes build_occurrences() and build_dpAll() have been called
    string kthWithMultiplicity(long long k) {
        string res;
        int v = 0;
        while (k > 0) {
            for (int c = 0; c < 26; ++c) {
                int to = st[v].next[c];
                if (to == -1) continue;
                long long block = st[to].occ + dpAll[to];
                if (k > block) {
                    k -= block;
                } else {
                    res.push_back('a' + c);
                    if (k <= st[to].occ) {
                        // the k-th is exactly this substring, no extension
                        return res;
                    }
                    k -= st[to].occ;
                    v = to;
                    break;
                }
            }
        }
        return res;
    }
};

// -------- CSES 2109: Substring Order II --------
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string s;
    long long k;
    cin >> s >> k;

    SuffixAutomaton sa((int)s.size());
    sa.build(s);
    sa.build_occurrences();
    sa.build_dpAll();

    // total substrings with multiplicity:
    // should equal n * (n + 1LL) / 2
    // long long total = sa.dpAll[0];

    string ans = sa.kthWithMultiplicity(k);
    cout << ans << "\n";
    return 0;
}
