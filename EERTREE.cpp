#include <bits/stdc++.h>
using namespace std;

// ============================================================
// EERTREE ESSENTIAL CORE
// ------------------------------------------------------------
// Handles: building from a string, counting distinct palindromes,
// total palindromes (with multiplicity), and accessing per-node info.
//
// Alphabet: 'a'..'z' (change ALPH + char_to_idx if needed).
// Build: O(n)
//
// Nodes:
//   1: imaginary root with len = -1
//   2: empty string root with len = 0
// Real palindromes are nodes 3..sz
// ============================================================

struct EertreeCore {
    static const int ALPH = 26;

    static int char_to_idx(char c) {
        return c - 'a';
    }

    struct Node {
        int next[ALPH];  // transitions by characters
        int link;        // suffix link
        int len;         // palindrome length
        int occ;         // occurrences (after finalize)
        int firstPos;    // last position (1-based in s)
        int num;         // number of palindromes in suffix chain

        Node() {
            memset(next, 0, sizeof(next));
            link = 0;
            len = 0;
            occ = 0;
            firstPos = -1;
            num = 0;
        }
    };

    vector<Node> tree;
    vector<char> s;      // s[0] = dummy, real string from s[1..N]
    string orig;
    int N = 0;

    int last;            // node of longest pal-suffix of current prefix
    int sz;              // number of used nodes
    long long totalPal;  // total pal substrings (with multiplicity)
    int distinctPal;     // distinct pal substrings
    vector<int> lastAtPos; // node index of longest pal-suffix at each pos (1-based)

    void init(int maxLen) {
        int maxNodes = maxLen + 3;
        tree.assign(maxNodes, Node());
        s.clear();
        s.reserve(maxLen + 1);
        s.push_back('#'); // dummy at index 0

        // root (len=-1), root0 (len=0)
        tree[1].len = -1; tree[1].link = 1;
        tree[2].len = 0;  tree[2].link = 1;

        last = 2;
        sz = 2;
        totalPal = 0;
        distinctPal = 0;

        lastAtPos.assign(maxLen + 1, 2);
    }

    void add_char_internal(char c, int pos) {
        int cId = char_to_idx(c);
        int cur = last;

        while (true) {
            int curlen = tree[cur].len;
            int mirror = pos - 1 - curlen;
            if (mirror >= 1 && s[mirror] == c) break;
            cur = tree[cur].link;
        }

        if (tree[cur].next[cId]) {
            // already exists
            last = tree[cur].next[cId];
            tree[last].occ++;
        } else {
            // create new node
            last = ++sz;
            tree[sz].len = tree[cur].len + 2;
            tree[sz].firstPos = pos;
            tree[sz].occ = 1;
            tree[cur].next[cId] = sz;

            if (tree[sz].len == 1) {
                tree[sz].link = 2;
            } else {
                int p = tree[cur].link;
                while (true) {
                    int plen = tree[p].len;
                    int mirror = pos - 1 - plen;
                    if (mirror >= 1 && s[mirror] == c) {
                        tree[sz].link = tree[p].next[cId];
                        break;
                    }
                    p = tree[p].link;
                }
            }

            tree[sz].num = tree[tree[sz].link].num + 1;
        }

        lastAtPos[pos] = last;
    }

    void finalize_tree() {
        // push occurrences up suffix links
        for (int v = sz; v >= 3; --v) {
            tree[tree[v].link].occ += tree[v].occ;
        }
        distinctPal = max(0, sz - 2);
        totalPal = 0;
        for (int v = 3; v <= sz; ++v) totalPal += tree[v].occ;
    }

    void build(const string &str) {
        orig = str;
        N = (int)orig.size();
        if (N == 0) {
            init(0);
            return;
        }
        init(N);
        for (int i = 0; i < N; ++i) {
            s.push_back(orig[i]);
            int pos = (int)s.size() - 1;
            add_char_internal(orig[i], pos);
        }
        finalize_tree();
    }

    // ===== Core queries =====
    int distinct_palindromes() const { return distinctPal; }
    long long total_palindromes() const { return totalPal; }

    // node of longest palindromic suffix at position i (0-based in orig)
    int node_longest_suffix_at(int i) const {
        int pos1 = i + 1;
        if (pos1 < 1 || pos1 >= (int)lastAtPos.size()) return 2;
        return lastAtPos[pos1];
    }
};
// ============================================================
// EERTREE FEATURES LAYER
// ------------------------------------------------------------
// Built on top of EertreeCore.
// Adds:
//   - longest palindromic substring
//   - palindromic counts per prefix / per position
//   - count by length
//   - enumerate_palindromes()
//   - counts with constraints
//   - minimal palindromic partition (DP)
// ============================================================

struct Eertree : EertreeCore {
    // longest palindrome
    int bestLen = 0;
    int bestPos = 0; // 1-based end position in s

    // prefix info:
    // pal_end_at[i]: #palindromes ending at i
    // pal_pref[i]:   #palindromes in prefix [0..i]
    vector<long long> pal_end_at;
    vector<long long> pal_pref;

    // count by length
    vector<long long> count_by_length;

    // partition DP stuff
    bool hasPartitionDP = false;
    bool partitionDPEnabled = false;
    vector<int> dp, bestPrev, bestNode, seriesAns;
    int minPartsWhole = -1;

    // seriesLink/diff needed for fast pal-partition DP
    vector<int> diff, seriesLink;

    // override/init hook
    void init_features(int maxLen) {
        bestLen = 0;
        bestPos = 0;
        pal_end_at.assign(maxLen, 0);
        pal_pref.assign(maxLen, 0);
        count_by_length.assign(maxLen + 1, 0);
        diff.assign((int)tree.size(), 0);
        seriesLink.assign((int)tree.size(), 0);
    }

    // we wrap the core init to also prepare feature arrays
    void init(int maxLen) {
        EertreeCore::init(maxLen);
        init_features(maxLen);
    }

    // we wrap the add_char_internal to update features too
    void add_char_with_features(char c, int pos) {
        int before_last = last;
        EertreeCore::add_char_internal(c, pos);
        int v = last;

        // if new node created, its index might be > before_last
        if (v > before_last && tree[v].len > 0) {
            // new palindrome node
            count_by_length[tree[v].len]++;
            // setup diff + seriesLink (for partition DP)
            int link_v = tree[v].link;
            diff[v] = tree[v].len - tree[link_v].len;
            if (diff[v] == diff[link_v]) seriesLink[v] = seriesLink[link_v];
            else seriesLink[v] = link_v;
        }

        // update longest palindrome
        if (tree[v].len > bestLen) {
            bestLen = tree[v].len;
            bestPos = pos;
        }

        // prefix counts
        int i0 = pos - 1;
        long long newEnd = tree[v].num;
        pal_end_at[i0] = newEnd;
        pal_pref[i0] = newEnd + (i0 > 0 ? pal_pref[i0 - 1] : 0);

        // partition DP incremental update
        if (partitionDPEnabled) update_partition_dp(pos);
    }

    // DP update for minimal pal partition
    void update_partition_dp(int pos) {
        const int INF = (int)1e9;
        dp[pos] = INF;
        for (int v = last; tree[v].len > 0; v = seriesLink[v]) {
            int sl = seriesLink[v];
            int prev = pos - (tree[sl].len + diff[v]);
            if (prev < 0) continue;

            seriesAns[v] = dp[prev];
            if (diff[v] == diff[tree[v].link]) {
                seriesAns[v] = min(seriesAns[v], seriesAns[tree[v].link]);
            }
            int cand = seriesAns[v] + 1;
            if (cand < dp[pos]) {
                dp[pos] = cand;
                bestPrev[pos] = pos - tree[v].len;
                bestNode[pos] = v;
            }
        }
    }

    // finalize hook: call core finalize then features remain consistent
    void finalize() {
        EertreeCore::finalize_tree();
    }

    // build WITHOUT partition DP
    void build(const string &str) {
        orig = str;
        N = (int)orig.size();
        if (N == 0) {
            init(0);
            finalize();
            hasPartitionDP = false;
            minPartsWhole = 0;
            return;
        }
        init(N);
        partitionDPEnabled = false;
        hasPartitionDP = false;

        for (int i = 0; i < N; ++i) {
            s.push_back(orig[i]);
            int pos = (int)s.size() - 1;
            add_char_with_features(orig[i], pos);
        }
        finalize();
    }

    // build WITH partition DP
    void build_with_partition_dp(const string &str) {
        orig = str;
        N = (int)orig.size();
        if (N == 0) {
            init(0);
            finalize();
            hasPartitionDP = true;
            dp.assign(1, 0);
            minPartsWhole = 0;
            return;
        }
        init(N);

        const int INF = (int)1e9;
        partitionDPEnabled = true;
        hasPartitionDP = true;

        dp.assign(N + 1, INF);
        dp[0] = 0;
        bestPrev.assign(N + 1, -1);
        bestNode.assign(N + 1, -1);
        seriesAns.assign((int)tree.size(), INF);

        for (int i = 0; i < N; ++i) {
            s.push_back(orig[i]);
            int pos = (int)s.size() - 1;
            add_char_with_features(orig[i], pos);
        }
        finalize();
        minPartsWhole = dp[N];
        partitionDPEnabled = false;
    }

    // ================= FEATURES =================

    // longest palindrome [l, r] 0-based
    pair<int,int> longest_pal_interval() const {
        if (bestLen == 0) return {-1, -1};
        int endPos1 = bestPos;
        int startPos1 = endPos1 - bestLen + 1;
        int l = startPos1 - 1;
        int r = endPos1 - 1;
        return {l, r};
    }

    string longest_pal_string() const {
        auto [l, r] = longest_pal_interval();
        if (l < 0) return "";
        return orig.substr(l, r - l + 1);
    }

    // per-prefix stats:
    //  pal_end_at[i], pal_pref[i]

    struct PalInfo {
        int start;
        int length;
        int occ;
    };

    // enumerate distinct palindromes
    vector<PalInfo> enumerate_palindromes() const {
        vector<PalInfo> res;
        for (int v = 3; v <= sz; ++v) {
            int L = tree[v].len;
            int endPos1 = tree[v].firstPos;
            int startPos1 = endPos1 - L + 1;
            int start0 = startPos1 - 1;
            if (start0 < 0 || start0 + L > N) continue;
            res.push_back({start0, L, tree[v].occ});
        }
        return res;
    }

    // count distinct palindromes with occ >= K
    int count_distinct_with_occ_at_least(int K) const {
        int cnt = 0;
        for (int v = 3; v <= sz; ++v)
            if (tree[v].occ >= K) cnt++;
        return cnt;
    }

    // count distinct palindromes with length in [L, R]
    int count_distinct_with_length_between(int L, int R) const {
        int cnt = 0;
        for (int v = 3; v <= sz; ++v) {
            int len = tree[v].len;
            if (len >= L && len <= R) cnt++;
        }
        return cnt;
    }

    // partition DP queries
    int min_pal_parts() const {
        if (!hasPartitionDP) return -1;
        return minPartsWhole;
    }

    int min_pal_cuts() const {
        if (!hasPartitionDP) return -1;
        return max(0, minPartsWhole - 1);
    }

    const vector<int>& min_pal_partition_dp() const {
        return dp;
    }

    vector<pair<int,int>> get_pal_partition_segments() const {
        vector<pair<int,int>> segs;
        if (!hasPartitionDP) return segs;
        int pos = N;
        while (pos > 0) {
            int prev = bestPrev[pos];
            if (prev < 0) break;
            int l = prev;
            int r = pos - 1;
            segs.push_back({l, r});
            pos = prev;
        }
        reverse(segs.begin(), segs.end());
        return segs;
    }
};
