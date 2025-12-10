#include <bits/stdc++.h>
using namespace std;

struct SuffixTree {
    // ---------- Node structure ----------
    struct Node {
        int start;          // edge [start..end]
        int end;            // inclusive, INF means leaf -> use leafEnd
        int suffixLink;
        int next[256];
        int suffixIndex;    // for leaves: starting index of suffix

        Node(int s = -1, int e = -1) {
            start = s;
            end = e;
            suffixLink = -1;
            suffixIndex = -1;
            memset(next, -1, sizeof(next));
        }
    };

    static const int INF = 1000000007;

    string s;              // text + '$'
    vector<Node> st;       // nodes
    int root;

    // Ukkonen state
    int activeNode;
    int activeEdge;
    int activeLength;
    int remainingSuffixCount;
    int leafEnd;
    int lastNewNode;

    // ---------- ctor ----------
    SuffixTree(const string &str = "") {
        if (!str.empty()) build(str);
    }

    int newNode(int start, int end) {
        st.emplace_back(start, end);
        return (int)st.size() - 1;
    }

    int edgeLength(int v) const {
        if (v == root) return 0;
        int realEnd = (st[v].end == INF ? leafEnd : st[v].end);
        return realEnd - st[v].start + 1;
    }

    bool walkDown(int v) {
        int len = edgeLength(v);
        if (activeLength >= len) {
            activeEdge += len;
            activeLength -= len;
            activeNode = v;
            return true;
        }
        return false;
    }

    void extend(int pos) {
        leafEnd = pos;
        remainingSuffixCount++;
        lastNewNode = -1;
        unsigned char cur = (unsigned char)s[pos];

        while (remainingSuffixCount > 0) {
            if (activeLength == 0) activeEdge = pos;
            unsigned char edgeChar = (unsigned char)s[activeEdge];
            int nextIdx = st[activeNode].next[edgeChar];

            if (nextIdx == -1) {
                // create new leaf
                int leaf = newNode(pos, INF);
                st[leaf].suffixIndex = pos - remainingSuffixCount + 1;
                st[activeNode].next[edgeChar] = leaf;

                if (lastNewNode != -1) {
                    st[lastNewNode].suffixLink = activeNode;
                    lastNewNode = -1;
                }
            } else {
                int nxt = nextIdx;
                if (walkDown(nxt)) continue;

                int edgePos = st[nxt].start + activeLength;
                unsigned char edgeCharAtPos = (unsigned char)s[edgePos];

                if (edgeCharAtPos == cur) {
                    // already on edge, just move activeLength
                    if (lastNewNode != -1 && activeNode != root) {
                        st[lastNewNode].suffixLink = activeNode;
                        lastNewNode = -1;
                    }
                    activeLength++;
                    break;
                }

                // split edge
                int split = newNode(st[nxt].start,
                                    st[nxt].start + activeLength - 1);
                st[activeNode].next[edgeChar] = split;

                int leaf = newNode(pos, INF);
                st[leaf].suffixIndex = pos - remainingSuffixCount + 1;
                st[split].next[cur] = leaf;

                st[nxt].start += activeLength;
                st[split].next[(unsigned char)s[st[nxt].start]] = nxt;

                if (lastNewNode != -1)
                    st[lastNewNode].suffixLink = split;
                lastNewNode = split;
            }

            remainingSuffixCount--;
            if (activeNode == root && activeLength > 0) {
                activeLength--;
                activeEdge = pos - remainingSuffixCount + 1;
            } else if (activeNode != root) {
                activeNode = (st[activeNode].suffixLink == -1 ?
                              root : st[activeNode].suffixLink);
            }
        }
    }

    void build(const string &str) {
        s = str;
        if (s.empty()) return;
        if (s.back() != '$') s.push_back('$');

        st.clear();
        root = newNode(-1, -1);
        st[root].suffixLink = -1;

        activeNode = root;
        activeEdge = -1;
        activeLength = 0;
        remainingSuffixCount = 0;
        leafEnd = -1;
        lastNewNode = -1;

        for (int i = 0; i < (int)s.size(); ++i) {
            extend(i);
        }
    }

    // =========== helper: find node for pattern ===========
    int findNode(const string &pat) const {
        if (pat.empty()) return root;
        int v = root;
        int i = 0, n = (int)pat.size();
        while (i < n) {
            unsigned char c = (unsigned char)pat[i];
            int e = st[v].next[c];
            if (e == -1) return -1;
            int stPos = st[e].start;
            int enPos = (st[e].end == INF ? leafEnd : st[e].end);
            for (int k = stPos; k <= enPos && i < n; ++k, ++i) {
                if (s[k] != pat[i]) return -1;
            }
            v = e;
        }
        return v;
    }

    // ================= FEATURE: contains ==================
    bool contains(const string &pat) const {
        return findNode(pat) != -1;
    }

    // ============ FEATURE: occurrences ====================
    vector<int> occurrences(const string &pat) const {
        vector<int> res;
        int node = findNode(pat);
        if (node == -1) return res;

        function<void(int)> dfs = [&](int v) {
            if (st[v].suffixIndex != -1) {
                res.push_back(st[v].suffixIndex);
            }
            for (int c = 0; c < 256; ++c) {
                int to = st[v].next[c];
                if (to != -1) dfs(to);
            }
        };
        dfs(node);

        sort(res.begin(), res.end());
        int origLen = (s.back() == '$' ? (int)s.size() - 1 : (int)s.size());
        vector<int> filtered;
        for (int p : res) {
            if (p < origLen) filtered.push_back(p);
        }
        return filtered;
    }

    int count_occurrences(const string &pat) const {
        return (int)occurrences(pat).size();
    }

    // ===== FEATURE: Longest Repeated Substring (LRS) ======
    pair<int,int> dfsLRS(int v, int depth, int &bestLen, int &bestStart) const {
        int leaves = 0;
        int sample = -1;

        if (st[v].suffixIndex != -1) {
            leaves = 1;
            sample = st[v].suffixIndex;
        }

        for (int c = 0; c < 256; ++c) {
            int to = st[v].next[c];
            if (to == -1) continue;
            int realEnd = (st[to].end == INF ? leafEnd : st[to].end);
            int edgeLen = realEnd - st[to].start + 1;
            auto child = dfsLRS(to, depth + edgeLen, bestLen, bestStart);
            if (child.first > 0 && sample == -1)
                sample = child.second;
            leaves += child.first;
        }

        if (leaves >= 2 && depth > bestLen && sample != -1) {
            bestLen = depth;
            bestStart = sample;
        }
        return {leaves, sample};
    }

    string longest_repeated_substring() const {
        int bestLen = 0, bestStart = -1;
        dfsLRS(root, 0, bestLen, bestStart);

        int origLen = (s.back() == '$' ? (int)s.size() - 1 : (int)s.size());
        if (bestLen == 0 || bestStart < 0 || bestStart + bestLen > origLen)
            return "";
        return s.substr(bestStart, bestLen);
    }

    // ===== FEATURE: number of distinct substrings =========
    //
    // Count distinct substrings = sum over edges of edgeLength(e),
    // ignoring anything that involves the '$' terminal.
    //
    long long dfsDistinct(int v, int depth) const {
        long long res = 0;
        for (int c = 0; c < 256; ++c) {
            int to = st[v].next[c];
            if (to == -1) continue;
            int realEnd = (st[to].end == INF ? leafEnd : st[to].end);
            int edgeLen = realEnd - st[to].start + 1;

            // adjust if edge includes the '$' at the end
            int addLen = edgeLen;
            // if '$' lies inside this edge, cut it off
            int posDollar = (int)s.size() - 1;
            if (st[to].start <= posDollar && posDollar <= realEnd) {
                addLen = max(0, posDollar - st[to].start);
            }

            res += addLen;
            res += dfsDistinct(to, depth + edgeLen);
        }
        return res;
    }

    long long count_distinct_substrings() const {
        if (s.empty()) return 0;
        return dfsDistinct(root, 0);
    }
};
