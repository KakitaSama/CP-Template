#include <bits/stdc++.h>
// ==================== EULER & HAMILTON PACK ====================

#include <bits/stdc++.h>
using namespace std;

// ------------------------------------------------------------
// Eulerian paths & circuits
// ------------------------------------------------------------

enum class EulerType { NONE = 0, PATH = 1, CIRCUIT = 2 };

struct EulerResult {
    EulerType type;
    vector<int> path;  // vertex sequence of length m+1
};

// ---------- Undirected Euler (multigraph allowed) ----------
// Vertices: 0..n-1
// Input: n, edges as (u, v)
// Output: {type, path}, where path length = m+1 if exists.
EulerResult eulerUndirected(int n, const vector<pair<int,int>>& edges) {
    int m = (int)edges.size();
    vector<int> degree(n, 0);

    if (m == 0) {
        // convention: circuit of single isolated vertex 0 (if exists)
        return {EulerType::CIRCUIT, {0}};
    }

    // Compute degrees
    for (auto [u, v] : edges) {
        degree[u]++;
        degree[v]++;
    }

    // Build adjacency with edge IDs
    vector<vector<pair<int,int>>> g(n);
    for (int i = 0; i < m; ++i) {
        auto [u, v] = edges[i];
        g[u].push_back({v, i});
        g[v].push_back({u, i});
    }

    // Check connectivity among vertices with degree > 0
    int startVertex = -1;
    for (int i = 0; i < n; ++i) {
        if (degree[i] > 0) {
            startVertex = i;
            break;
        }
    }

    if (startVertex != -1) {
        vector<char> vis(n, 0);
        stack<int> st;
        st.push(startVertex);
        vis[startVertex] = 1;

        while (!st.empty()) {
            int v = st.top();
            st.pop();
            for (auto [to, id] : g[v]) {
                if (!vis[to] && degree[to] > 0) {
                    vis[to] = 1;
                    st.push(to);
                }
            }
        }

        for (int i = 0; i < n; ++i) {
            if (degree[i] > 0 && !vis[i]) {
                return {EulerType::NONE, {}};
            }
        }
    }

    // Count vertices of odd degree
    int oddCount = 0;
    int trailStart = startVertex;
    for (int i = 0; i < n; ++i) {
        if (degree[i] & 1) {
            oddCount++;
            trailStart = i;
        }
    }

    if (!(oddCount == 0 || oddCount == 2)) {
        return {EulerType::NONE, {}};
    }

    // Hierholzer
    vector<char> used(m, 0);
    vector<int> edgeIndex(n, 0);
    vector<int> stackVertices;
    vector<int> path;

    stackVertices.push_back(trailStart);

    while (!stackVertices.empty()) {
        int v = stackVertices.back();
        auto &adjList = g[v];

        while (edgeIndex[v] < (int)adjList.size() &&
               used[adjList[edgeIndex[v]].second]) {
            edgeIndex[v]++;
        }

        if (edgeIndex[v] == (int)adjList.size()) {
            // dead end -> add to path
            path.push_back(v);
            stackVertices.pop_back();
        } else {
            auto [to, id] = adjList[edgeIndex[v]++];
            if (used[id]) continue;
            used[id] = 1;
            stackVertices.push_back(to);
        }
    }

    if ((int)path.size() != m + 1) {
        return {EulerType::NONE, {}};
    }

    reverse(path.begin(), path.end());
    EulerType type = (oddCount ? EulerType::PATH : EulerType::CIRCUIT);
    return {type, path};
}

// ---------- Directed Euler ----------
// Conditions (for vertices with non-zero degree):
//   - weakly connected
//   - CIRCUIT: in[v] == out[v] for all v
//   - PATH   : exactly one vertex with out - in = 1 (start)
//              and one with in - out = 1 (end)
// Returns {type, path} with path length = m+1 if exists.
EulerResult eulerDirected(int n, const vector<pair<int,int>>& edges) {
    int m = (int)edges.size();
    if (m == 0) {
        return {EulerType::CIRCUIT, {0}};
    }

    vector<int> inDeg(n, 0), outDeg(n, 0);
    for (auto [u, v] : edges) {
        outDeg[u]++;
        inDeg[v]++;
    }

    // Weak connectivity check on underlying undirected graph
    vector<vector<int>> undirected(n);
    for (auto [u, v] : edges) {
        undirected[u].push_back(v);
        undirected[v].push_back(u);
    }

    int start0 = -1;
    for (int i = 0; i < n; ++i) {
        if (inDeg[i] + outDeg[i] > 0) {
            start0 = i;
            break;
        }
    }

    if (start0 != -1) {
        vector<char> vis(n, 0);
        stack<int> st;
        st.push(start0);
        vis[start0] = 1;

        while (!st.empty()) {
            int v = st.top();
            st.pop();
            for (int to : undirected[v]) {
                if (!vis[to] && (inDeg[to] + outDeg[to] > 0)) {
                    vis[to] = 1;
                    st.push(to);
                }
            }
        }

        for (int i = 0; i < n; ++i) {
            if ((inDeg[i] + outDeg[i] > 0) && !vis[i]) {
                return {EulerType::NONE, {}};
            }
        }
    }

    // Check degree conditions
    int start = -1, end = -1;
    for (int v = 0; v < n; ++v) {
        int diff = outDeg[v] - inDeg[v];
        if (diff == 1) {
            if (start != -1) return {EulerType::NONE, {}};
            start = v;
        } else if (diff == -1) {
            if (end != -1) return {EulerType::NONE, {}};
            end = v;
        } else if (diff != 0) {
            return {EulerType::NONE, {}};
        }
    }

    if ((start == -1) ^ (end == -1)) {
        // one of them is -1 but not both
        return {EulerType::NONE, {}};
    }

    EulerType type = EulerType::CIRCUIT;
    if (start == -1) {
        // Circuit case: any vertex with outgoing edge can be start
        for (int i = 0; i < n; ++i) {
            if (outDeg[i] > 0) {
                start = i;
                break;
            }
        }
    } else {
        type = EulerType::PATH;
    }

    // Hierholzer on directed graph
    struct Edge { int to, id; };
    vector<vector<Edge>> g(n);
    for (int i = 0; i < m; ++i) {
        auto [u, v] = edges[i];
        g[u].push_back({v, i});
    }

    vector<char> used(m, 0);
    vector<int> it(n, 0), st, path;
    st.push_back(start);

    while (!st.empty()) {
        int v = st.back();
        auto &adj = g[v];

        while (it[v] < (int)adj.size() && used[adj[it[v]].id]) {
            it[v]++;
        }

        if (it[v] == (int)adj.size()) {
            path.push_back(v);
            st.pop_back();
        } else {
            auto e = adj[it[v]++];
            if (used[e.id]) continue;
            used[e.id] = 1;
            st.push_back(e.to);
        }
    }

    if ((int)path.size() != m + 1) {
        return {EulerType::NONE, {}};
    }

    reverse(path.begin(), path.end());
    return {type, path};
}

// ------------------------------------------------------------
// Hamiltonian paths & cycles
// ------------------------------------------------------------

using Graph = vector<vector<int>>;

// -------- Hamiltonian PATH (bitmask DP, n <= ~20–22) --------
// g: adjacency list of an undirected simple graph
// out: one Hamiltonian path if exists
bool hamiltonPath(const Graph &g, vector<int> &out) {
    int n = (int)g.size();
    if (n == 0) return false;

    int maxMask = 1 << n;
    vector<vector<char>> dp(maxMask, vector<char>(n, 0));
    vector<vector<int>> parent(maxMask, vector<int>(n, -1));

    for (int v = 0; v < n; ++v) {
        dp[1 << v][v] = 1;
    }

    for (int mask = 1; mask < maxMask; ++mask) {
        for (int last = 0; last < n; ++last) {
            if (!dp[mask][last]) continue;
            for (int to : g[last]) {
                if (mask & (1 << to)) continue;
                int newMask = mask | (1 << to);
                if (!dp[newMask][to]) {
                    dp[newMask][to] = 1;
                    parent[newMask][to] = last;
                }
            }
        }
    }

    int fullMask = maxMask - 1;
    int endVertex = -1;
    for (int v = 0; v < n; ++v) {
        if (dp[fullMask][v]) {
            endVertex = v;
            break;
        }
    }

    if (endVertex == -1) return false;

    // Reconstruct path
    out.clear();
    int mask = fullMask;
    int cur = endVertex;
    while (cur != -1) {
        out.push_back(cur);
        int p = parent[mask][cur];
        mask ^= (1 << cur);
        cur = p;
    }
    reverse(out.begin(), out.end());
    return true;
}

// -------- Hamiltonian CYCLE (bitmask DP) --------
// g: undirected simple graph
// cyc: returns cycle (first == last)
bool hamiltonCycle(const Graph &g, vector<int> &cyc) {
    int n = (int)g.size();
    if (n == 0) return false;

    vector<vector<char>> adj(n, vector<char>(n, 0));
    for (int u = 0; u < n; ++u) {
        for (int v : g[u]) adj[u][v] = 1;
    }

    int maxMask = 1 << n;
    for (int start = 0; start < n; ++start) {
        vector<vector<char>> dp(maxMask, vector<char>(n, 0));
        vector<vector<int>> parent(maxMask, vector<int>(n, -1));

        dp[1 << start][start] = 1;

        for (int mask = 1; mask < maxMask; ++mask) {
            for (int last = 0; last < n; ++last) {
                if (!dp[mask][last]) continue;
                for (int to : g[last]) {
                    if (mask & (1 << to)) continue;
                    int newMask = mask | (1 << to);
                    if (!dp[newMask][to]) {
                        dp[newMask][to] = 1;
                        parent[newMask][to] = last;
                    }
                }
            }
        }

        int fullMask = maxMask - 1;
        for (int last = 0; last < n; ++last) {
            if (last == start) continue;
            if (!dp[fullMask][last]) continue;
            if (!adj[last][start])  continue; // must close cycle

            // Reconstruct cycle
            cyc.clear();
            int mask = fullMask;
            int cur = last;
            vector<int> tmp;
            while (cur != -1) {
                tmp.push_back(cur);
                int p = parent[mask][cur];
                mask ^= (1 << cur);
                cur = p;
            }
            reverse(tmp.begin(), tmp.end());
            cyc = tmp;
            cyc.push_back(start);  // close the cycle
            return true;
        }
    }
    return false;
}

// -------- Sufficient conditions for Hamiltonian cycle --------

// Dirac's condition:
// For an undirected simple graph with n >= 3,
// if deg(v) >= n/2 for all v, then G has a Hamiltonian cycle.
bool hasHamiltonCycleDirac(const Graph &g) {
    int n = (int)g.size();
    if (n < 3) return false;
    for (int i = 0; i < n; ++i) {
        if ((int)g[i].size() < (n + 1) / 2) return false;
    }
    return true;
}

// Ore's condition:
// For an undirected simple graph with n >= 3,
// if for every non-adjacent u != v, deg(u) + deg(v) >= n,
// then G has a Hamiltonian cycle.
bool hasHamiltonCycleOre(const Graph &g) {
    int n = (int)g.size();
    if (n < 3) return false;

    vector<vector<char>> adj(n, vector<char>(n, 0));
    for (int u = 0; u < n; ++u) {
        for (int v : g[u]) adj[u][v] = 1;
    }

    for (int u = 0; u < n; ++u) {
        for (int v = u + 1; v < n; ++v) {
            if (!adj[u][v] && !adj[v][u]) {
                if ((int)g[u].size() + (int)g[v].size() < n) {
                    return false;
                }
            }
        }
    }
    return true;
}

// -------- Hamiltonian PATH in DAG --------
// Returns true and fills `path` with a Hamiltonian path if
// there exists a topological ordering where every consecutive
// pair (path[i] -> path[i+1]) has a directed edge.
bool dagHamiltonPath(const Graph &g, vector<int> &path) {
    int n = (int)g.size();
    vector<int> indeg(n, 0);
    for (int u = 0; u < n; ++u) {
        for (int v : g[u]) indeg[v]++;
    }

    queue<int> q;
    for (int i = 0; i < n; ++i) {
        if (indeg[i] == 0) q.push(i);
    }

    vector<int> topo;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        topo.push_back(u);
        for (int v : g[u]) {
            if (--indeg[v] == 0) q.push(v);
        }
    }

    if ((int)topo.size() != n) {
        // not a DAG
        return false;
    }

    vector<vector<char>> adj(n, vector<char>(n, 0));
    for (int u = 0; u < n; ++u) {
        for (int v : g[u]) adj[u][v] = 1;
    }

    for (int i = 0; i + 1 < n; ++i) {
        if (!adj[topo[i]][topo[i + 1]]) return false;
    }

    path = topo;
    return true;
}

// ==================== /PACK ====================

// -------------------- Example usage --------------------
/*
int main() {
    int n;
    vector<pair<int,int>> edges;

    // Undirected Euler
    auto e1 = eulerUndirected(n, edges);
    if (e1.type == EulerType::NONE) {
        // no Euler path/circuit
    } else {
        // e1.path holds vertices of Euler path/circuit
    }

    // Directed Euler
    auto e2 = eulerDirected(n, edges);

    // Hamilton path / cycle
    Graph g(n);
    vector<int> hp, hc;
    bool hasHP = hamiltonPath(g, hp);
    bool hasHC = hamiltonCycle(g, hc);

    // Dirac / Ore check
    bool diracOK = hasHamiltonCycleDirac(g);
    bool oreOK   = hasHamiltonCycleOre(g);

    // DAG Hamilton path
    vector<int> dagHP;
    bool dagOK = dagHamiltonPath(g, dagHP);

    return 0;
}
*/
