struct FlowEdge {
    int v, u;
    long long cap, flow = 0;
    FlowEdge(int v, int u, long long cap) : v(v), u(u), cap(cap) {}
};
 
struct Dinic {
    const long long INF = (long long)1e18;
    int n, s, t, m = 0;
    vector<FlowEdge> edges;
    vector<vector<int>> adj;
    vector<int> level, ptr;
    queue<int> q;
 
    Dinic(int n, int s, int t) : n(n), s(s), t(t), adj(n), level(n), ptr(n) {}
 
    void add_edge(int v, int u, long long cap) {
        edges.emplace_back(v, u, cap);
        edges.emplace_back(u, v, 0);
        adj[v].push_back(m);
        adj[u].push_back(m + 1);
        m += 2;
    }
    void add_undirected(int a, int b, long long c){ add_edge(a,b,c); add_edge(b,a,c); }
 
    bool bfs() {
        while (!q.empty()) {
            int v = q.front(); q.pop();
            for (int id : adj[v]) {
                if (edges[id].cap == edges[id].flow) continue;
                int u = edges[id].u;
                if (level[u] != -1) continue;
                level[u] = level[v] + 1;
                q.push(u);
            }
        }
        return level[t] != -1;
    }
 
    long long dfs(int v, long long pushed) {
        if (!pushed) return 0;
        if (v == t)  return pushed;
        for (int &cid = ptr[v]; cid < (int)adj[v].size(); ++cid) {
            int id = adj[v][cid], u = edges[id].u;
            if (level[v] + 1 != level[u]) continue;
            long long avail = edges[id].cap - edges[id].flow;
            if (avail <= 0) continue;
            long long tr = dfs(u, min(pushed, avail));
            if (!tr) continue;
            edges[id].flow += tr;
            edges[id ^ 1].flow -= tr;
            return tr;
        }
        return 0;
    }
 
    long long flow() {
        long long F = 0;
        for (;;) {
            fill(level.begin(), level.end(), -1);
            level[s] = 0;
            q = queue<int>(); q.push(s);
            if (!bfs()) break;
            fill(ptr.begin(), ptr.end(), 0);
            while (long long pushed = dfs(s, INF)) F += pushed;
        }
        return F;
    }
 
    vector<char> reachable_from_s_in_residual() {
        vector<char> vis(n, 0);
        queue<int> qq; qq.push(s); vis[s] = 1;
        while (!qq.empty()) {
            int v = qq.front(); qq.pop();
            for (int id : adj[v]) {
                auto &E = edges[id];
                if (!vis[E.u] && E.cap > E.flow) { vis[E.u] = 1; qq.push(E.u); }
            }
        }
        return vis;
    }
 
    vector<tuple<int,int,long long>> min_cut_edges() {
        auto inS = reachable_from_s_in_residual();
        vector<tuple<int,int,long long>> cut;
        for (int i = 0; i < (int)edges.size(); i += 2) {
            auto &E = edges[i];
            if (E.cap > 0 && inS[E.v] && !inS[E.u]) cut.emplace_back(E.v, E.u, E.cap);
        }
        return cut;
    }
 
    vector<tuple<int,int,long long>> positive_flow_edges() {
        vector<tuple<int,int,long long>> res;
        for (int i = 0; i < (int)edges.size(); i += 2) {
            auto &E = edges[i];
            long long f = E.flow;
            if (f > 0) res.emplace_back(E.v, E.u, f);
        }
        return res;
    }
};
