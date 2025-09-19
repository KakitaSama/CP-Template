int n,m;
vector<vector<int>> adj, adj_rev;
vector<bool> used;
vector<int> order, component;
vector<int> ct;
void dfs1(int v) {
    used[v] = true;
 
    for (auto u : adj[v])
        if (!used[u])
            dfs1(u);
 
    order.push_back(v);
}
 
void dfs2(int v) {
    used[v] = true;
    component.push_back(v);
    for (auto u : adj_rev[v])
        if (!used[u])
            dfs2(u);
}
