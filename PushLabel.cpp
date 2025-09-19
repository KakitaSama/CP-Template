#include <bits/stdc++.h>
using namespace std;

struct PushRelabel {
    struct E { int to, rev; long long cap, cap0; };  // cap0 = original forward capacity
    int n;
    vector<vector<E>> g;
    vector<long long> ex;
    vector<int> h, it;

    PushRelabel(int N=0){ init(N); }
    void init(int N){ n=N; g.assign(n,{}); }

    // directed edge u->v with capacity c
    void add_edge(int u,int v,long long c){
        E a{v, (int)g[v].size(), c, c};
        E b{u, (int)g[u].size(), 0, 0};
        g[u].push_back(a); g[v].push_back(b);
    }
    // convenience for undirected unit edges (used in Police Chase)
    void add_undirected(int u,int v,long long c){ add_edge(u,v,c); add_edge(v,u,c); }

    void push(int u, E &e){
        if (ex[u]==0 || e.cap==0 || h[u]!=h[e.to]+1) return;
        long long f = min(ex[u], e.cap);
        e.cap -= f;
        g[e.to][e.rev].cap += f;
        ex[u] -= f;
        ex[e.to] += f;
    }
    void relabel(int u){
        int mh = INT_MAX;
        for (auto &e: g[u]) if (e.cap>0) mh = min(mh, h[e.to]);
        if (mh<INT_MAX) h[u] = mh + 1;
    }
    void discharge(int u){
        while (ex[u]>0){
            if (it[u] == (int)g[u].size()){
                relabel(u);
                it[u] = 0;
            } else {
                push(u, g[u][it[u]]);
                if (g[u][it[u]].cap==0 || h[u]!=h[g[u][it[u]].to]+1) ++it[u];
            }
        }
    }
    long long maxflow(int s,int t){
        ex.assign(n,0); h.assign(n,0); it.assign(n,0);
        h[s]=n;
        for (auto &e: g[s]){ // preflow
            long long f=e.cap;
            if (f>0){ e.cap-=f; g[e.to][e.rev].cap+=f; ex[e.to]+=f; }
        }
        vector<int> vs; vs.reserve(n-2);
        for (int i=0;i<n;i++) if (i!=s && i!=t) vs.push_back(i);

        // relabel-to-front
        for (size_t i=0;i<vs.size();){
            int u=vs[i];
            int old = h[u];
            discharge(u);
            if (h[u]>old){ rotate(vs.begin(), vs.begin()+i, vs.begin()+i+1); i=0; }
            else ++i;
        }
        long long F=0;
        for (auto &e: g[s]) F += g[e.to][e.rev].cap; // total pushed from s
        return F;
    }

    // ---- Min-cut helpers (INSIDE struct) ----

    // S-side: nodes reachable from s in residual graph (cap>0)
    vector<char> mincut_side(int s) const {
        vector<char> vis(n,0);
        queue<int> q; q.push(s); vis[s]=1;
        while(!q.empty()){
            int v=q.front(); q.pop();
            for (auto const &e: g[v]) if (!vis[e.to] && e.cap>0){
                vis[e.to]=1; q.push(e.to);
            }
        }
        return vis;
    }

    // Directed edges crossing the min cut (from S to T), with original capacity
    vector<tuple<int,int,long long>> mincut_edges(int s) const {
        auto inS = mincut_side(s);
        vector<tuple<int,int,long long>> cut;
        for (int v=0; v<n; ++v) if (inS[v]){
            for (auto const &e: g[v]){
                if (e.cap0>0 && !inS[e.to]) cut.emplace_back(v, e.to, e.cap0);
            }
        }
        return cut;
    }

    // (Optional) forward edges with positive flow
    vector<tuple<int,int,long long>> flow_edges() const {
        vector<tuple<int,int,long long>> res;
        for (int u=0; u<n; ++u){
            for (auto const &e: g[u]){
                if (e.cap0>0){
                    long long sent = g[e.to][e.rev].cap; // reverse residual = flow sent
                    if (sent>0) res.emplace_back(u,e.to,sent);
                }
            }
        }
        return res;
    }
};

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n,m; 
    if(!(cin>>n>>m)) return 0;

    PushRelabel pr(n);
    for(int i=0;i<m;++i){
        int a,b; cin>>a>>b; --a; --b;
        pr.add_undirected(a,b,1);  // each street: capacity 1 in both directions
    }

    int s=0, t=n-1;
    long long k = pr.maxflow(s,t);

    auto cut = pr.mincut_edges(s);     // edges crossing S -> T
    cout << cut.size() << '\n';
    // For Police Chase (unit caps), |cut| == k. Print first k just in case.
    int printed=0;
    for (auto [u,v,c] : cut){
        if (printed==k) break;
        cout << u+1 << ' ' << v+1 << '\n';
        ++printed;
    }
    return 0;
}
