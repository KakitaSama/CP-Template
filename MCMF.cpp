#include <bits/stdc++.h>
using namespace std;

// -------- Minimum-Cost Max-Flow (compact, 0-indexed) --------
struct MCMF {
    struct E{ int to, rev; long long cap, cost, cap0; }; // cap0 = original forward cap
    int n; vector<vector<E>> g;
    MCMF(int N=0){ init(N); }
    void init(int N){ n=N; g.assign(n,{}); }
    void add_edge(int u,int v,long long cap,long long cost){
        E a{v, (int)g[v].size(), cap,  cost, cap};
        E b{u, (int)g[u].size(), 0,   -cost, 0};
        g[u].push_back(a); g[v].push_back(b);
    }
    pair<long long,long long> min_cost_flow(int s,int t,long long need=(long long)4e18){
        const long long INF = (long long)4e18;
        long long F=0, C=0;
        vector<long long> pot(n,0), dist(n);
        vector<int> pv(n), pe(n);
        while(F<need){
            fill(dist.begin(),dist.end(),INF); dist[s]=0;
            using P=pair<long long,int>;
            priority_queue<P, vector<P>, greater<P>> pq; pq.push({0,s});
            while(!pq.empty()){
                auto [d,u]=pq.top(); pq.pop(); if(d!=dist[u]) continue;
                for(int i=0;i<(int)g[u].size();++i){
                    auto &e=g[u][i]; if(e.cap<=0) continue;
                    long long nd = d + e.cost + pot[u] - pot[e.to];
                    if(nd < dist[e.to]){
                        dist[e.to]=nd; pv[e.to]=u; pe[e.to]=i; pq.push({nd, e.to});
                    }
                }
            }
            if(dist[t]==INF) break;
            for(int i=0;i<n;++i) if(dist[i]<INF) pot[i]+=dist[i];
            long long add = need - F;
            for(int v=t; v!=s; v=pv[v]) add = min(add, g[pv[v]][pe[v]].cap);
            for(int v=t; v!=s; v=pv[v]){
                auto &e=g[pv[v]][pe[v]];
                auto &r=g[e.to][e.rev];
                e.cap -= add; r.cap += add;
            }
            F += add; C += add * pot[t];
        }
        return {F,C};
    }
    // ---- OPTIONAL: list forward edges that carry positive flow (u, v, flow, cost) ----
    // NOTE: call this immediately after min_cost_flow(), BEFORE you modify residuals (e.g., path extraction).
    vector<tuple<int,int,long long,long long>> chosen_edges() const {
        vector<tuple<int,int,long long,long long>> res;
        for(int u=0; u<n; ++u){
            for(auto const &e : g[u]){
                if(e.cap0>0){
                    long long sent = g[e.to][e.rev].cap; // reverse residual equals flow sent
                    if(sent>0) res.emplace_back(u, e.to, sent, e.cost);
                }
            }
        }
        return res;
    }
};

// -------- CSES 1698: Distinct Routes (edge-disjoint paths) --------
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m; 
    if(!(cin >> n >> m)) return 0;
    MCMF f(n);
    for(int i=0;i<m;++i){
        int a,b; cin >> a >> b; --a; --b;
        f.add_edge(a,b,1,0);      // capacity 1, zero cost
    }
    int s = 0, t = n-1;

    auto [flow, cost] = f.min_cost_flow(s, t);
    cout << flow << "\n";

    // OPTIONAL: you can fetch the used arcs here (before path extraction)
    // auto used = f.chosen_edges();  // vector of (u,v,flow,cost)

    // Extract 'flow' distinct paths: follow forward edges with reverse residual > 0.
    for(long long k=0; k<flow; ++k){
        vector<int> par(n, -1), pari(n, -1);
        queue<int> q; q.push(s); par[s] = s;
        while(!q.empty() && par[t]==-1){
            int u = q.front(); q.pop();
            for(int i=0;i<(int)f.g[u].size();++i){
                auto &e = f.g[u][i];
                if(e.cap0==0) continue; // skip artificial reverse arcs
                if(par[e.to]==-1 && f.g[e.to][e.rev].cap > 0){ // edge carried flow
                    par[e.to] = u; pari[e.to] = i; q.push(e.to);
                    if(e.to==t) break;
                }
            }
        }
        vector<int> path;
        for(int v=t; v!=s; v=par[v]){
            path.push_back(v);
            auto &e = f.g[par[v]][pari[v]];
            auto &r = f.g[e.to][e.rev];
            r.cap -= 1;  // consume the unit of flow so the edge won't be reused
            e.cap += 1;
        }
        path.push_back(s);
        reverse(path.begin(), path.end());

        cout << (int)path.size() << "\n";
        for(int i=0;i<(int)path.size();++i){
            if(i) cout << ' ';
            cout << (path[i]+1);
        }
        cout << "\n";
    }
    return 0;
}
