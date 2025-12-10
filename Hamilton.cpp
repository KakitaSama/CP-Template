#include <bits/stdc++.h>
// ==================== EULER & HAMILTON PACK ====================
namespace Euler {
    enum Type { NONE=0, PATH=1, CIRCUIT=2 };

    // ---------- Undirected Euler (multiedges OK) ----------
    // Input: n, edges as (u,v). Returns {type, vertex-sequence} (len = m+1 if exists).
    pair<Type, vector<int>> undirected(int n, const vector<pair<int,int>>& E){
        int m=E.size(); vector<int> deg(n,0);
        if(m==0){ return {CIRCUIT, {0}}; }
        for (auto [u,v]:E) deg[u]++, deg[v]++;
        // check connectivity on vertices with deg>0
        vector<vector<pair<int,int>>> g(n);
        for(int i=0;i<m;++i){ auto [u,v]=E[i]; g[u].push_back({v,i}); g[v].push_back({u,i}); }
        int s=-1; for(int i=0;i<n;++i) if(deg[i]){ s=i; break; }
        vector<int> st; st.reserve(n); vector<char> vis(n,0);
        if(s!=-1){
            stack<int> S; S.push(s); vis[s]=1;
            while(!S.empty()){ int v=S.top(); S.pop(); st.push_back(v);
                for(auto [to,id]:g[v]) if(!vis[to] && deg[to]) vis[to]=1, S.push(to);
            }
            for(int i=0;i<n;++i) if(deg[i] && !vis[i]) return {NONE,{}};
        }
        int odd=0, stv=s;
        for(int i=0;i<n;++i) if(deg[i]&1) { ++odd; stv=i; }
        if(!(odd==0 || odd==2)) return {NONE,{}};
        if(odd==0) { /* circuit */ } else { /* trail starts at odd */ }
        vector<char> used(m,0);
        vector<int> it(n,0), stackv, path; stackv.push_back(stv);
        while(!stackv.empty()){
            int v=stackv.back();
            while(it[v]< (int)g[v].size() && used[g[v][it[v]].second]) ++it[v];
            if(it[v]==(int)g[v].size()){ path.push_back(v); stackv.pop_back(); }
            else {
                auto [to,id]=g[v][it[v]++]; if(used[id]) continue;
                used[id]=1; stackv.push_back(to);
            }
        }
        if((int)path.size()!=m+1) return {NONE,{}};
        reverse(path.begin(), path.end());
        return {odd?PATH:CIRCUIT, path};
    }

    // ---------- Directed Euler ----------
    // Conditions:
    //  - weakly connected over vertices with deg>0
    //  - CIRCUIT: in[v]==out[v] for all v
    //  - PATH   : exactly one (out-in)=+1 (start) and one (in-out)=+1 (end), others equal
    // Returns {type, vertex-sequence}.
    pair<Type, vector<int>> directed(int n, const vector<pair<int,int>>& E){
        int m=E.size(); if(m==0) return {CIRCUIT,{0}};
        vector<int> in(n,0), out(n,0); for(auto [u,v]:E){ ++out[u]; ++in[v]; }
        // weak connectivity
        vector<vector<int>> ug(n);
        for(auto [u,v]:E){ ug[u].push_back(v); ug[v].push_back(u); }
        int s0=-1; for(int i=0;i<n;++i) if(in[i]+out[i]) { s0=i; break; }
        vector<char> vis(n,0); if(s0!=-1){
            stack<int> S; S.push(s0); vis[s0]=1;
            while(!S.empty()){ int v=S.top(); S.pop();
                for(int to:ug[v]) if(!vis[to] && (in[to]+out[to])) vis[to]=1, S.push(to);
            }
            for(int i=0;i<n;++i) if((in[i]+out[i]) && !vis[i]) return {NONE,{}};
        }
        int s=-1, t=-1;
        for(int v=0; v<n; ++v){
            int d = out[v]-in[v];
            if(d==1){ if(s!=-1) return {NONE,{}}; s=v; }
            else if(d==-1){ if(t!=-1) return {NONE,{}}; t=v; }
            else if(d!=0) return {NONE,{}};
        }
        Type tp = CIRCUIT;
        if((s==-1) ^ (t==-1)) return {NONE,{}};
        if(s!=-1) tp=PATH; else { // circuit
            for(int i=0;i<n;++i) if(out[i]) { s=i; break; }
        }
        // Hierholzer
        struct Edge{int to; int id;};
        vector<vector<Edge>> g(n);
        for(int i=0;i<m;++i){ auto [u,v]=E[i]; g[u].push_back({v,i}); }
        vector<int> it(n,0), st, path; vector<char> used(m,0);
        st.push_back(s);
        while(!st.empty()){
            int v=st.back();
            while(it[v]<(int)g[v].size() && used[g[v][it[v]].id]) ++it[v];
            if(it[v]==(int)g[v].size()){ path.push_back(v); st.pop_back(); }
            else {
                auto e = g[v][it[v]++]; if(used[e.id]) continue;
                used[e.id]=1; st.push_back(e.to);
            }
        }
        if((int)path.size()!=m+1) return {NONE,{}};
        reverse(path.begin(), path.end());
        return {tp, path};
    }
}

namespace Hamilton {
    // -------- Small-n (<=20~22) bitmask DP for Hamiltonian PATH --------
    // g: adjacency list (n x vector<int>), simple graph.
    bool path(const vector<vector<int>>& g, vector<int>& out){
        int n=g.size(); if(n==0) return false;
        vector<vector<char>> dp(1<<n, vector<char>(n,0));
        vector<vector<int>> par(1<<n, vector<int>(n,-1));
        for(int v=0; v<n; ++v) dp[1<<v][v]=1;
        for(int mask=1; mask<(1<<n); ++mask){
            for(int last=0; last<n; ++last) if(dp[mask][last]){
                for(int to: g[last]) if(!(mask>>to & 1)){
                    int m2=mask|1<<to;
                    if(!dp[m2][to]) dp[m2][to]=1, par[m2][to]=last;
                }
            }
        }
        int full=(1<<n)-1, end=-1;
        for(int v=0; v<n; ++v) if(dp[full][v]) { end=v; break; }
        if(end==-1) return false;
        out.clear(); int mask=full, cur=end;
        while(cur!=-1){ out.push_back(cur); int p=par[mask][cur]; mask^=1<<cur; cur=p; }
        reverse(out.begin(), out.end()); return true;
    }

    // -------- Hamiltonian CYCLE via DP (try each start) --------
    bool cycle(const vector<vector<int>>& g, vector<int>& cyc){
        int n=g.size(); if(n==0) return false;
        vector<vector<char>> adj(n, vector<char>(n,0));
        for(int u=0;u<n;++u) for(int v: g[u]) adj[u][v]=1;
        for(int s=0; s<n; ++s){
            vector<vector<char>> dp(1<<n, vector<char>(n,0));
            vector<vector<int>> par(1<<n, vector<int>(n,-1));
            dp[1<<s][s]=1;
            for(int mask=1; mask<(1<<n); ++mask){
                for(int last=0; last<n; ++last) if(dp[mask][last]){
                    for(int to: g[last]) if(!(mask>>to & 1)){
                        int m2=mask|1<<to;
                        if(!dp[m2][to]) dp[m2][to]=1, par[m2][to]=last;
                    }
                }
            }
            int full=(1<<n)-1;
            for(int last=0; last<n; ++last) if(last!=s && dp[full][last] && adj[last][s]){
                cyc.clear(); int mask=full, cur=last;
                cyc.push_back(s); // will append at end
                vector<int> tmp; while(cur!=-1){ tmp.push_back(cur); int p=par[mask][cur]; mask^=1<<cur; cur=p; }
                reverse(tmp.begin(), tmp.end());
                cyc = tmp; cyc.push_back(s); return true;
            }
        }
        return false;
    }

    // -------- Sufficient conditions (quick checks) --------
    // Dirac: simple undirected graph, n>=3, if deg(v) >= n/2 for all v -> has Hamiltonian cycle.
    bool dirac(const vector<vector<int>>& g){
        int n=g.size(); if(n<3) return false;
        for(int i=0;i<n;++i) if((int)g[i].size() < (n+1)/2) return false;
        return true;
    }
    // Ore: for every non-adjacent u!=v, deg(u)+deg(v) >= n -> Hamiltonian cycle.
    bool ore(const vector<vector<int>>& g){
        int n=g.size(); if(n<3) return false;
        vector<vector<char>> adj(n, vector<char>(n,0));
        for(int u=0;u<n;++u) for(int v: g[u]) adj[u][v]=1;
        for(int u=0;u<n;++u) for(int v=u+1; v<n; ++v){
            if(!adj[u][v] && !adj[v][u]){
                if((int)g[u].size() + (int)g[v].size() < n) return false;
            }
        }
        return true;
    }
    // DAG Hamiltonian path existence: topological order where every consecutive pair has an edge.
    // Return path if exists.
    bool dag_path(const vector<vector<int>>& g, vector<int>& path){
        int n=g.size(); vector<int> indeg(n,0);
        for(int u=0;u<n;++u) for(int v: g[u]) indeg[v]++;
        queue<int> q; for(int i=0;i<n;++i) if(!indeg[i]) q.push(i);
        vector<int> topo;
        while(!q.empty()){ int u=q.front(); q.pop(); topo.push_back(u);
            for(int v: g[u]) if(--indeg[v]==0) q.push(v);
        }
        if((int)topo.size()!=n) return false;
        vector<vector<char>> adj(n, vector<char>(n,0));
        for(int u=0;u<n;++u) for(int v: g[u]) adj[u][v]=1;
        for(int i=0;i+1<n;++i) if(!adj[topo[i]][topo[i+1]]) return false;
        path=topo; return true;
    }
}
// ==================== /PACK ====================
auto [tp, path] = Euler::undirected(n, edges); // edges: vector<pair<int,int>>
if (tp==Euler::NONE) { /* no euler trail/circuit */ }
else { /* path has vertices of Euler PATH/CIRCUIT */ }


auto [tp, path] = Euler::directed(n, directed_edges);


vector<int> hp; bool okP = Hamilton::path(g, hp);    // hp is Hamiltonian path
vector<int> hc; bool okC = Hamilton::cycle(g, hc);   // hc is cycle (last == first)


bool diracOK = Hamilton::dirac(g);
bool oreOK   = Hamilton::ore(g);
vector<int> dagHP; bool ok = Hamilton::dag_path(g, dagHP);
