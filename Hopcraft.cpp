struct HopcroftKarp {
    int nL, nR;
    vector<vector<int>> adj;     // adj[uL] = list of neighbors vR
    vector<int> dist, matchL, matchR; // matchL[u] = matched v or -1; matchR[v] = matched u or -1

    HopcroftKarp(int _nL=0, int _nR=0) { init(_nL,_nR); }
    void init(int _nL,int _nR){
        nL=_nL; nR=_nR;
        adj.assign(nL,{});
        dist.assign(nL,0);
        matchL.assign(nL,-1);
        matchR.assign(nR,-1);
    }
    // Add edge from left uL (0..nL-1) to right vR (0..nR-1)
    void add_edge(int uL,int vR){ adj[uL].push_back(vR); }

    // BFS builds level graph of alternating paths (left-unmatched as sources)
    bool bfs(){
        queue<int> q;
        const int INF=1e9;
        bool foundFree=false;
        for(int u=0;u<nL;++u){
            if(matchL[u]==-1){ dist[u]=0; q.push(u); }
            else dist[u]=INF;
        }
        while(!q.empty()){
            int u=q.front(); q.pop();
            for(int v: adj[u]){
                int u2=matchR[v];
                if(u2==-1) foundFree=true;           // reachable free right -> some augmenting path exists
                else if(dist[u2]==INF){               // follow matching edge R->L
                    dist[u2]=dist[u]+1; q.push(u2);
                }
            }
        }
        return foundFree;
    }

    // DFS tries to find augmenting paths along the BFS layers
    bool dfs(int u){
        for(int v: adj[u]){
            int u2=matchR[v];
            if(u2==-1 || (dist[u2]==dist[u]+1 && dfs(u2))){
                matchL[u]=v; matchR[v]=u; return true;
            }
        }
        dist[u]=INT_MAX; // dead end; prune in future DFS
        return false;
    }

    // Run Hopcroft–Karp to compute maximum matching size
    int max_matching(){
        int m=0;
        while(bfs()){
            for(int u=0;u<nL;++u)
                if(matchL[u]==-1 && dfs(u)) ++m;
        }
        return m;
    }

    // --- Helpers to extract MVC / MIS (Kőnig) ---
    // Minimum Vertex Cover construction after max_matching():
    // Do alternating BFS from ALL unmatched left vertices along:
    //   - non-matching edges L->R
    //   - matching edges R->L
    // Then MVC = (L \ VisL) U (R ∩ VisR)
    pair<vector<int>,vector<int>> min_vertex_cover() const {
        vector<char> visL(nL,0), visR(nR,0);
        queue<int> q;
        for(int u=0;u<nL;++u) if(matchL[u]==-1) visL[u]=1, q.push(u);
        while(!q.empty()){
            int u=q.front(); q.pop();
            for(int v: adj[u]){
                if(matchL[u]!=v && !visR[v]){          // follow non-matching L->R
                    visR[v]=1;
                    int u2=matchR[v];
                    if(u2!=-1 && !visL[u2]) visL[u2]=1, q.push(u2); // follow matching R->L
                }
            }
        }
        vector<int> coverL, coverR;
        for(int u=0;u<nL;++u) if(!visL[u]) coverL.push_back(u);
        for(int v=0;v<nR;++v) if( visR[v]) coverR.push_back(v);
        return {coverL, coverR};
    }

    // Maximum Independent Set = complement of MVC in bipartite graphs.
    // Using the same alternating BFS visitation:
    // MIS = (L ∩ VisL) U (R \ VisR)
    pair<vector<int>,vector<int>> max_independent_set() const {
        vector<char> visL(nL,0), visR(nR,0);
        queue<int> q;
        for(int u=0;u<nL;++u) if(matchL[u]==-1) visL[u]=1, q.push(u);
        while(!q.empty()){
            int u=q.front(); q.pop();
            for(int v: adj[u]){
                if(matchL[u]!=v && !visR[v]){
                    visR[v]=1;
                    int u2=matchR[v];
                    if(u2!=-1 && !visL[u2]) visL[u2]=1, q.push(u2);
                }
            }
        }
        vector<int> indepL, indepR;
        for(int u=0;u<nL;++u) if( visL[u]) indepL.push_back(u);
        for(int v=0;v<nR;++v) if(!visR[v]) indepR.push_back(v);
        return {indepL, indepR};
    }
};

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // INPUT FORMAT for CSES 1696 "School Dance":
    // nL nR k
    // then k lines of (boy, girl) 1-indexed pairs meaning an edge from boy->girl
    int nL, nR, k;
    if(!(cin >> nL >> nR >> k)) return 0;

    HopcroftKarp hk(nL, nR);

    // Read edges and convert to 0-index for the template
    for(int i=0;i<k;++i){
        int a,b; cin >> a >> b; // boys:1..nL, girls:1..nR
        hk.add_edge(a-1, b-1);
    }

    // 1) Compute maximum matching
    int M = hk.max_matching();

    // 2) Standard CSES output: print matching size and the pairs
    cout << M << "\n";
    for(int u=0; u<nL; ++u)
        if(hk.matchL[u] != -1)
            cout << (u+1) << " " << (hk.matchL[u]+1) << "\n"; // back to 1-index

    // 3) (Optional) If you want to SEE MVC / MIS while testing locally:
    //    compile with:  g++ -O2 -DPRINT_MVC_MIS file.cpp
#ifdef PRINT_MVC_MIS
    auto [coverL, coverR] = hk.min_vertex_cover();
    auto [indepL, indepR] = hk.max_independent_set();

    cerr << "MVC size = " << M << "\n"; // equals matching size
    cerr << "MVC L (" << coverL.size() << "):";
    for(int u: coverL) cerr << " " << (u+1);
    cerr << "\nMVC R (" << coverR.size() << "):";
    for(int v: coverR) cerr << " " << (v+1);
    cerr << "\n";

    cerr << "MIS size = " << (nL + nR - M) << "\n";
    cerr << "MIS L (" << indepL.size() << "):";
    for(int u: indepL) cerr << " " << (u+1);
    cerr << "\nMIS R (" << indepR.size() << "):";
    for(int v: indepR) cerr << " " << (v+1);
    cerr << "\n";
#endif
    return 0;
}
