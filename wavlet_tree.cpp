// Wavelet Tree (static). Values are ints. Positions: 1..n inside WT (wrap for 0-index arrays).
struct Wavelet {
    int lo, hi; Wavelet *L=nullptr, *R=nullptr;
    vector<int> b;              // b[i] = #elements sent to left in prefix [1..i]
    vector<long long> ps;       // ps[i] = sum of elements sent left in prefix
    // Build on [from,to), values in [x,y]
    Wavelet(vector<int>::iterator from, vector<int>::iterator to, int x, int y): lo(x), hi(y){
        if(from>=to || lo==hi){ b={0}; ps={0}; return; }
        int mid = lo + ((hi-lo)>>1);
        auto f = [mid](int v){ return v<=mid; };
        int n = (int)distance(from,to);
        b.resize(n+1); ps.resize(n+1);
        vector<int> Ls; Ls.reserve(n); vector<int> Rs; Rs.reserve(n);
        int cnt=0; long long s=0;
        for(int i=0;i<n;++i){
            int v = *(from+i);
            if(f(v)){ Ls.push_back(v); ++cnt; s+=v; } else Rs.push_back(v);
            b[i+1]=cnt; ps[i+1]=s;
        }
        if(!Ls.empty()) L = new Wavelet(Ls.begin(), Ls.end(), lo, mid);
        if(!Rs.empty()) R = new Wavelet(Rs.begin(), Rs.end(), mid+1, hi);
    }
    // Convenience ctor from vector<int> a (0-index input → WT uses 1..n)
    Wavelet(vector<int> a): Wavelet(a.begin(), a.end(),
           *min_element(a.begin(),a.end()), *max_element(a.begin(),a.end())) {}

    ~Wavelet(){ delete L; delete R; }

    // kth smallest in [l,r], 1-indexed positions
    int kth(int l,int r,int k){
        if(l>r) return INT_MIN;
        if(lo==hi) return lo;
        int inL = b[r]-b[l-1];
        if(k<=inL) return L? L->kth(b[l-1]+1, b[r], k) : INT_MIN;
        return R? R->kth(l - b[l-1], r - b[r], k - inL) : INT_MIN;
    }
    // count <= x in [l,r]
    int lte(int l,int r,int x){
        if(l>r || x<lo) return 0;
        if(hi<=x) return r-l+1;
        int left = L? L->lte(b[l-1]+1, b[r], x) : 0;
        int right= R? R->lte(l - b[l-1], r - b[r], x) : 0;
        return left + right;
    }
    // count == x in [l,r]
    int cnt(int l,int r,int x){
        if(l>r || x<lo || x>hi) return 0;
        if(lo==hi) return r-l+1;
        int mid = lo + ((hi-lo)>>1);
        if(x<=mid) return L? L->cnt(b[l-1]+1, b[r], x) : 0;
        return R? R->cnt(l - b[l-1], r - b[r], x) : 0;
    }
    // count of a<=val<=b in [l,r]
    int freq(int l,int r,int a,int b){ if(a>b) return 0; return lte(l,r,b) - lte(l,r,a-1); }
    // kth largest
    int kthLargest(int l,int r,int k){ return kth(l,r, r-l+1 - k + 1); }

    // sum of k smallest in [l,r]
    long long sumk(int l,int r,int k){
        if(l>r || k<=0) return 0;
        if(lo==hi) return 1LL*k*lo;
        int inL = b[r]-b[l-1];
        long long sL = ps[r]-ps[l-1];
        if(k<=inL) return L? L->sumk(b[l-1]+1, b[r], k) : 0;
        return sL + (R? R->sumk(l - b[l-1], r - b[r], k - inL) : 0);
    }
    // sum of values <= x in [l,r]
    long long sumLTE(int l,int r,int x){
        if(l>r || x<lo) return 0;
        if(hi<=x) return sumk(l,r, r-l+1); // sum of all
        int mid = lo + ((hi-lo)>>1);
        long long left = L? L->sumLTE(b[l-1]+1, b[r], x) : 0;
        if(x<=mid) return left;
        long long sL = ps[r]-ps[l-1];
        long long right = R? R->sumLTE(l - b[l-1], r - b[r], x) : 0;
        return sL + right;
    }
    // predecessor (max value <= x) in [l,r] ; INT_MIN if none
    int pred(int l,int r,int x){ int c = lte(l,r,x); return c? kth(l,r,c) : INT_MIN; }
    // successor (min value >= x) in [l,r] ; INT_MAX if none
    int succ(int l,int r,int x){ int c = lte(l,r,x-1); return (c==r-l+1)? INT_MAX : kth(l,r,c+1); }
};

// given vector<int> a (0-index); we query on positions [1..n]
vector<int> a = /* ... */; int n = a.size();
Wavelet wt(a);

// examples:
int l=2, r=10;                  // 1-indexed positions
int ksmall = wt.kth(l,r,3);     // 3rd smallest in [l,r]
int howMany = wt.lte(l,r,100);  // count ≤ 100 in [l,r]
int cntX = wt.cnt(l,r,42);      // frequency of 42 in [l,r]
int inRange = wt.freq(l,r,10,50);      // count in [10..50]
int klarg = wt.kthLargest(l,r,2);      // 2nd largest
long long sumSmall = wt.sumk(l,r,5);   // sum of 5 smallest
long long sumLeX = wt.sumLTE(l,r,100); // sum of values ≤ 100
int p = wt.pred(l,r,37);        // predecessor of 37
int s = wt.succ(l,r,37);        // successor of 37
