// Segment Tree Beats (sum, max, add, chmin, mod). 0-indexed.
struct SegBeatsMod {
    struct Nd {
        long long sum, mx, smax; int mxc, len; long long add;
        Nd(long long v=0,int L=1): sum(v), mx(v), smax(LLONG_MIN), mxc(1), len(L), add(0) {}
    };
    int n; vector<Nd> st;
    SegBeatsMod(int N=0){init(N);}
    void init(int N){ n=N; st.assign(4*n, Nd()); }
    Nd merge(const Nd& A, const Nd& B){
        Nd C; C.len=A.len+B.len; C.sum=A.sum+B.sum; C.add=0;
        if (A.mx==B.mx){ C.mx=A.mx; C.mxc=A.mxc+B.mxc; C.smax=max(A.smax,B.smax); }
        else if (A.mx>B.mx){ C.mx=A.mx; C.mxc=A.mxc; C.smax=max(A.smax,B.mx); }
        else { C.mx=B.mx; C.mxc=B.mxc; C.smax=max(B.smax,A.mx); }
        return C;
    }
    void build(const vector<long long>& a,int p,int l,int r){
        if(l==r){ st[p]=Nd(a[l],1); return; }
        int m=(l+r)>>1; build(a,p<<1,l,m); build(a,p<<1|1,m+1,r); st[p]=merge(st[p<<1],st[p<<1|1]);
    }
    void build(const vector<long long>& a){ init((int)a.size()); if(n) build(a,1,0,n-1); }

    inline void apply_add(int p,long long d){
        Nd &x=st[p]; x.sum+=d*x.len; x.mx+=d; if(x.smax!=LLONG_MIN) x.smax+=d; x.add+=d;
    }
    inline void apply_chmin(int p,long long x){ // pre: x < st[p].mx and x >= st[p].smax
        Nd &a=st[p]; a.sum += (x - a.mx)*a.mxc; a.mx=x;
    }
    void push(int p){
        Nd &a=st[p]; if(a.len==1) { a.add=0; return; }
        if(a.add){ apply_add(p<<1,a.add); apply_add(p<<1|1,a.add); a.add=0; }
        if(st[p<<1].mx > a.mx) apply_chmin(p<<1, a.mx);
        if(st[p<<1|1].mx > a.mx) apply_chmin(p<<1|1, a.mx);
    }

    void r_chmin(int p,int l,int r,int ql,int qr,long long x){
        if(qr<l||r<ql || st[p].mx<=x) return;
        if(ql<=l&&r<=qr && st[p].smax < x){ apply_chmin(p,x); return; }
        push(p); int m=(l+r)>>1;
        r_chmin(p<<1,l,m,ql,qr,x); r_chmin(p<<1|1,m+1,r,ql,qr,x);
        st[p]=merge(st[p<<1],st[p<<1|1]);
    }
    void r_add(int p,int l,int r,int ql,int qr,long long d){
        if(qr<l||r<ql) return;
        if(ql<=l&&r<=qr){ apply_add(p,d); return; }
        push(p); int m=(l+r)>>1;
        r_add(p<<1,l,m,ql,qr,d); r_add(p<<1|1,m+1,r,ql,qr,d);
        st[p]=merge(st[p<<1],st[p<<1|1]);
    }
    // Range modulo: safe (possibly slower) version; pushes to leaves when needed.
    void r_mod(int p,int l,int r,int ql,int qr,long long x){
        if(qr<l||r<ql || st[p].mx < x) return;
        if(l==r){
            long long v = st[p].mx % x;
            st[p]=Nd(v,1);
            return;
        }
        push(p); int m=(l+r)>>1;
        r_mod(p<<1,l,m,ql,qr,x); r_mod(p<<1|1,m+1,r,ql,qr,x);
        st[p]=merge(st[p<<1],st[p<<1|1]);
    }
    long long q_sum(int p,int l,int r,int ql,int qr){
        if(qr<l||r<ql) return 0;
        if(ql<=l&&r<=qr) return st[p].sum;
        push(p); int m=(l+r)>>1;
        return q_sum(p<<1,l,m,ql,qr)+q_sum(p<<1|1,m+1,r,ql,qr);
    }
    long long q_max(int p,int l,int r,int ql,int qr){
        if(qr<l||r<ql) return LLONG_MIN;
        if(ql<=l&&r<=qr) return st[p].mx;
        push(p); int m=(l+r)>>1;
        return max(q_max(p<<1,l,m,ql,qr), q_max(p<<1|1,m+1,r,ql,qr));
    }

    // public wrappers
    void chmin(int l,int r,long long x){ if(n) r_chmin(1,0,n-1,l,r,x); }
    void add  (int l,int r,long long d){ if(n) r_add  (1,0,n-1,l,r,d); }
    void mod  (int l,int r,long long x){ if(n) r_mod  (1,0,n-1,l,r,x); }
    long long sum(int l,int r){ return n? q_sum(1,0,n-1,l,r):0; }
    long long mx (int l,int r){ return n? q_max(1,0,n-1,l,r):LLONG_MIN; }
};
