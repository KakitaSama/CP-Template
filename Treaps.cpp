#include <bits/stdc++.h>
using namespace std;
using ll = long long;

struct Treap {
    struct Node {
        ll val, sum;
        int pr, sz;
        bool rev;
        ll la_mul, la_add;      // (optional) affine lazy tags
        Node *l, *r;
        Node(ll v)
            : val(v), sum(v), pr((int)rng()), sz(1),
              rev(false), la_mul(1), la_add(0), l(nullptr), r(nullptr) {}
    };

    Node* root = nullptr;

    // ---- rng ---- (optional: can use rand() instead if time pressed)
    static uint64_t rng_state;
    static inline uint32_t rng() {
        rng_state ^= rng_state << 7;
        rng_state ^= rng_state >> 9;
        return (uint32_t)rng_state;
    }

    // ---- core utils ---- (essential)
    static inline int _sz(Node* t){ return t ? t->sz : 0; }
    static inline ll  _sum(Node* t){ return t ? t->sum : 0LL; }

    // (optional) apply affine
    static inline void _apply_affine(Node* t, ll m, ll b){
        if(!t) return;
        t->val = m * t->val + b;
        t->sum = m * t->sum + b * _sz(t);
        t->la_mul = m * t->la_mul;
        t->la_add = m * t->la_add + b;
    }

    // (optional) reverse range
    static inline void _apply_rev(Node* t){
        if(!t) return;
        t->rev ^= 1;
        swap(t->l, t->r);
    }

    static inline void _push(Node* t){ // essential if you use lazy
        if(!t) return;
        if(t->rev){
            _apply_rev(t->l);
            _apply_rev(t->r);
            t->rev = false;
        }
        if(t->la_mul != 1 || t->la_add != 0){
            _apply_affine(t->l, t->la_mul, t->la_add);
            _apply_affine(t->r, t->la_mul, t->la_add);
            t->la_mul = 1; t->la_add = 0;
        }
    }

    static inline void _pull(Node* t){ // essential
        if(!t) return;
        t->sz  = 1 + _sz(t->l) + _sz(t->r);
        t->sum = t->val + _sum(t->l) + _sum(t->r);
    }

    // ---- split / merge ---- (essential)
    static void _split(Node* t, Node*& a, Node*& b, int k){
        if(!t){ a=b=nullptr; return; }
        _push(t);
        if(k <= _sz(t->l)){
            _split(t->l, a, t->l, k);
            b = t;
            _pull(b);
        }else{
            _split(t->r, t->r, b, k - _sz(t->l) - 1);
            a = t;
            _pull(a);
        }
    }

    static Node* _merge(Node* a, Node* b){
        if(!a||!b) return a? a : b;
        if(a->pr > b->pr){
            _push(a);
            a->r = _merge(a->r, b);
            _pull(a);
            return a;
        }else{
            _push(b);
            b->l = _merge(a, b->l);
            _pull(b);
            return b;
        }
    }

    // ---- public api ----
    int size() const { return _sz(root); }  // optional
    void clear(){ _destroy(root); root=nullptr; } // optional

    void insert(int pos, ll x){ // essential
        Node *a, *b;
        _split(root, a, b, pos);
        root = _merge(_merge(a, new Node(x)), b);
    }

    void erase(int pos){ // essential
        Node *a, *b, *c;
        _split(root, a, b, pos);
        _split(b, b, c, 1);
        _destroy(b);
        root = _merge(a, c);
    }

    ll get(int pos){ // optional
        Node *a, *b, *c;
        _split(root, a, b, pos);
        _split(b, b, c, 1);
        ll res = 0;
        if(b){ _push(b); res = b->val; }
        root = _merge(_merge(a, b), c);
        return res;
    }

    void set(int pos, ll x){ // optional
        Node *a, *b, *c;
        _split(root, a, b, pos);
        _split(b, b, c, 1);
        if(b){ _push(b); b->val = x; _pull(b); }
        else b = new Node(x);
        root = _merge(_merge(a, b), c);
    }

    ll range_sum(int l, int r){ // essential
        if(l>=r) return 0;
        Node *a, *b, *c;
        _split(root, a, b, l);
        _split(b, b, c, r-l);
        ll res = _sum(b);
        root = _merge(_merge(a, b), c);
        return res;
    }

    void range_add(int l, int r, ll add){ // optional
        Node *a, *b, *c;
        _split(root, a, b, l);
        _split(b, b, c, r-l);
        _apply_affine(b, 1, add);
        root = _merge(_merge(a, b), c);
    }

    void range_mul(int l, int r, ll mul){ // optional
        Node *a, *b, *c;
        _split(root, a, b, l);
        _split(b, b, c, r-l);
        _apply_affine(b, mul, 0);
        root = _merge(_merge(a, b), c);
    }

    void range_affine(int l, int r, ll m, ll b){ // optional
        Node *a, *b2, *c;
        _split(root, a, b2, l);
        _split(b2, b2, c, r-l);
        _apply_affine(b2, m, b);
        root = _merge(_merge(a, b2), c);
    }

    void range_reverse(int l, int r){ // optional
        Node *a, *b, *c;
        _split(root, a, b, l);
        _split(b, b, c, r-l);
        _apply_rev(b);
        root = _merge(_merge(a, b), c);
    }

private:
    static void _destroy(Node* t){ // optional
        if(!t) return;
        _destroy(t->l);
        _destroy(t->r);
        delete t;
    }
};

uint64_t Treap::rng_state =
    (uint64_t)chrono::steady_clock::now().time_since_epoch().count();
// ------------------------------ demo ------------------------------
// int main(){
//     ios::sync_with_stdio(false);
//     cin.tie(nullptr);
//
//     int n,q; cin>>n>>q;
//     Treap T;
//     vector<ll> a(n);
//     for(int i=0;i<n;i++) cin>>a[i];
//     T.build(a);
//
//     while(q--){
//         int tp; cin>>tp;
//         if(tp==0){ int i; ll x; cin>>i>>x; T.insert(i,x); }
//         else if(tp==1){ int i; cin>>i; T.erase(i); }
//         else if(tp==2){ int l,r; cin>>l>>r; T.range_reverse(l,r); }
//         else if(tp==3){ int l,r; ll w,b; cin>>l>>r>>w>>b; T.range_affine(l,r,w,b); }
//         else if(tp==4){ int l,r; cin>>l>>r; cout<<T.range_sum(l,r)<<"\n"; }
//         else if(tp==5){ int i; ll x; cin>>i>>x; T.set(i,x); }
//         else if(tp==6){ int i; cout<<T.get(i)<<"\n"; }
//     }
// }
