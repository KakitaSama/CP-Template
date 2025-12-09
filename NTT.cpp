// ---------- Fast & safe NTT (mod 998244353, G = 3) ----------
const int MOD = 998244353;
const int G   = 3;

inline int addmod(int a, int b) {
    a += b;
    if (a >= MOD) a -= MOD;
    return a;
}

inline int submod(int a, int b) {
    a -= b;
    if (a < 0) a += MOD;
    return a;
}

inline int mulmod(long long a, long long b) {
    return int((a * b) % MOD);
}

int mod_pow(long long a, long long e) {
    long long r = 1;
    while (e > 0) {
        if (e & 1) r = (r * a) % MOD;
        a = (a * a) % MOD;
        e >>= 1;
    }
    return int(r);
}

int mod_inv(int x) {
    return mod_pow(x, MOD - 2);
}

// internals
int ntt_base = 1;
vector<int> roots = {0, 1};
vector<int> iroots = {0, 1};
vector<int> rev;

void ensure_base(int nbase) {
    if (nbase <= ntt_base) {
        // rev might still be too small for this n, handle that in ntt()
        return;
    }

    roots.resize(1 << nbase);
    iroots.resize(1 << nbase);

    while (ntt_base < nbase) {
        int z  = mod_pow(G, (MOD - 1) >> (ntt_base + 1)); // primitive 2^(base+1)-th root
        int iz = mod_inv(z);
        int len = 1 << (ntt_base - 1);
        // fill second half from first half
        for (int i = len; i < (1 << ntt_base); i++) {
            roots[2 * i]      = roots[i];
            roots[2 * i + 1]  = mulmod(roots[i], z);
            iroots[2 * i]     = iroots[i];
            iroots[2 * i + 1] = mulmod(iroots[i], iz);
        }
        ntt_base++;
    }
}

void ntt(vector<int> &a, bool invert) {
    int n = (int)a.size();
    if (n == 1) return;  // <-- avoids using rev/roots when not needed

    int lg = 0;
    while ((1 << lg) < n) lg++;
    ensure_base(lg);

    // build rev for this n if needed
    if ((int)rev.size() != n) {
        rev.assign(n, 0);
        for (int i = 1; i < n; i++) {
            rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (lg - 1));
        }
    }

    for (int i = 0; i < n; i++) {
        if (i < rev[i]) swap(a[i], a[rev[i]]);
    }

    for (int len = 1; len < n; len <<= 1) {
        for (int i = 0; i < n; i += (len << 1)) {
            for (int j = 0; j < len; j++) {
                int u = a[i + j];
                int v = mulmod(a[i + j + len],
                               invert ? iroots[len + j] : roots[len + j]);
                a[i + j]       = addmod(u, v);
                a[i + j + len] = submod(u, v);
            }
        }
    }

    if (invert) {
        int inv_n = mod_inv(n);
        for (int &x : a)
            x = mulmod(x, inv_n);
    }
}

vector<int> convolution(vector<int> A, vector<int> B) {
    if (A.empty() || B.empty()) return {};
    int n1 = (int)A.size();
    int n2 = (int)B.size();
    int n = 1;
    while (n < n1 + n2 - 1) n <<= 1;
    A.resize(n);
    B.resize(n);

    ntt(A, false);
    ntt(B, false);
    for (int i = 0; i < n; i++)
        A[i] = mulmod(A[i], B[i]);
    ntt(A, true);
    A.resize(n1 + n2 - 1);
    return A;
}
// ---------- end NTT ----------


// Example usage

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(0);
    #ifndef ONLINE_JUDGE
    setIO();
    #endif
    int n,m,s;
    cin>>n>>m;
    vector<int> v(m+1,1);
    vector<int> res(m+1,0);
    for(int i = 1 ; i <= m ; i += 2)
    	res[i] = 1;
    int p = n;
    n--;
    swap(res,v);
    while(n > 0)
    {
    	if(n&1)
    		res = convolution(res,v);
    	v = convolution(v,v);
    	if(v.size() > m+1)
    		v.resize(m+1);
    	if(res.size() > m+1)
    		res.resize(m+1);
    	n >>= 1;
    }
    long long ans = 0;
    for(int i = 0; i <= m ; i++)
    {
    	ans = (ans + res[i]);
    	if(ans >= mod)
    		ans -= mod;
    }
    for(long long i = 1 ; i <= p ; i++)
    {
    	ans = (ans*i)%mod;
    }
    cout<<ans<<endl;
}
