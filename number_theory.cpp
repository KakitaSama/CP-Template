/***************** NUMBER THEORY MEGA TEMPLATE *****************/
#include <bits/stdc++.h>
using namespace std;

using ll  = long long;
using ull = unsigned long long;

/*-------------------------------------------------------------
 * 1. BASIC GCD, EXTENDED GCD, MODULAR ARITHMETIC
 *-----------------------------------------------------------*/

// gcd
ll gcdll(ll a, ll b) {
    while (b) {
        a %= b;
        swap(a, b);
    }
    return a;
}

// extended gcd: ax + by = g = gcd(a, b)
ll extgcd(ll a, ll b, ll &x, ll &y) {
    if (b == 0) {
        x = 1; y = 0;
        return a;
    }
    ll x1, y1;
    ll g = extgcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

// fast pow: a^e
ll powll(ll a, long long e) {
    ll r = 1;
    while (e > 0) {
        if (e & 1) r = r * a;
        a = a * a;
        e >>= 1;
    }
    return r;
}

// fast pow mod: a^e mod m  (works even if m not prime)
ll modpow(ll a, long long e, ll m) {
    a %= m;
    if (a < 0) a += m;
    ll r = 1 % m;
    while (e > 0) {
        if (e & 1) r = (__int128)r * a % m;
        a = (__int128)a * a % m;
        e >>= 1;
    }
    return r;
}

// modular inverse if gcd(a,m)=1 (using extgcd)
ll modinv_general(ll a, ll m) {
    ll x, y;
    ll g = extgcd(a, m, x, y);
    if (g != 1 && g != -1) return -1; // no inverse
    x %= m;
    if (x < 0) x += m;
    return x;
}

// If MOD is prime, we can use Fermat: a^(MOD-2) mod MOD.
// We'll later define a global MOD for convenience.
const int MOD = 1000000007; // change if needed
int mod_add(int a, int b) { a += b; if (a >= MOD) a -= MOD; return a; }
int mod_sub(int a, int b) { a -= b; if (a < 0) a += MOD; return a; }
int mod_mul(long long a, long long b) { return int((a * b) % MOD); }
int mod_pow_int(long long a, long long e) {
    long long r = 1;
    while (e > 0) {
        if (e & 1) r = r * a % MOD;
        a = a * a % MOD;
        e >>= 1;
    }
    return int(r);
}
int mod_inv_int(int a) { return mod_pow_int(a, MOD - 2); } // MOD prime

/*-------------------------------------------------------------
 * 2. SIEVE (BASIC) & LINEAR SIEVE (PRIMES, SPF, PHI, MU)
 *-----------------------------------------------------------*/

// Basic sieve (O(n log log n)) — just primes
struct BasicSieve {
    int n;
    vector<int> primes;
    vector<bool> is_composite;

    BasicSieve(int n_) { init(n_); }

    void init(int n_) {
        n = n_;
        is_composite.assign(n + 1, false);
        primes.clear();
        for (int i = 2; i <= n; i++) {
            if (!is_composite[i]) primes.push_back(i);
            for (int p : primes) {
                if (1LL * p * i > n) break;
                is_composite[p * i] = true;
                if (i % p == 0) break;
            }
        }
    }
};

// Linear sieve: primes, SPF, phi, mu in O(n)
struct LinearSieve {
    int n;
    vector<int> primes;
    vector<int> spf;    // smallest prime factor
    vector<int> phi;    // Euler totient
    vector<int> mu;     // Möbius

    LinearSieve(int n_ = 0) {
        if (n_ > 0) init(n_);
    }

    void init(int n_) {
        n = n_;
        primes.clear();
        spf.assign(n + 1, 0);
        phi.assign(n + 1, 0);
        mu.assign(n + 1, 0);

        phi[1] = 1;
        mu[1] = 1;
        for (int i = 2; i <= n; i++) {
            // if i is prime
            if (!spf[i]) {
                spf[i] = i;
                primes.push_back(i);
                phi[i] = i - 1;
                mu[i] = -1;
            }
            for (int p : primes) {
                if (1LL * p * i > n) break;
                spf[p * i] = p;
                if (i % p == 0) {
                    // p divides i
                    phi[p * i] = phi[i] * p;
                    mu[p * i] = 0;
                    break;
                } else {
                    phi[p * i] = phi[i] * (p - 1);
                    mu[p * i] = -mu[i];
                }
            }
        }
    }

    // factorization using SPF (for n <= this->n)
    vector<pair<int,int>> factorize(int x) const {
        vector<pair<int,int>> res;
        while (x > 1) {
            int p = spf[x];
            int cnt = 0;
            while (x % p == 0) x /= p, cnt++;
            res.push_back({p, cnt});
        }
        return res;
    }

    // list all divisors using factorization
    vector<int> divisors(int x) const {
        auto f = factorize(x);
        vector<int> divs = {1};
        for (auto [p, e] : f) {
            int sz = (int)divs.size();
            ll cur = 1;
            vector<int> add;
            for (int i = 1; i <= e; i++) {
                cur *= p;
                for (int j = 0; j < sz; j++) {
                    add.push_back((int)(divs[j] * cur));
                }
            }
            divs.insert(divs.end(), add.begin(), add.end());
        }
        sort(divs.begin(), divs.end());
        return divs;
    }
};

/*-------------------------------------------------------------
 * 3. COMBINATORICS MOD (FACTORIAL, INV FACT, nCr)
 *   - works for n up to MAXN with MOD prime
 *-----------------------------------------------------------*/

struct CombMod {
    int N;
    vector<int> fact, invfact, inv; // inv[i] = modular inverse of i

    CombMod(int n = 0) {
        if (n > 0) init(n);
    }

    void init(int n) {
        N = n;
        fact.assign(N + 1, 1);
        invfact.assign(N + 1, 1);
        inv.assign(N + 1, 1);

        for (int i = 2; i <= N; i++) {
            fact[i] = mod_mul(fact[i-1], i);
            inv[i] = mod_mul(MOD - MOD / i, inv[MOD % i]); // linear inverse
            invfact[i] = mod_mul(invfact[i-1], inv[i]);
        }
    }

    int C(int n, int k) const {
        if (k < 0 || k > n) return 0;
        return mod_mul(fact[n], mod_mul(invfact[k], invfact[n-k]));
    }

    int P(int n, int k) const {
        if (k < 0 || k > n) return 0;
        return mod_mul(fact[n], invfact[n-k]);
    }
};

/*-------------------------------------------------------------
 * 4. CRT (CHINESE REMAINDER THEOREM)
 *   General case: may handle not-coprime moduli
 *   Solve:
 *      x ≡ r1 (mod m1)
 *      x ≡ r2 (mod m2)
 *   returns (r, m) with x ≡ r (mod m)
 *   or (0, -1) if no solution
 *-----------------------------------------------------------*/

pair<ll,ll> crt_pair(ll r1, ll m1, ll r2, ll m2) {
    // normalize
    r1 %= m1; if (r1 < 0) r1 += m1;
    r2 %= m2; if (r2 < 0) r2 += m2;

    ll x, y;
    ll g = extgcd(m1, m2, x, y); // x*m1 + y*m2 = g

    if ((r2 - r1) % g != 0) {
        return {0, -1}; // no solution
    }

    ll lcm = m1 / g * m2;
    // solve m1 * k ≡ (r2 - r1) (mod m2)
    ll k = (ll)((__int128)(r2 - r1) / g * x % (m2 / g));
    ll res = r1 + m1 * k;
    res %= lcm;
    if (res < 0) res += lcm;
    return {res, lcm};
}

// CRT for many congruences: [(r0,m0), (r1,m1), ...]
pair<ll,ll> crt(const vector<ll> &r, const vector<ll> &m) {
    ll r0 = 0, m0 = 1; // x ≡ 0 (mod 1)
    for (size_t i = 0; i < r.size(); i++) {
        auto [r1, m1] = crt_pair(r0, m0, r[i], m[i]);
        if (m1 == -1) return {0, -1};
        r0 = r1;
        m0 = m1;
    }
    return {r0, m0};
}

/*-------------------------------------------------------------
 * 5. MILLER–RABIN PRIMALITY TEST (64-bit)
 *-----------------------------------------------------------*/

ll mul_mod_64(ll a, ll b, ll m) {
    return (ll)((__int128)a * b % m);
}

ll pow_mod_64(ll a, ll d, ll m) {
    ll r = 1;
    while (d > 0) {
        if (d & 1) r = mul_mod_64(r, a, m);
        a = mul_mod_64(a, a, m);
        d >>= 1;
    }
    return r;
}

// deterministic Miller-Rabin for 64-bit
bool isPrime64(ll n) {
    if (n < 2) return false;
    // small primes
    for (ll p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29}) {
        if (n % p == 0) return n == p;
    }
    ll d = n - 1, s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        s++;
    }

    auto check = [&](ll a) {
        if (a % n == 0) return true;
        ll x = pow_mod_64(a, d, n);
        if (x == 1 || x == n - 1) return true;
        for (ll r = 1; r < s; r++) {
            x = mul_mod_64(x, x, n);
            if (x == n - 1) return true;
        }
        return false;
    };

    // bases that are enough for 64-bit determinism
    for (ll a : {2, 325, 9375, 28178, 450775, 9780504, 1795265022}) {
        if (!check(a)) return false;
    }
    return true;
}

/*-------------------------------------------------------------
 * 6. POLLARD–RHO 64-bit FACTORIZATION
 *-----------------------------------------------------------*/

mt19937_64 rng_chrono((uint64_t)chrono::steady_clock::now().time_since_epoch().count());

ll pollard_f(ll x, ll c, ll mod) {
    return (mul_mod_64(x, x, mod) + c) % mod;
}

ll pollard_rho(ll n) {
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;
    ll c = (ll)uniform_int_distribution<ll>(1, n-1)(rng_chrono);
    ll x = (ll)uniform_int_distribution<ll>(0, n-1)(rng_chrono);
    ll y = x;
    ll d = 1;
    while (d == 1) {
        x = pollard_f(x, c, n);
        y = pollard_f(pollard_f(y, c, n), c, n);
        d = gcdll(llabs(x - y), n);
        if (d == n) return pollard_rho(n);
    }
    return d;
}

// factor n into prime powers (unordered)
void factor64(ll n, map<ll,int> &fac) {
    if (n == 1) return;
    if (isPrime64(n)) {
        fac[n]++;
        return;
    }
    ll d = pollard_rho(n);
    factor64(d, fac);
    factor64(n / d, fac);
}

/*-------------------------------------------------------------
 * EXAMPLE USAGE in main()
 *-----------------------------------------------------------*/

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Example 1: compute gcd, extended gcd
    /*
    ll a = 30, b = 18;
    ll x, y;
    ll g = extgcd(a, b, x, y);
    cerr << "gcd(" << a << "," << b << ")=" << g << "  x=" << x << " y=" << y << "\n";
    */

    // Example 2: linear sieve
    /*
    LinearSieve ls(1000000);
    cerr << "phi[10]=" << ls.phi[10] << ", mu[10]=" << ls.mu[10] << "\n";
    auto fac10 = ls.factorize(360);
    for (auto [p,e] : fac10) cerr << p << "^" << e << " ";
    cerr << "\n";
    */

    // Example 3: combinatorics
    /*
    CombMod C(1000000);
    cerr << "C(5,2) = " << C.C(5,2) << "\n";
    */

    // Example 4: CRT
    /*
    vector<ll> r = {2, 3};
    vector<ll> m = {3, 5};
    auto [R, M] = crt(r, m);
    cerr << "x ≡ 2 (mod 3), x ≡ 3 (mod 5) => x ≡ " << R << " (mod " << M << ")\n";
    */

    // Example 5: primality / factor64
    /*
    ll n = 1'000'000'007LL;
    cerr << n << " prime? " << (isPrime64(n) ? "YES" : "NO") << "\n";
    map<ll,int> fac;
    factor64(1234567891011LL, fac);
    for (auto &kv : fac) cerr << kv.first << "^" << kv.second << " ";
    cerr << "\n";
    */

    return 0;
}
/***************** END NUMBER THEORY MEGA TEMPLATE *****************/
