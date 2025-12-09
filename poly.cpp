#include <bits/stdc++.h>
using namespace std;

// ========================= MOD & BASIC OPS =========================
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

// ========================= NTT CORE =========================
// Adapted for reuse across many convolutions.
int ntt_base = 1;
vector<int> roots = {0, 1};
vector<int> iroots = {0, 1};
vector<int> rev;

void ensure_base(int nbase) {
    if (nbase <= ntt_base) return;
    roots.resize(1 << nbase);
    iroots.resize(1 << nbase);

    while (ntt_base < nbase) {
        int z  = mod_pow(G, (MOD - 1) >> (ntt_base + 1));
        int iz = mod_inv(z);
        int len = 1 << (ntt_base - 1);
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
    if (n == 1) return;

    int lg = 0;
    while ((1 << lg) < n) lg++;
    ensure_base(lg);

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

// ========================= POLYNOMIAL TYPE & BASIC OPS =========================
using Poly = vector<int>;

void normalize(Poly &a) {
    while (!a.empty() && a.back() == 0) a.pop_back();
}

Poly operator+(Poly a, const Poly &b) {
    if (b.size() > a.size()) a.resize(b.size());
    for (size_t i = 0; i < b.size(); i++)
        a[i] = addmod(a[i], b[i]);
    normalize(a);
    return a;
}

Poly operator-(Poly a, const Poly &b) {
    if (b.size() > a.size()) a.resize(b.size());
    for (size_t i = 0; i < b.size(); i++)
        a[i] = submod(a[i], b[i]);
    normalize(a);
    return a;
}

Poly operator*(const Poly &a, const Poly &b) {
    return convolution(a, b);
}

// scalar multiply
Poly operator*(Poly a, int k) {
    for (int &x : a)
        x = mulmod(x, k);
    normalize(a);
    return a;
}

// x^k * f(x)
Poly shift(const Poly &a, int k) {
    if (a.empty()) return {};
    Poly res(a.size() + k);
    for (size_t i = 0; i < a.size(); i++)
        res[i + k] = a[i];
    return res;
}

// ========================= DERIVATIVE & INTEGRAL =========================
Poly deriv(const Poly &a) {
    if (a.empty()) return {};
    Poly res(max(0, (int)a.size() - 1));
    for (size_t i = 1; i < a.size(); i++)
        res[i - 1] = mulmod((int)i, a[i]);
    return res;
}

Poly integ(const Poly &a) {
    Poly res(a.size() + 1);
    res[0] = 0;
    for (size_t i = 0; i < a.size(); i++) {
        res[i + 1] = mulmod(a[i], mod_inv((int)(i + 1)));
    }
    return res;
}

// ========================= FPS INVERSE =========================
// Compute g s.t. f * g ≡ 1 (mod x^n), assuming f[0] != 0
Poly inv(const Poly &f, int n) {
    Poly g(1);
    g[0] = mod_inv(f[0]);
    int m = 1;
    while (m < n) {
        m <<= 1;
        Poly f_cut(min((int)f.size(), m));
        for (int i = 0; i < (int)f_cut.size(); i++)
            f_cut[i] = f[i];

        Poly gg = g * g;   // g^2
        gg.resize(m);
        Poly t = gg * f_cut; // g^2 * f
        t.resize(m);

        g.resize(m);
        for (int i = 0; i < m; i++) {
            int val = (i < (int)t.size() ? t[i] : 0);
            g[i] = submod(addmod(g[i], g[i]), val);
        }
    }
    g.resize(n);
    return g;
}

// ========================= FPS LOG & EXP =========================
// f[0] must be 1
Poly log_series(const Poly &f, int n) {
    Poly df = deriv(f);
    Poly inv_f = inv(f, n);
    Poly res = df * inv_f;
    res.resize(n - 1);
    res = integ(res);
    res.resize(n);
    return res;
}

// f[0] must be 0
Poly exp_series(const Poly &f, int n) {
    Poly g(1);
    g[0] = 1;
    int m = 1;
    while (m < n) {
        m <<= 1;
        Poly g_cut = g;
        g_cut.resize(m);
        Poly lg = log_series(g_cut, m);
        Poly rhs(m);
        for (int i = 0; i < min(m, (int)f.size()); i++)
            rhs[i] = f[i];
        rhs[0] = addmod(rhs[0], 1);
        for (int i = 0; i < m; i++) {
            if (i < (int)lg.size())
                rhs[i] = submod(rhs[i], lg[i]);
        }
        Poly g2 = g_cut * rhs;
        g2.resize(m);
        g = g2;
    }
    g.resize(n);
    return g;
}

// ========================= FPS POWER =========================
// f^k mod x^n
Poly pow_series(Poly f, long long k, int n) {
    normalize(f);
    if (k == 0) {
        Poly res(1);
        res[0] = 1;
        return res;
    }
    if (f.empty()) return Poly(n);

    // factor f = x^m * g, g(0) != 0
    int m0 = 0;
    while (m0 < (int)f.size() && f[m0] == 0) m0++;
    if (m0 > 0) {
        long long shift_deg = 1LL * m0 * k;
        if (shift_deg >= n) return Poly(n);
        f.erase(f.begin(), f.begin() + m0);
        normalize(f);
    }

    int c0 = f[0];
    int inv_c0 = mod_inv(c0);
    for (int &x : f)
        x = mulmod(x, inv_c0);

    Poly h = log_series(f, n);
    int kmod = (int)(k % (MOD - 1));
    for (int &x : h)
        x = mulmod(x, kmod);

    Poly res = exp_series(h, n);
    int c0k = mod_pow(c0, kmod);
    for (int &x : res)
        x = mulmod(x, c0k);

    if (m0 > 0) {
        long long shift_deg = 1LL * m0 * k;
        if (shift_deg >= n) return Poly(n);
        res = shift(res, (int)shift_deg);
        if ((int)res.size() > n) res.resize(n);
    }
    res.resize(n);
    return res;
}

// ========================= BOSTAN–MORI (nth coeff of P/Q) =========================
// Helper: Q(-x)
Poly negate_x(const Poly &Q) {
    Poly R = Q;
    for (int i = 1; i < (int)R.size(); i += 2) {
        R[i] = (MOD - R[i]) % MOD;
    }
    return R;
}

// nth coefficient of P(x)/Q(x)
int bostan_mori(Poly P, Poly Q, long long n) {
    normalize(P);
    normalize(Q);
    while (n > 0) {
        Poly Qminus = negate_x(Q);
        Poly P2 = P * Qminus;
        P2.resize(Q.size() * 2 - 1);

        Poly Q2 = Q * Qminus;
        Q2.resize(Q.size() * 2 - 1);

        Poly P_new((P2.size() + 1) / 2);
        for (int i = (int)(n & 1); i < (int)P2.size(); i += 2)
            P_new[i >> 1] = P2[i];

        Poly Q_new((Q2.size() + 1) / 2);
        for (int i = 0; i < (int)Q2.size(); i += 2)
            Q_new[i >> 1] = Q2[i];

        P.swap(P_new);
        Q.swap(Q_new);
        normalize(P);
        normalize(Q);

        n >>= 1;
    }
    int invQ0 = mod_inv(Q[0]);
    return mulmod(P[0], invQ0);
}

// ========================= KITAMASA (nth term of linear recurrence) =========================
// Recurrence: a_n = c[0]*a_{n-1} + c[1]*a_{n-2} + ... + c[k-1]*a_{n-k]
// Input:
//   c  : size k, coefficients as above
//   init : size k, init[0]=a_0,...,init[k-1]=a_{k-1}
//   n : index (0-based)
// Returns a_n in O(k^2 log n).
//
// Internally, Kitamasa computes polynomial P such that:
//   a_n = sum_{i=0}^{k-1} P[i] * a_i
// where P is x^n mod characteristic polynomial.

vector<int> kitamasa_combine(const vector<int> &a, const vector<int> &b, const vector<int> &c) {
    int k = (int)c.size();
    vector<int> tmp(2 * k);
    for (int i = 0; i < k; i++) {
        if (!a[i]) continue;
        long long ai = a[i];
        for (int j = 0; j < k; j++) {
            tmp[i + j] = addmod(tmp[i + j], mulmod(ai, b[j]));
        }
    }
    // reduce tmp[>=k] using recurrence
    for (int i = 2 * k - 1; i >= k; i--) {
        if (!tmp[i]) continue;
        int val = tmp[i];
        for (int j = 0; j < k; j++) {
            tmp[i - 1 - j] = addmod(tmp[i - 1 - j], mulmod(val, c[j]));
        }
        tmp[i] = 0;
    }
    tmp.resize(k);
    return tmp;
}

int kitamasa(const vector<int> &c, const vector<int> &init, long long n) {
    int k = (int)c.size();
    if (n < (long long)k) {
        return init[(int)n];
    }
    vector<int> pol(k), base(k);
    pol[0] = 1;        // x^0
    // base represents x^1
    if (k > 1) base[1] = 1;
    else base[0] = c[0]; // deg=1 case, adjust but rarely needed

    long long e = n;
    while (e > 0) {
        if (e & 1LL) {
            pol = kitamasa_combine(pol, base, c);
        }
        base = kitamasa_combine(base, base, c);
        e >>= 1LL;
    }
    int ans = 0;
    for (int i = 0; i < k; i++) {
        ans = addmod(ans, mulmod(pol[i], init[i]));
    }
    return ans;
}

// ========================= END MEGA TEMPLATE =========================


// Example main (empty solve):
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Use the above functions in your solve().
    // Example: simple Fibonacci with Kitamasa:
    //
    // a_n = a_{n-1} + a_{n-2}, a_0 = 0, a_1 = 1
    //
    // vector<int> c = {1, 1};               // c0,c1
    // vector<int> init = {0, 1};           // a0,a1
    // long long n; cin >> n;
    // cout << kitamasa(c, init, n) << "\n";

    return 0;
}
