#include <bits/stdc++.h>
using namespace std;

// ========================= MOD & BASIC OPS =========================
// We work in mod 998244353, a standard NTT-friendly prime: MOD = c*2^23 + 1.
// G is a primitive root modulo MOD.
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

// Fast exponentiation: O(log e)
int mod_pow(long long a, long long e) {
    long long r = 1;
    while (e > 0) {
        if (e & 1) r = (r * a) % MOD;
        a = (a * a) % MOD;
        e >>= 1;
    }
    return int(r);
}

// Fermat inverse: x^(MOD-2) mod MOD. O(log MOD)
int mod_inv(int x) {
    return mod_pow(x, MOD - 2);
}

// ========================= NTT CORE =========================
// NTT (Number Theoretic Transform) = FFT over finite field.
// Complexity for transform of size N: O(N log N).
//
// We precompute roots of unity for all powers of two up to needed size.
// ensure_base grows the global base (roots, iroots) lazily.
//
// convolution(A, B) -> C = A * B (polynomial multiplication) in O((n1+n2) log (n1+n2)).
int ntt_base = 1;
vector<int> roots = {0, 1};
vector<int> iroots = {0, 1};
vector<int> rev;

void ensure_base(int nbase) {
    if (nbase <= ntt_base) return;
    roots.resize(1 << nbase);
    iroots.resize(1 << nbase);

    while (ntt_base < nbase) {
        // z is primitive 2^(base+1)-th root of unity
        int z  = mod_pow(G, (MOD - 1) >> (ntt_base + 1));
        int iz = mod_inv(z);
        int len = 1 << (ntt_base - 1);
        // From smaller base, fill double-sized base:
        // roots[2*i]   = roots[i]
        // roots[2*i+1] = roots[i] * z
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

    // Precompute bit-reversed indices for this n
    if ((int)rev.size() != n) {
        rev.assign(n, 0);
        for (int i = 1; i < n; i++) {
            rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (lg - 1));
        }
    }

    // Apply bit-reversal permutation
    for (int i = 0; i < n; i++) {
        if (i < rev[i]) swap(a[i], a[rev[i]]);
    }

    // Iterative Cooley–Tukey
    for (int len = 1; len < n; len <<= 1) {
        for (int i = 0; i < n; i += (len << 1)) {
            for (int j = 0; j < len; j++) {
                int u = a[i + j];
                // root index len+j corresponds to w^(j) for block size 2*len
                int v = mulmod(a[i + j + len],
                               invert ? iroots[len + j] : roots[len + j]);
                a[i + j]       = addmod(u, v);
                a[i + j + len] = submod(u, v);
            }
        }
    }

    // For inverse transform: divide by n
    if (invert) {
        int inv_n = mod_inv(n);
        for (int &x : a)
            x = mulmod(x, inv_n);
    }
}

// Polynomial multiplication via NTT; A, B are coefficient arrays.
// Complexity: O((|A| + |B|) log (|A| + |B|))
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
// Represent polynomial as Poly = vector<int>, where index = power of x.
using Poly = vector<int>;

// Remove leading zeros to keep degree minimal.
void normalize(Poly &a) {
    while (!a.empty() && a.back() == 0) a.pop_back();
}

// Pointwise add: O(n)
Poly operator+(Poly a, const Poly &b) {
    if (b.size() > a.size()) a.resize(b.size());
    for (size_t i = 0; i < b.size(); i++)
        a[i] = addmod(a[i], b[i]);
    normalize(a);
    return a;
}

// Pointwise subtract: O(n)
Poly operator-(Poly a, const Poly &b) {
    if (b.size() > a.size()) a.resize(b.size());
    for (size_t i = 0; i < b.size(); i++)
        a[i] = submod(a[i], b[i]);
    normalize(a);
    return a;
}

// Polynomial multiplication: O(n log n) via NTT
Poly operator*(const Poly &a, const Poly &b) {
    return convolution(a, b);
}

// Scalar multiply: O(n)
Poly operator*(Poly a, int k) {
    for (int &x : a)
        x = mulmod(x, k);
    normalize(a);
    return a;
}

// Shift: x^k * f(x) : just prepend k zeros. O(n)
Poly shift(const Poly &a, int k) {
    if (a.empty()) return {};
    Poly res(a.size() + k);
    for (size_t i = 0; i < a.size(); i++)
        res[i + k] = a[i];
    return res;
}

// ========================= DERIVATIVE & INTEGRAL =========================
// deriv(f)(x) = sum_{i>=1} i * f[i] * x^{i-1}  (O(n))
Poly deriv(const Poly &a) {
    if (a.empty()) return {};
    Poly res(max(0, (int)a.size() - 1));
    for (size_t i = 1; i < a.size(); i++)
        res[i - 1] = mulmod((int)i, a[i]);
    return res;
}

// integ(f)(x) = C + sum_{i>=0} f[i] * x^{i+1} / (i+1), here C=0.  (O(n))
Poly integ(const Poly &a) {
    Poly res(a.size() + 1);
    res[0] = 0;
    for (size_t i = 0; i < a.size(); i++) {
        res[i + 1] = mulmod(a[i], mod_inv((int)(i + 1)));
    }
    return res;
}

// ========================= FPS INVERSE =========================
// Compute g s.t. f * g ≡ 1 (mod x^n), assuming f[0] != 0.
// Uses Newton iteration: g_{k+1} = g_k * (2 - f*g_k) mod x^{2^k}.
// Complexity: O(M(n)) = about O(n log n).
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
            // g_new = g * (2 - f*g) = 2g - g*(f*g)
            g[i] = submod(addmod(g[i], g[i]), val);
        }
    }
    g.resize(n);
    return g;
}

// ========================= FPS LOG & EXP =========================
// log_series(f) = log f, assuming f[0] = 1.
// Uses f'/f then integrate. Complexity: dominated by inv => O(M(n)).
Poly log_series(const Poly &f, int n) {
    Poly df = deriv(f);
    Poly inv_f = inv(f, n);
    Poly res = df * inv_f;
    res.resize(n - 1);
    res = integ(res);
    res.resize(n);
    return res;
}

// exp_series(f) = exp f, assuming f[0] = 0.
// Uses Newton iteration: if g approximates exp(f), then log(g) ~ f.
// Complexity: O(M(n)).
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
        rhs[0] = addmod(rhs[0], 1); // 1 + f - log(g)
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
// pow_series(f, k, n) = f^k mod x^n.
// Steps:
//  1) factor out leading zeros: f = x^m * g, g(0) != 0
//  2) write g = c0 * h, with h(0) = 1
//  3) compute h^k via exp(k*log(h))
//  4) multiply back c0^k, shift by m*k.
// Complexity: O(M(n)).
Poly pow_series(Poly f, long long k, int n) {
    normalize(f);
    if (k == 0) {
        Poly res(1);
        res[0] = 1;
        return res;
    }
    if (f.empty()) return Poly(n);

    // Handle leading zeros: f(x) = x^m * g(x), g(0) != 0
    int m0 = 0;
    while (m0 < (int)f.size() && f[m0] == 0) m0++;
    if (m0 > 0) {
        long long shift_deg = 1LL * m0 * k;
        if (shift_deg >= n) return Poly(n); // all zeros
        f.erase(f.begin(), f.begin() + m0);
        normalize(f);
    }

    int c0 = f[0];
    int inv_c0 = mod_inv(c0);
    for (int &x : f)
        x = mulmod(x, inv_c0);

    Poly h = log_series(f, n);
    int kmod = (int)(k % (MOD - 1)); // exponent in log/exp modulo MOD-1
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
// We want coefficient [x^n] P(x) / Q(x), with deg P < deg Q, Q[0] != 0.
// Trick:
//  - Consider Q(x) and Q(-x). Multiply to eliminate either even or odd positions.
//  - Each step halves n by picking only even or only odd coefficients.
// Complexity: O(M(k) log n), where k = deg Q.
Poly negate_x(const Poly &Q) {
    Poly R = Q;
    // Q(-x) = sum q_i (-1)^i x^i: negate coefficients with odd index
    for (int i = 1; i < (int)R.size(); i += 2) {
        R[i] = (MOD - R[i]) % MOD;
    }
    return R;
}

int bostan_mori(Poly P, Poly Q, long long n) {
    normalize(P);
    normalize(Q);
    while (n > 0) {
        Poly Qminus = negate_x(Q);
        // P2(x) = P(x) * Q(-x)
        Poly P2 = P * Qminus;
        P2.resize(Q.size() * 2 - 1);

        // Q2(x) = Q(x) * Q(-x) (even polynomial)
        Poly Q2 = Q * Qminus;
        Q2.resize(Q.size() * 2 - 1);

        // If n is even, we want even coefficients; if odd, odd ones
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
    return mulmod(P[0], invQ0); // coefficient [x^0] P/Q
}

// ========================= KITAMASA (nth term of recurrence) =========================
// Recurrence: a_n = c[0]*a_{n-1} + c[1]*a_{n-2} + ... + c[k-1]*a_{n-k}.
// Input:
//   c    : size k, recurrence coefficients
//   init : size k, initial terms a_0..a_{k-1}
//   n    : query index (0-based)
// Kitamasa interprets recurrence as polynomial modulo P(x) = x^k - c[0]x^{k-1}-...-c[k-1].
// We compute x^n mod P(x) in basis {1, x, ..., x^{k-1}} by repeated squaring.
// Then a_n = sum pol[i] * a_i.
//
// Complexity: O(k^2 log n).

// Multiply two polynomials of degree < k and then reduce by recurrence.
vector<int> kitamasa_combine(const vector<int> &a, const vector<int> &b, const vector<int> &c) {
    int k = (int)c.size();
    vector<int> tmp(2 * k);
    // tmp = a(x) * b(x)
    for (int i = 0; i < k; i++) {
        if (!a[i]) continue;
        long long ai = a[i];
        for (int j = 0; j < k; j++) {
            tmp[i + j] = addmod(tmp[i + j], mulmod(ai, b[j]));
        }
    }
    // Reduce high-degree terms using x^k = c[0]x^{k-1} + ... + c[k-1]
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
    // pol represents x^0 initially
    vector<int> pol(k), base(k);
    pol[0] = 1;        // x^0

    // base represents x^1 in the basis (mod recurrence).
    // For k >= 2, x^1 corresponds to polynomial "x" => [0,1,0,...]
    if (k > 1) base[1] = 1;
    else base[0] = c[0]; // deg=1 edge case

    long long e = n;
    // Binary exponentiation of x^n mod characteristic polynomial.
    while (e > 0) {
        if (e & 1LL) {
            pol = kitamasa_combine(pol, base, c);
        }
        base = kitamasa_combine(base, base, c);
        e >>= 1LL;
    }
    // a_n = sum pol[i] * a_i
    int ans = 0;
    for (int i = 0; i < k; i++) {
        ans = addmod(ans, mulmod(pol[i], init[i]));
    }
    return ans;
}

// ========================= END MEGA TEMPLATE =========================

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Example: coefficient of x^n in 1/(1 + x + x^2) using Bostan–Mori:
    //
    // Poly P = {1};
    // Poly Q = {1, 1, 1};
    // long long n; cin >> n;
    // cout << bostan_mori(P, Q, n) << "\n";
    //
    // Example: same sequence via recurrence + Kitamasa:
    //
    // a_n + a_{n-1} + a_{n-2} = 0 => a_n = -a_{n-1} - a_{n-2}
    // c[0] = -1, c[1] = -1 (mod)
    // a_0 = 1, a_1 = -1
    //
    // vector<int> c = {MOD - 1, MOD - 1};
    // vector<int> init = {1, MOD - 1};
    // long long n; cin >> n;
    // cout << kitamasa(c, init, n) << "\n";

    return 0;
}
