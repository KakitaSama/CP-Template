/*************** DISCRETE LOGARITHM & DISCRETE ROOT TOOLS ***************/
using ll = long long;

/*
 * Helpers: gcd, extended gcd, modpow, modinv_general
 * (You probably already have these in your mega template; if yes, reuse.)
 */

// gcd
ll gcdll(ll a, ll b) {
    while (b) {
        a %= b;
        swap(a, b);
    }
    return a;
}

// extended gcd: ax + by = g
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

// fast exponentiation mod m: a^e % m (works for any m)
ll modpow_ll(ll a, long long e, ll m) {
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

// modular inverse: a^{-1} mod m, assuming gcd(a,m)=1
// (If no inverse exists, returns -1.)
ll modinv_general(ll a, ll m) {
    ll x, y;
    ll g = extgcd(a, m, x, y);
    if (g != 1 && g != -1) return -1;
    x %= m;
    if (x < 0) x += m;
    return x;
}

/*
 * ------------------------ DISCRETE LOGARITHM ------------------------
 *
 * Solve a^x ≡ b (mod m) for integer x >= 0 (if it exists).
 * Returns:
 *   - smallest non-negative x if a solution exists
 *   - -1 if no solution
 *
 * Algorithm: Baby-Step Giant-Step (BSGS) with gcd handling.
 *
 * Idea (coprime case gcd(a,m)=1):
 *   - Let n = ceil(sqrt(m)).
 *   - Precompute "baby steps": a^0, a^1, ..., a^{n-1} and store (value -> exponent j).
 *   - Let factor = a^{-n} (mod m).
 *   - Then for i from 0..n:
 *       we want a^{i*n + j} = b
 *       rearrange: a^j = b * (a^{-n})^i
 *       so we check if cur = b * factor^i appears in the baby-step table.
 *
 * For gcd(a,m) != 1, we repeatedly divide out gcd(a,m) while adjusting b and m.
 * Reference idea: CP-algorithms style "discrete logarithm (modulo not necessarily prime)".
 *
 * Time complexity: O( sqrt(m) ) memory and operations.
 *   m up to ~1e9 is OK (sqrt ~ 3e4), bigger mod can be slow.
 */
ll discrete_log(ll a, ll b, ll m) {
    a %= m; if (a < 0) a += m;
    b %= m; if (b < 0) b += m;

    if (m == 1) return 0;
    if (b == 1 % m) return 0; // x=0 is a trivial solution

    ll cnt = 0;   // how many times we factor out gcd(a,m)
    ll cur = 1;   // keeps track of a^cnt / gcd factors
    ll g;

    // Reduce the equation while gcd(a,m) > 1:
    // if g = gcd(a,m), then:
    //   a^x ≡ b (mod m)
    // => (a/g)^x ≡ b/g (mod m/g), but only if b % g == 0
    while ((g = gcdll(a, m)) > 1) {
        if (b % g != 0) {
            // No solution if b is not divisible by gcd
            return -1;
        }
        m /= g;
        b /= g;
        cur = (cur * (a / g)) % m;
        cnt++;
        if (cur == b) {
            // We found a solution with x = cnt
            return cnt;
        }
    }

    // Now gcd(a, m) == 1
    // Use standard BSGS for a^x = b * cur^{-1} (mod m)
    ll inv_cur = modinv_general(cur, m);
    if (inv_cur == -1) return -1; // shouldn't happen if gcd(a,m)=1 and we managed cur

    b = (b * inv_cur) % m; // we solve a^x = b (mod m) now

    // Baby-step giant-step
    ll n = (ll) sqrt((long double)m) + 1;

    // Baby steps: table[a^j] = j for j in [0, n-1]
    unordered_map<ll,ll> table;
    table.reserve((size_t)n * 2);
    table.max_load_factor(0.7f);

    ll value = 1;
    for (ll j = 0; j < n; ++j) {
        // Only store first occurrence of each value (minimal exponent)
        if (!table.count(value))
            table[value] = j;
        value = (__int128)value * a % m;
    }

    // Giant step factor = a^{-n} mod m
    ll a_inv = modinv_general(a, m);
    if (a_inv == -1) return -1; // no inverse, no solution in coprime case
    ll factor = modpow_ll(a_inv, n, m);

    // Giant steps: b, b*factor, b*factor^2, ...
    ll gamma = b;
    for (ll i = 0; i <= n; ++i) {
        // We want: a^j = gamma
        auto it = table.find(gamma);
        if (it != table.end()) {
            // Found: x = i*n + j  (plus the cnt offsets we removed at the beginning)
            ll x = i * n + it->second + cnt;
            return x;
        }
        gamma = (__int128)gamma * factor % m;
    }

    return -1; // no solution
}

/*
 * ------------------------ PRIMITIVE ROOT (mod prime) ------------------------
 *
 * Find a primitive root g modulo prime p.
 *   - g is a generator of the multiplicative group (Z/pZ)^*, which has size phi(p)=p-1.
 *   - This means that g^k produces all non-zero residues mod p as k runs from 0..p-2.
 *
 * How:
 *   - Factor p-1 = ∏ q_i^{e_i}.
 *   - For candidate g = 2,3,4,...:
 *       check that for all prime factors q_i:
 *          g^{(p-1)/q_i} != 1 (mod p).
 *     The first that passes is a primitive root.
 *
 * Time: O( sqrt(p-1) * log p ) for factoring p-1 (using trial division or something better),
 *       plus some extra log p for power checks.
 */
vector<ll> factorize_ll_simple(ll n) {
    vector<ll> fac;
    for (ll i = 2; i * i <= n; ++i) {
        if (n % i == 0) {
            fac.push_back(i);
            while (n % i == 0) n /= i;
        }
    }
    if (n > 1) fac.push_back(n);
    return fac;
}

ll primitive_root_prime(ll p) {
    // p must be prime
    ll phi = p - 1;
    vector<ll> fac = factorize_ll_simple(phi);
    for (ll g = 2; g < p; ++g) {
        bool ok = true;
        for (ll q : fac) {
            if (modpow_ll(g, phi / q, p) == 1) {
                ok = false;
                break;
            }
        }
        if (ok) return g;
    }
    return -1; // should never happen for prime p
}

/*
 * ------------------------ DISCRETE ROOT (mod prime) ------------------------
 *
 * Solve x^k ≡ a (mod p) where:
 *   - p is prime,
 *   - k >= 1,
 *   - 0 <= a < p.
 *
 * Idea:
 *   - If a == 0, solution is x = 0 (only if k > 0).
 *   - Let g be a primitive root mod p, so every non-zero x can be written x = g^t.
 *       x^k ≡ a  =>  (g^t)^k ≡ a  =>  g^{t*k} ≡ a.
 *   - Let y be such that g^y ≡ a (discrete log of a base g).
 *       Then we need t*k ≡ y (mod p-1)   [since g has order p-1].
 *   - Solve k * t ≡ y (mod p-1) using linear congruence:
 *       k * t - y = (p-1) * m.
 *   - If gcd(k, p-1) doesn't divide y, no solution.
 *   - Otherwise we get one solution t0; all solutions differ by (p-1)/gcd.
 *   - Any solution x = g^t0 is a k-th root of a mod p.
 *
 * This function returns ONE solution x if exists, or -1 if none.
 *
 * Complexity:
 *   - primitive_root_prime: O(sqrt(p-1) * log p)
 *   - discrete_log via BSGS: O(sqrt(p) * log p)
 *   So overall roughly O(sqrt(p)) for one discrete root.
 */

ll discrete_root_prime(ll k, ll a, ll p) {
    a %= p;
    if (a < 0) a += p;
    if (a == 0) {
        // x^k = 0 has solution x=0 for any k >= 1
        return 0;
    }

    // 1) Get primitive root g of mod p
    ll g = primitive_root_prime(p);
    if (g == -1) return -1; // shouldn't happen for prime p

    // 2) Write a = g^y mod p => y = discrete_log(g, a, p)
    ll y = discrete_log(g, a, p);
    if (y == -1) return -1; // should not happen if g is generator and a != 0

    // Solve k * t ≡ y (mod p-1)
    ll phi = p - 1;
    ll x0, y0;
    ll d = extgcd(k, phi, x0, y0); // x0*k + y0*phi = d = gcd(k, phi)

    if (y % d != 0) {
        // No solution if gcd(k, p-1) doesn't divide y
        return -1;
    }

    // Particular solution:
    // t ≡ x0 * (y/d)  (mod phi/d)
    ll t_mod = phi / d;
    ll t = (x0 % t_mod + t_mod) % t_mod; // make x0 positive mod t_mod
    t = (__int128)t * (y / d) % t_mod;

    // x = g^t mod p is ONE solution
    ll x = modpow_ll(g, t, p);
    return x;
}

/*************** END DISCRETE LOG & ROOT TOOLS ***************/
