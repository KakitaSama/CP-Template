#include <bits/stdc++.h>
// =============================================================
// [NAME]   : SOS DP + Slope Trick + Alien Trick mega template
// [AUTHOR] : Duukh ICPC Template
// [LANG]   : C++17
// [USAGE]  : Copy-paste and plug into your solution.
// =============================================================

using namespace std;
using ll  = long long;
using ld  = long double;
const ll INFLL = (ll)4e18;

// =============================================================
//                        SOS DP
//    (Sum Over Subsets / Sum Over Supersets on bitmasks)
// =============================================================
// ESSENTIAL IDEA:
//   - We have an array f[mask] (size = 2^K).
//   - "sum over submasks" transform:
//        g[mask] = sum_{sub ⊆ mask} f[sub]
//   - "sum over supersets" transform:
//        g[mask] = sum_{sup ⊇ mask} f[sup]
//   - Implemented in O(K * 2^K) using DP over bits.
// =============================================================
namespace SOS {

    // Build "sum over submasks" in-place
    // After this, f[mask] = sum_{sub ⊆ mask} original_f[sub]
    // Complexity: O(K * 2^K)
    void sum_over_submasks(vector<ll> &f, int K) {
        int N = 1 << K;
        for (int bit = 0; bit < K; ++bit) {
            for (int mask = 0; mask < N; ++mask) {
                if (mask & (1 << bit)) {
                    f[mask] += f[mask ^ (1 << bit)];
                }
            }
        }
    }

    // Build "sum over supersets" in-place
    // After this, f[mask] = sum_{sup ⊇ mask} original_f[sup]
    // Complexity: O(K * 2^K)
    void sum_over_supersets(vector<ll> &f, int K) {
        int N = 1 << K;
        for (int bit = 0; bit < K; ++bit) {
            for (int mask = 0; mask < N; ++mask) {
                if ((mask & (1 << bit)) == 0) {
                    f[mask] += f[mask | (1 << bit)];
                }
            }
        }
    }

    // EXTRA: iterate over all submasks of a mask in O(#submasks)
    //  for (int sub = mask; ; sub = (sub - 1) & mask) { ... ; if (!sub) break; }
    // This is just a reminder pattern.
}

// =============================================================
//                        SLOPE TRICK
// =============================================================
// We maintain a convex piecewise-linear function f(x) whose slope
// changes only by ±1 at integer points. Typical operations:
//
//   - add_abs(a):   f(x) += |x - a|
//   - add_a_minus_x_plus(a): f(x) += max(0, a - x)
//   - add_x_minus_a_plus(a): f(x) += max(0, x - a)
//   - shift(dx):    f(x) becomes f(x - dx)
//   - add_const(c): f(x) += c
//
// The structure keeps:
//   - min_f = minimum value of f(x)
//   - two heaps L (max-heap) and R (min-heap) which store "break points"
//   - addL, addR: lazy shifts applied to all elements in L/R
//
// USAGE:
//   SlopeTrick st;
//   st.add_abs(a1);
//   st.add_abs(a2);
//   ll ans = st.get_min();
//
// Classic DP: many CF/AtCoder convex DP problems.
// =============================================================
struct SlopeTrick {
    // minimal value of f(x)
    ll min_f;
    // left: max-heap of "border positions"
    priority_queue<ll> L;
    // right: min-heap
    priority_queue<ll, vector<ll>, greater<ll>> R;
    // lazy offsets
    ll addL, addR;

    SlopeTrick() {
        min_f = 0;
        addL = addR = 0;
    }

    // ------------- internal helpers ------------- //
    inline void push_L(ll x) { L.push(x - addL); }
    inline void push_R(ll x) { R.push(x - addR); }

    inline ll top_L() { return L.empty() ? (ll)-4e18 : L.top() + addL; }
    inline ll top_R() { return R.empty() ? (ll)+4e18 : R.top() + addR; }

    inline void pop_L() { L.pop(); }
    inline void pop_R() { R.pop(); }

    // ------------- basic operations ------------- //

    // f(x) += c
    void add_const(ll c) { min_f += c; }

    // Shift argument: f(x) ← f(x - dx)
    // This is equivalent to shifting all breakpoints by +dx
    void shift(ll dx) {
        addL += dx;
        addR += dx;
    }

    // ESSENTIAL:
    // f(x) += |x - a|
    // Careful balancing of heaps.
    void add_abs(ll a) {
        if (L.empty() && R.empty()) {
            push_L(a);
            push_R(a);
            return;
        }

        ll l = top_L();
        ll r = top_R();

        if (a < l) {
            min_f += l - a;
            pop_L();
            push_L(a);
            push_R(l);
        } else if (a > r) {
            min_f += a - r;
            pop_R();
            push_R(a);
            push_L(r);
        } else {
            // a is between l and r, no min change
            push_L(a);
            push_R(a);
        }
    }

    // EXTRA:
    // f(x) += max(0, a - x)  (aka (a - x)_+ )
    // This enforces that the minimum tends to be on the right of "a".
    void add_a_minus_x_plus(ll a) {
        // ensure all breakpoints on the left side do not exceed 'a'
        ll l = top_L();
        if (l > a) {
            min_f += l - a;
            pop_L();
            push_L(a);
            push_R(l);
        } else {
            push_R(a);
        }
    }

    // EXTRA:
    // f(x) += max(0, x - a)  (aka (x - a)_+ )
    // This enforces that the minimum tends to be on the left of "a".
    void add_x_minus_a_plus(ll a) {
        ll r = top_R();
        if (r < a) {
            min_f += a - r;
            pop_R();
            push_R(a);
            push_L(r);
        } else {
            push_L(a);
        }
    }

    // Query minimum value of f(x)
    ll get_min() const {
        return min_f;
    }

    // Optional: get an x where min is attained (any value in [top_L, top_R])
    ll argmin() {
        ll l = top_L();
        ll r = top_R();
        if (l == (ll)-4e18 && r == (ll)+4e18) return 0; // undefined
        if (l <= r) return l;  // any point in [l,r] is fine; choose l
        return r;
    }
};

// =============================================================
//                        ALIEN TRICK
// =============================================================
// Alien trick is a parametric optimization technique for problems like:
//
//   "Divide into at most K segments / choose at most K items"
//   with some cost per segment/item. Often we:
//
//   - Introduce λ ("penalty per segment").
//   - Solve a DP that maximizes (original_value - λ * (#segments)).
//   - The DP returns (best_value, used_segments).
//   - By binary searching λ, we force used_segments to equal K.
//   - Then answer = best_value + λ * K.
//
// Generic pattern (maximization):
//
//   DP(λ) -> pair<value, segments_used>
//   if segments_used > K: λ too small  (segments are cheap -> we use many)
//   if segments_used < K: λ too large  (segments expensive -> we use few)
//   binary search λ until segments_used == K or closest.
//
// Below is a "template skeleton" for the pattern.
// You must fill the recurrence & the "gain" depending on the problem.
// =============================================================

// Example structure for DP result.
struct AlienRes {
    ll val;   // best value of (original_value - λ * (#segments))
    int cnt;  // number of segments/items used
};

// PROBLEM-SPECIFIC:
//   Global / captured data used by DP, e.g. arrays, prefix sums, etc.
vector<ll> W;  // example array, fill from main
int N_GLOBAL;
int K_GLOBAL;

// Example "gain" between (j -> i). Customize!
// For example, value of segment (j+1 .. i).
// Here we just set something dummy (like sum of W[j+1..i]).
ll gain_segment(int j, int i, const vector<ll> &pref) {
    // segment value = sum(W[j+1..i]) as an example
    return pref[i] - pref[j];
}

// ESSENTIAL: DP for a fixed λ
AlienRes alien_dp(ll lambda, const vector<ll> &pref) {
    // Here we do a simple O(n^2) DP for demonstration:
    // dp[i] = best value up to i
    // cnt[i] = number of segments used for that best value.
    // Transition:
    //   dp[i] = max over j < i: dp[j] + gain_segment(j, i) - lambda
    //   cnt[i] = cnt[j] + 1
    // In practice you will combine alien trick with Convex Hull, Divide&Conquer, etc.

    int n = N_GLOBAL;
    vector<ll> dp(n + 1, -INFLL);
    vector<int> cnt(n + 1, 0);

    dp[0] = 0;
    cnt[0] = 0;

    for (int i = 1; i <= n; ++i) {
        ll bestVal = -INFLL;
        int bestCnt = 0;
        for (int j = 0; j < i; ++j) {
            ll candidate = dp[j] + gain_segment(j, i, pref) - lambda;
            if (candidate > bestVal) {
                bestVal = candidate;
                bestCnt = cnt[j] + 1;
            }
        }
        dp[i] = bestVal;
        cnt[i] = bestCnt;
    }

    AlienRes res;
    res.val = dp[n];
    res.cnt = cnt[n];
    return res;
}

// Binary search on λ to get exactly K_GLOBAL segments (or closest).
//   lambda_low, lambda_high: search range for penalty (problem-dependent).
//   The function returns the optimal value of the *original* objective.
// NOTE: you MUST set N_GLOBAL, K_GLOBAL, and build pref before using.
ll alien_solve(ll lambda_low, ll lambda_high) {
    int n = N_GLOBAL;
    vector<ll> pref(n + 1, 0);
    for (int i = 1; i <= n; ++i) pref[i] = pref[i - 1] + W[i];

    // We'll binary search for λ such that cnt(λ) <= K_GLOBAL
    // (for maximization: segments decrease as λ increases).
    for (int it = 0; it < 60; ++it) { // sufficient for 64-bit
        ll mid = (lambda_low + lambda_high) / 2;
        AlienRes res = alien_dp(mid, pref);
        if (res.cnt > K_GLOBAL) {
            // too many segments => λ too small
            lambda_low = mid;
        } else {
            // too few segments => λ too big
            lambda_high = mid;
        }
    }

    // After we find λ*, recompute DP at λ* and restore original objective:
    ll lambda_star = lambda_high;
    AlienRes res = alien_dp(lambda_star, pref);

    // original_value = val + λ * cnt
    // If we insist on exactly K_GLOBAL segments, we usually do:
    //   ans = res.val + lambda_star * K_GLOBAL;
    // but cnt may differ by 1 depending on monotonicity, so adjust if needed.
    ll answer = res.val + lambda_star * (ll)K_GLOBAL;
    return answer;
}

// =============================================================
//                          MAIN (example)
// =============================================================
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // ---- EXAMPLES OF HOW TO USE (comment out in contests) ----

    // 1) SOS DP example:
    /*
    int K = 3;
    vector<ll> f(1 << K);
    // fill f[mask] ...
    SOS::sum_over_submasks(f, K);   // f[mask] = sum_{sub ⊆ mask} orig_f[sub]
    SOS::sum_over_supersets(f, K);  // f[mask] = sum_{sup ⊇ mask} orig_f[sup]
    */

    // 2) Slope Trick example:
    /*
    SlopeTrick st;
    vector<ll> a = {1, 5, 2};
    for (ll x : a) st.add_abs(x);
    cout << "Min of sum |x - a_i| is " << st.get_min() << "\n";
    */

    // 3) Alien Trick example:
    /*
    N_GLOBAL = 5;
    K_GLOBAL = 2;
    W.assign(N_GLOBAL + 1, 0);
    // index from 1: W[1..N_GLOBAL]
    W[1] = 1; W[2] = 2; W[3] = -3; W[4] = 4; W[5] = 5;
    ll ans = alien_solve(-1000, 1000); // search range depends on problem
    cout << "Alien trick answer = " << ans << "\n";
    */

    return 0;
}
