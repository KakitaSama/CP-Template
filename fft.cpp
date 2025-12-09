#include <bits/stdc++.h>
using namespace std;

using cd = complex<double>;
const double PI = acos(-1.0);

// Iterative FFT: O(n log n), in-place, no recursion
void fft(vector<cd> & a, bool invert) {
    int n = (int)a.size();

    // ----- Bit-reversal permutation -----
    // Rearrange elements in bit-reversed order
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;

        if (i < j)
            swap(a[i], a[j]);
    }

    // ----- Cooley–Tukey iterative -----
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        cd wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            cd w(1);
            for (int j = 0; j < len / 2; j++) {
                cd u = a[i + j];
                cd v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    // For inverse: divide by n at the end (more cache-friendly than dividing every stage)
    if (invert) {
        for (int i = 0; i < n; i++)
            a[i] /= n;
    }
}

// Convolution of two integer vectors (can be large n)
vector<long long> multiply(const vector<long long> &a, const vector<long long> &b) {
    vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n1 = (int)a.size(), n2 = (int)b.size();
    int need = n1 + n2 - 1;
    int n = 1;
    while (n < need) n <<= 1;
    fa.resize(n);
    fb.resize(n);

    fft(fa, false);
    fft(fb, false);
    for (int i = 0; i < n; i++)
        fa[i] *= fb[i];
    fft(fa, true);

    vector<long long> res(need);
    for (int i = 0; i < need; i++)
        res[i] = (long long)llround(fa[i].real());
    return res;
}
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m;
    cin >> n >> m;
    vector<long long> a(n), b(m);
    for (int i = 0; i < n; ++i) cin >> a[i];
    for (int i = 0; i < m; ++i) cin >> b[i];

    vector<long long> c = multiply(a, b); // size = n + m - 1

    for (int i = 0; i < (int)c.size(); ++i)
        cout << c[i] << " \n"[i + 1 == (int)c.size()];
}
