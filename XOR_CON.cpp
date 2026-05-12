using ll = long long;

// ================= XOR FWHT =================
void fwht_xor(std::vector<ll>& a, bool inverse) {
    int n = (int)a.size();
    for (int len = 1; len < n; len <<= 1) {
        for (int i = 0; i < n; i += (len << 1)) {
            for (int j = 0; j < len; j++) {
                ll u = a[i + j];
                ll v = a[i + j + len];
                a[i + j] = u + v;
                a[i + j + len] = u - v;
            }
        }
    }
    if (inverse) {
        for (ll &x : a) x /= n;
    }
}

// ================= OR FWHT =================
void fwht_or(std::vector<ll>& a, bool inverse) {
    int n = (int)a.size();
    for (int len = 1; len < n; len <<= 1) {
        for (int i = 0; i < n; i += (len << 1)) {
            for (int j = 0; j < len; j++) {
                if (!inverse)
                    a[i + j + len] += a[i + j];
                else
                    a[i + j + len] -= a[i + j];
            }
        }
    }
}

// ================= AND FWHT =================
void fwht_and(std::vector<ll>& a, bool inverse) {
    int n = (int)a.size();
    for (int len = 1; len < n; len <<= 1) {
        for (int i = 0; i < n; i += (len << 1)) {
            for (int j = 0; j < len; j++) {
                if (!inverse)
                    a[i + j] += a[i + j + len];
                else
                    a[i + j] -= a[i + j + len];
            }
        }
    }
}

// ================= Generic Convolution Helpers =================
template<typename Transform>
std::vector<ll> convolution(std::vector<ll> a,
                            std::vector<ll> b,
                            Transform transform) {
    int n = 1;
    while (n < (int)std::max(a.size(), b.size())) n <<= 1;
    a.resize(n);
    b.resize(n);

    transform(a, false);
    transform(b, false);

    for (int i = 0; i < n; i++) a[i] *= b[i];

    transform(a, true);
    return a;
}

// Convenience wrappers:
std::vector<ll> xor_convolution(std::vector<ll> a, std::vector<ll> b) {
    return convolution(a, b, fwht_xor);
}

std::vector<ll> or_convolution(std::vector<ll> a, std::vector<ll> b) {
    return convolution(a, b, fwht_or);
}

std::vector<ll> and_convolution(std::vector<ll> a, std::vector<ll> b) {
    return convolution(a, b, fwht_and);
}
int main(){
    ios::sync_with_stdio(false);
    cin.tie(0);
    #ifndef ONLINE_JUDGE 
    setIO();
    #endif
    int n;
    cin>>n;

    vector<ll> a((1<<20));
    int x;
    vector<int> v(n+1);
    long long S = 0;
    bool sp = false;
    a[0] = 1;
    for(int i = 1; i <= n ; i++)
    {
        cin>>x;
        S ^= x;
        a[S]++;
        if(a[S] > 1 || S == 0)
            sp = true;
    }
    vector<ll> d = xor_convolution(a, a);
    vector<int> ans;
    if(sp)
    ans.push_back(0);
    for(int i = 1 ; i < (1<<20) ; i++)
    {
        if(d[i])
            ans.push_back(i);
    }
    cout<<ans.size()<<endl;
    for(int x : ans)
        cout<<x<<" ";
    cout<<endl;
}
