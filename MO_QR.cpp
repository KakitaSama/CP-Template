#include <bits/stdc++.h>
using namespace std;

struct Query {
	int l, r, idx;
};

int main() {
	int n;
	cin >> n;
	vector<int> v(n);
	for (int i = 0; i < n; i++) { cin >> v[i]; }

	int q;
	cin >> q;
	vector<Query> queries;
	for (int i = 0; i < q; i++) {
		int x, y;
		cin >> x >> y;
		queries.push_back({--x, --y, i});
	}

	int block_size = (int)sqrt(n);
	auto mo_cmp = [&](Query a, Query b) {
		int block_a = a.l / block_size;
		int block_b = b.l / block_size;
		if (block_a == block_b) { return a.r < b.r; }
		return block_a < block_b;
	};
	sort(queries.begin(), queries.end(), mo_cmp);

	int different_values = 0;
	vector<int> values(VALMAX);
	auto remove = [&](int idx) {
		values[v[idx]]--;
		if (values[v[idx]] == 0) { different_values--; }
	};
	auto add = [&](int idx) {
		values[v[idx]]++;
		if (values[v[idx]] == 1) { different_values++; }
	};

	int mo_left = -1;
	int mo_right = -1;
	vector<int> ans(q);
	for (int i = 0; i < q; i++) {
		int left = queries[i].l;
		int right = queries[i].r;

		while (mo_left < left) { remove(mo_left++); }
		while (mo_left > left) { add(--mo_left); }
		while (mo_right < right) { add(++mo_right); }
		while (mo_right > right) { remove(mo_right--); }

		ans[queries[i].idx] = different_values;
	}

	for (int i = 0; i < q; i++) { cout << ans[i] << '\n'; }
}
