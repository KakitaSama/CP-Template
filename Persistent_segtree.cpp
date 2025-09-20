#include <bits/stdc++.h>
using namespace std;

using ll = long long;

class PersistentSegtree {
  private:
	struct Node {
		ll sum = 0;
		int l = 0, r = 0;
	};

	const int n;
	vector<Node> tree;
	int timer = 1;

	Node join(int l, int r) { return Node{tree[l].sum + tree[r].sum, l, r}; }

	int build(int tl, int tr, const vector<int> &arr) {
		if (tl == tr) {
			tree[timer] = {arr[tl], 0, 0};
			return timer++;
		}

		int mid = (tl + tr) / 2;
		tree[timer] = join(build(tl, mid, arr), build(mid + 1, tr, arr));

		return timer++;
	}

	int set(int v, int pos, int val, int tl, int tr) {
		if (tl == tr) {
			tree[timer] = {val, 0, 0};
			return timer++;
		}

		int mid = (tl + tr) / 2;
		if (pos <= mid) {
			tree[timer] = join(set(tree[v].l, pos, val, tl, mid), tree[v].r);
		} else {
			tree[timer] = join(tree[v].l, set(tree[v].r, pos, val, mid + 1, tr));
		}

		return timer++;
	}

	ll range_sum(int v, int ql, int qr, int tl, int tr) {
		if (qr < tl || tr < ql) { return 0ll; }
		if (ql <= tl && tr <= qr) { return tree[v].sum; }

		int mid = (tl + tr) / 2;
		return range_sum(tree[v].l, ql, qr, tl, mid) +
		       range_sum(tree[v].r, ql, qr, mid + 1, tr);
	}

  public:
	PersistentSegtree(int n, int MX_NODES) : n(n), tree(MX_NODES) {}

	int build(const vector<int> &arr) { return build(0, n - 1, arr); }

	int set(int root, int pos, int val) { return set(root, pos, val, 0, n - 1); }

	ll range_sum(int root, int l, int r) { return range_sum(root, l, r, 0, n - 1); }

	int add_copy(int root) {
		tree[timer] = tree[root];
		return timer++;
	}
};

int main() {
	int n, q;
	cin >> n >> q;
	vector<int> a(n);
	for (int &i : a) { cin >> i; }

	const int MX_NODES = 2 * n + q * (2 + __lg(n));
	PersistentSegtree st(n, MX_NODES);

	vector<int> roots = {st.build(a)};

	for (int t = 0; t < q; t++) {
		int type, k;
		cin >> type >> k;
		k--;

		if (type == 1) {
			int pos, val;
			cin >> pos >> val;
			pos--;
			roots[k] = st.set(roots[k], pos, val);
		} else if (type == 2) {
			int a, b;
			cin >> a >> b;
			a--, b--;
			cout << st.range_sum(roots[k], a, b) << '\n';
		} else if (type == 3) {
			roots.push_back(st.add_copy(roots[k]));
		}
	}
}
