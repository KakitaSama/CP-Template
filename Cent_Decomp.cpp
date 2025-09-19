	const int INF = 1e9;

vector<vector<int>> adj;
vector<int> subtree_size;
// min_dist[v] := the minimal distance between v and a red node
vector<int> min_dist;
vector<bool> is_removed;
vector<vector<pair<int, int>>> ancestors;

int get_subtree_size(int node, int parent = -1) {
	subtree_size[node] = 1;
	for (int child : adj[node]) {
		if (child == parent || is_removed[child]) { continue; }
		subtree_size[node] += get_subtree_size(child, node);
	}
	return subtree_size[node];
}

int get_centroid(int node, int tree_size, int parent = -1) {
	for (int child : adj[node]) {
		if (child == parent || is_removed[child]) { continue; }
		if (subtree_size[child] * 2 > tree_size) {
			return get_centroid(child, tree_size, node);
		}
	}
	return node;
}

/**
 * Calculate the distance between current `node` and the `centroid` it belongs
 * to. The distances between a node and all its centroid ancestors are stored
 * in the vector `ancestors`.
 * @param cur_dist the distance between `node` and `centroid`
 */
void get_dists(int node, int centroid, int parent = -1, int cur_dist = 1) {
	for (int child : adj[node]) {
		if (child == parent || is_removed[child]) { continue; }
		cur_dist++;
		get_dists(child, centroid, node, cur_dist);
		cur_dist--;
	}
	ancestors[node].push_back({centroid, cur_dist});
}

void build_centroid_decomp(int node = 0) {
	int centroid = get_centroid(node, get_subtree_size(node));

	/*
	 * For all nodes in the subtree rooted at `centroid`, calculate their
	 * distances to the centroid
	 */
	for (int child : adj[centroid]) {
		if (is_removed[child]) { continue; }
		get_dists(child, centroid, centroid);
	}

	is_removed[centroid] = true;
	for (int child : adj[centroid]) {
		if (is_removed[child]) { continue; }
		// build the centroid decomposition for all child components
		build_centroid_decomp(child);
	}
}

/**
 * Paint `node` red by updating all of its ancestors' minimal distances
 * to a red node
 */
void paint(int node) {
	for (auto &[ancestor, dist] : ancestors[node]) {
		min_dist[ancestor] = min(min_dist[ancestor], dist);
	}
	min_dist[node] = 0;
}

/** Print the minimal distance between `node` to a red node */
void query(int node) {
	int ans = min_dist[node];
	for (auto &[ancestor, dist] : ancestors[node]) {
		if (!dist) { continue; }
		/*
		 * The distance between `node` and a red painted node is the sum of
		 * the distance from `node` to one of its ancestors (`dist`) and the
		 * distance from this ancestor to the nearest red node
		 * (`min_dist[ancestor]`).
		 */
		ans = min(ans, dist + min_dist[ancestor]);
	}

	cout << ans << "\n";
}
