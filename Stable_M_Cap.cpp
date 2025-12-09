#include <bits/stdc++.h>
using namespace std;

// Stable matching with capacities: customers -> restaurants (SWERC 2020 L)

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N, M;
    if (!(cin >> N >> M)) return 0;

    // 0-based indices internally
    vector<long long> cap(M);
    for (int j = 0; j < M; ++j) cin >> cap[j];

    // consume end of line before using getline
    string line;
    getline(cin, line);

    // customers' preference lists (0..N-1 customers, 0..M-1 restaurants)
    vector<vector<int>> app_pref(N);
    for (int i = 0; i < N; ++i) {
        getline(cin, line);
        stringstream ss(line);
        int r;
        while (ss >> r) {
            // input restaurants are 1..M
            app_pref[i].push_back(r - 1);
        }
        // Problem guarantees each customer makes at least one reservation
    }

    // restaurants' ranking of customers
    // inst_rank[r][c] = rank of customer c for restaurant r (0 = best)
    // use sparse structure since total edges <= 1e6
    vector<unordered_map<int,int>> inst_rank(M);
    for (int r = 0; r < M; ++r) {
        getline(cin, line);
        stringstream ss(line);
        int x;
        if (!(ss >> x)) {
            // empty line => no reservations
            continue;
        }
        if (x == 0) {
            // restaurant r has no customers
            continue;
        } else {
            vector<int> customers;
            customers.push_back(x);
            while (ss >> x) customers.push_back(x);

            for (int pos = 0; pos < (int)customers.size(); ++pos) {
                int c = customers[pos] - 1; // customers are 1..N in input
                inst_rank[r][c] = pos;
            }
        }
    }

    // Gale–Shapley with capacities, customers propose

    vector<int> match_app(N, -1);          // customer -> restaurant
    vector<int> current_size(M, 0);        // how many customers currently assigned

    // For each restaurant: max-heap by rank (worst on top).
    using PII = pair<int,int>;             // (rank, customer)
    vector<priority_queue<PII>> pq(M);

    // For each customer: next restaurant index in app_pref to propose to
    vector<int> next_choice(N, 0);

    queue<int> free_customers;
    for (int c = 0; c < N; ++c) free_customers.push(c);

    while (!free_customers.empty()) {
        int c = free_customers.front();
        free_customers.pop();

        // If c has exhausted all options, stays unmatched
        if (next_choice[c] >= (int)app_pref[c].size()) continue;

        int r = app_pref[c][next_choice[c]++]; // restaurant index

        auto it = inst_rank[r].find(c);
        if (it == inst_rank[r].end()) {
            // Should not happen in valid input, but be safe:
            // restaurant didn't rank this customer -> treat as worst and reject.
            if (next_choice[c] < (int)app_pref[c].size())
                free_customers.push(c);
            continue;
        }
        int rank_c = it->second;

        if (current_size[r] < cap[r]) {
            // Restaurant r has free capacity: accept c
            match_app[c] = r;
            current_size[r]++;
            pq[r].push({rank_c, c});
        } else if (cap[r] > 0) {
            // Restaurant r is full, check if it prefers c over its worst current customer
            auto worst = pq[r].top();
            int worst_rank = worst.first;
            int worst_c    = worst.second;

            if (rank_c < worst_rank) {
                // r prefers c over worst_c: replace
                pq[r].pop();
                match_app[worst_c] = -1;
                free_customers.push(worst_c);

                match_app[c] = r;
                pq[r].push({rank_c, c});
                // current_size[r] unchanged
            } else {
                // r rejects c; c remains free if he has more options
                if (next_choice[c] < (int)app_pref[c].size())
                    free_customers.push(c);
            }
        } else {
            // cap[r] == 0 (degenerate, but allowed by statement? it says 1 ≤ c_i ≤ N, so not)
            // Still, handle gracefully:
            if (next_choice[c] < (int)app_pref[c].size())
                free_customers.push(c);
        }
    }

    // Collect matched customers (1-based ids), sorted
    vector<int> assigned;
    assigned.reserve(N);
    for (int c = 0; c < N; ++c) {
        if (match_app[c] != -1) assigned.push_back(c + 1);
    }
    sort(assigned.begin(), assigned.end());

    for (int id : assigned)
        cout << id << '\n';

    return 0;
}
