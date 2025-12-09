// Example integer function to minimize
long long f_int(int x) {
    // define your function for integer x, e.g.:
    // return 1LL * (x - 7) * (x - 7) + 5;
    return 0;
}

// Ternary search on integers to find MINIMUM of f_int on [L, R]
int ternary_search_int(int L, int R) {
    while (R - L > 3) { // keep shrinking until interval is small
        int m1 = L + (R - L) / 3;
        int m2 = R - (R - L) / 3;

        long long f1 = f_int(m1);
        long long f2 = f_int(m2);

        if (f1 < f2) {
            // Minimum in [L, m2]
            R = m2;
        } else {
            // Minimum in [m1, R]
            L = m1;
        }
    }

    // Now brute-force in [L, R]
    int bestX = L;
    long long bestVal = f_int(L);
    for (int x = L + 1; x <= R; ++x) {
        long long val = f_int(x);
        if (val < bestVal) {
            bestVal = val;
            bestX = x;
        }
    }
    return bestX; // argmin
}

// Example function to minimize
double f(double x) {
    // define your function here, e.g.:
    // return (x - 3.5) * (x - 3.5) + 10; // minimum near x = 3.5
    // In a real problem, this might use global data.
    return 0.0;
}

// Ternary search on doubles to find approximate MINIMUM of f
double ternary_search_double(double L, double R) {
    // Number of iterations controls precision.
    // 100–200 iterations is usually enough for 1e-9 precision.
    for (int it = 0; it < 200; ++it) {
        double m1 = L + (R - L) / 3.0;
        double m2 = R - (R - L) / 3.0;

        if (f(m1) < f(m2)) {
            // Minimum in [L, m2]
            R = m2;
        } else {
            // Minimum in [m1, R]
            L = m1;
        }
    }
    return (L + R) / 2.0; // approximate argmin
}
