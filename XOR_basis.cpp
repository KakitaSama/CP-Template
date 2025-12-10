#include <bits/stdc++.h>
#define i23 long long

vector<i23> basis;

void add(i23 x) {
	for (int i = 0; i < (int)basis.size();
	     i++) {  // reduce x using the current basis vectors
		x = min(x, x ^ basis[i]);
	}

	if (x != 0) { basis.push_back(x); }
}
