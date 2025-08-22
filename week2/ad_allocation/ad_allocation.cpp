#include <algorithm>
#include <iostream>
#include <vector>
#include <cstdio>
using namespace std;

typedef vector<vector<double>> matrix;

void build_tableau(const int n, const int m, const matrix& A, const vector<double>& b, const vector<double>& c, matrix& tableau, vector<int>& basis) {
    // Put values of matrix A into simplex tableau
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            tableau[i][j] = A[i][j];
        }
    }

    // Put identity matrix for slack variables into simplex tableau
    for (int j = n; j < m + n; ++j) {
        tableau[j - n][j] = 1;
    }

    // Put c into simplex tableau
    for (int j = 0; j < n; ++j) {
        tableau[m][j] = c[j];
    }

    // Put b into simplex tableau
    for (int i = 0; i < m; ++i) {
        tableau[i][m + n] = b[i];
    }

    // Store initial basic variables
    for (int i = n; i < n + m; ++i) {
        basis[i - n] = i;
    }
}

void build_aux_tableau(const int n, const int m, const int n_artificial, const matrix& A, const vector<double>& b, matrix& tableau, vector<int>& basis) {
    // Put values of matrix A into simplex tableau
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            tableau[i][j] = A[i][j];
        }
    }

    // Put identity matrix for slack variables into simplex tableau
    for (int j = n; j < m + n; ++j) {
        tableau[j - n][j] = 1;
    }

    // Put auxiliar objective function into simplex tableau
    for (int j = m + n; j < m + n + n_artificial; ++j) {
        tableau[m][j] = -1;
    }

    // Put b into simplex tableau
    for (int i = 0; i < m; ++i) {
        tableau[i][m + n + n_artificial] = b[i];
    }

    // Put artificial variables into simple tableau. Store initial basic variables
    int count_artificial = 0;
    for (int i = 0; i < m; ++i) {
        if (b[i] >= 0) {
            basis[i] = i + n;
        } else {
            basis[i] = m + n + count_artificial;
            tableau[i][m + n + count_artificial++] = -1;
        }
    }
}

int bland_rule(const int n, const int m, const matrix& tableau) {
    // Applies Bland's rule. Returns the index of first column whose objective function coefficient is positive.
    // Bland's Rule guarantees that the Simplex Algorithm will not get into a cycle.

    for (int q = 0; q < m + n; ++q) {
        if (tableau[m][q] > 0) {
            return q;
        }
    }

    // If there are no positive coefficients, it means that we have the optimal
    return -1;
}

int min_ratio_rule(const int q, const int n, const int m, const matrix& tableau) {
    // Criterion to determine the leaving variable (the basic variable that will become non-basic) during each iteration.
    // Its primary purpose is to ensure that the new basic feasible solution remains feasible,
    // meaning all basic variables remain non-negative.

    // Start considering p = -1
    int p = -1;

    // Given the column q, iterate over all the rows
    for (int i = 0; i < m; i++) {
        // Consider only positive entries
        if (tableau[i][q] <= 0) {
            continue;
        // Use i as first value (to have something to compare later)
        } else if (p == -1) {
            p = i;
        // Compare the current best ratio with the ratio of the current row.
        // Keep the lowest value.
        } else if (tableau[i][m + n] / tableau[i][q] < tableau[p][m + n] / tableau[p][q])
            p = i;
    }

    // If it returns -1, it means that the problem is unbounded
    return p;
}

void pivot(const int p, const int q, const int n, const int m, matrix& tableau, vector<int>& basis) {
    // Function to perform pivot operation

    // Scale the entries. The coefficient A[p][j]/A[p][q] is the new coefficient in (p, j), but normalized.
    // We need to amplify it by A[i][q] and substract, to be able to make the coefficient A[i][q] equal to zero
    for (int i = 0; i <= m; i++) {
        for (int j = 0; j <= m + n; j++) {
            if (i != p && j != q) {
                tableau[i][j] -= tableau[p][j] * tableau[i][q] / tableau[p][q];
            }
        }
    }

    // Make all other elements in the pivot column equal to 0, for all rows of A
    for (int i = 0; i <= m; i++) {
        if (i != p) {
            tableau[i][q] = 0;
        }
    }

    // Scale row p by dividing by the coefficient [p][q], for all columns in A and I
    for (int j = 0; j <= m + n; j++) {
        if (j != q) {
            tableau[p][j] /= tableau[p][q];
        }
    }

    // After all these operations the coefficient [p][q] should be 1
    tableau[p][q] = 1;

    // Update basis variables
    basis[p] = q;
}

bool b_has_negatives(const vector<double>& b) {
    // Function to check if there is any inequality "greater than"
    for (int i = 0; i < b.size(); ++i) {
        if (b[i] < 0) return true;
    }
    return false;
}

void solve(const int n, const int m, matrix& tableau, vector<int>& basis, bool& is_bounded) {
    // Implements Simplex algorithm
    while (true) {
        // Get entering column using Bland's rule (first positive)
        int q = bland_rule(n, m, tableau);

        // Optimal solution found (no positive coefficients in objective row)
        if (q == -1) {
            return;
        }

        // Get leaving row using minimum ratio rule
        int p = min_ratio_rule(q, n, m, tableau);

        // The problem is unbounded
        if (p == -1) {
            is_bounded = false;
            return;
        }

        // Perform pivot operation on row p, column q
        pivot(p, q, n, m, tableau, basis);
    }
}

void two_phases(const int n, const int m, const matrix& A, const vector<double>& b, const vector<double>& c, matrix& tableau, vector<int>& basis, bool& is_feasible) {
    // Vector to keep track of the artificial variables
    vector<int> artificial_vars;

    // Iterate through b to find inequalities "greater than"
    for (int i = 0; i < m; ++i) {
        if (b[i] < 0) {
            // Keep track of the inequalities where the artificial variable is basic
            artificial_vars.push_back(i);
        }
    }

    // Initialize auxiliar data structure
    int n_artificial = artificial_vars.size();
    matrix new_tableau(m + 1, vector<double>(m + n + n_artificial + 1, 0));

    // Build auxiliar tableau
    build_aux_tableau(n, m, n_artificial, A, b, new_tableau, basis);

    // Remove the basic variables (artificial variables) from the objective function
    for (const auto& artificial_var : artificial_vars) {
        for (int j = 0; j <= m + n + n_artificial; ++j) {
            // Normalize row containing artificial variable and substract from objective function row
            new_tableau[artificial_var][j] *= -1;
            new_tableau[m][j] += new_tableau[artificial_var][j];
        }
    }

    // Solve the auxiliar problem using Simplex algorithm
    bool is_bounded = true;
    solve(n + n_artificial, m, new_tableau, basis, is_bounded);

    // Check if the problem is feasible
    if (abs(new_tableau[m][m + n + n_artificial]) >= 0.0001) {
        is_feasible = false;
        return;
    }

    // Copy results from auxiliar tableau to the original tableau
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m + n; ++j) {
            tableau[i][j] = new_tableau[i][j];
        }
        tableau[i][m + n] = new_tableau[i][m + n + n_artificial];
    }
   
    // Remove basic variables from objective function.
    for (int i = 0; i < m; ++i) {
        // Note: Using the condition basis[i] < m + n is necessary in case that an artificial
        // variable is still basic (then we get out of bounds of the columns in tableau).
        // It might happen when the two-phase simplex method is reaching an optimal objective
        // value of Z = 0 and some artificial variables might still be basic variables but have a value of zero.
        // This might happen because of Degeneracy.
        if (tableau[m][basis[i]] != 0 && basis[i] < m + n) {
            for (int j = 0; j <= m + n; ++j) {
                if (j != basis[i]) {
                    tableau[m][j] -= tableau[i][j] * tableau[m][basis[i]];
                }
            }
            tableau[m][basis[i]] = 0;
        }
    }

    return;
}

pair<int, vector<double>> allocate_ads(const int n, const int m, const matrix A, const vector<double> b, const vector<double> c) {
    // Initialize data structures
    matrix tableau(m + 1, vector<double>(m + n + 1, 0)); // Matrix to store simplex tableau
    vector<int> basis(m); // Vector to store current basic variables
    vector<double> x(n, 0); // Vector to store the solution
    bool is_feasible = true;

    // Build simplex tableau
    build_tableau(n, m, A, b, c, tableau, basis);

    // Check if there is any inequality "greater than"
    if (b_has_negatives(b)) {
        // Get a basic feasible solution using Two-Phase Simplex Method
        two_phases(n, m, A, b, c, tableau, basis, is_feasible);
        if (!is_feasible) {
            return {-1, x};
        }
    }

    // Solve the problem using simplex algorithm
    bool is_bounded = true;
    solve(n, m , tableau, basis, is_bounded);

    // Check if the solution is unbounded
    if (!is_bounded) {
        return {1, x};
    }

    // Retrieve the solution from final simplex tableau. If the variable is basic, its value is in the RHS.
    // If the variable is not basic, its value is zero.
    for (int i = 0; i < m; ++i) {
        if (basis[i] < n) {
            x[basis[i]] = tableau[i][m + n];
        }
    }
    return {0, x};
}

int main() {
    // n variables, m restrictions
    int n, m;
    cin >> m >> n;
    matrix A(m, vector<double>(n));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }
    vector<double> b(m);
    for (int i = 0; i < m; i++) {
        cin >> b[i];
    }
    vector<double> c(n);
    for (int i = 0; i < n; i++) {
        cin >> c[i];
    }

    pair<int, vector<double>> ans = allocate_ads(n, m, A, b, c);

    switch (ans.first) {
        case -1: 
            printf("No solution\n");
            break;
        case 0: 
            printf("Bounded solution\n");
            for (int i = 0; i < n; i++) {
                printf("%.18f%c", ans.second[i], " \n"[i + 1 == m]);
            }
            break;
        case 1:
            printf("Infinity\n");
            break;
    }
    return 0;
}
