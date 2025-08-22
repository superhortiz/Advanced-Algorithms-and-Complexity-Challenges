#include <ios>
#include <iostream>
#include <vector>

using namespace std;

struct ConvertILPToSat {
    vector<vector<int>> A;
    vector<int> b;

    ConvertILPToSat(int n, int m) : A(n, vector<int>(m)), b(n)
    {}

    void printEquisatisfiableSatFormula() {
        vector<vector<int>> clauses;
        int n_ineq = A.size();
        int n_vars = A[0].size();
        vector<int> non_zero_coeff, curr_clause;

        // Iterate over inequalities
        for (int i = 0; i < n_ineq; ++i) {
            non_zero_coeff.clear();

            // Get the non zero coefficients. It is guaranteed that there will be at most
            // 3 different variables with non-zero coefficients in each inequality
            for (int j = 0; j < n_vars; ++j) {
                if (A[i][j] != 0) {
                    non_zero_coeff.push_back(j);
                }
            }

            // Get the number of non zero coefficients in the current inequality
            int p = non_zero_coeff.size();

            // This evaluates the case when all coefficients in the inequation i are 0
            if (p == 0) {
                // This implies that Ai * xi = 0 <= bi is unsatisfiable
                if (b[i] < 0) {
                    // This will include an empty clause, making the problem unsatisfiable
                    clauses.push_back({});
                    break;
                } else {
                    // There is nothing else to evaluate, continue to the next inequality
                    continue;
                }
            }

            // In this case we have just one coefficient. We evaluate when the associated variable is {0} {1}
            if (p == 1) {
                int a_index = non_zero_coeff[0];

                for (int x_a = 0; x_a <= 1; ++x_a) {
                    if (A[i][a_index] * x_a > b[i]) {
                        if (x_a == 0) clauses.push_back({a_index + 1});
                        else clauses.push_back({-(a_index + 1)});
                    }
                }

            // In this case we have two coefficient. We evaluate when the associated variables are {0 0} {0 1} {1 0} {1 1}
            } else if (p == 2) {
                int a_index = non_zero_coeff[0];
                int b_index = non_zero_coeff[1];

                for (int x_a = 0; x_a <= 1; ++x_a) {
                    for (int x_b = 0; x_b <= 1; ++x_b) {
                        if (A[i][a_index] * x_a + A[i][b_index] * x_b > b[i]) {
                            curr_clause.clear();
                            if (x_a == 0) curr_clause.push_back(a_index + 1); else curr_clause.push_back(-(a_index + 1));
                            if (x_b == 0) curr_clause.push_back(b_index + 1); else curr_clause.push_back(-(b_index + 1));
                            clauses.push_back(curr_clause);
                        }
                    }
                }

            // In this case we have three coefficient. We evaluate when the associated variables are
            // {0 0 0} {0 0 1} {0 1 0} {0 1 1} {1 0 0} {1 0 1} {1 1 0} {1 1 1}
            } else {
                int a_index = non_zero_coeff[0];
                int b_index = non_zero_coeff[1];
                int c_index = non_zero_coeff[2];

                for (int x_a = 0; x_a <= 1; ++x_a) {
                    for (int x_b = 0; x_b <= 1; ++x_b) {
                        for (int x_c = 0; x_c <= 1; ++x_c) {
                            if (A[i][a_index] * x_a + A[i][b_index] * x_b + A[i][c_index] * x_c  > b[i]) {
                                curr_clause.clear();
                                if (x_a == 0) curr_clause.push_back(a_index + 1); else curr_clause.push_back(-(a_index + 1));
                                if (x_b == 0) curr_clause.push_back(b_index + 1); else curr_clause.push_back(-(b_index + 1));
                                if (x_c == 0) curr_clause.push_back(c_index + 1); else curr_clause.push_back(-(c_index + 1));
                                clauses.push_back(curr_clause);
                            }
                        }
                    }
                }
            }
        }

        // This handles the case when the original ILP problem is trivially satisfiable.
        // Any combination of 0s and 1s for your variables would be a valid solution.
        // We add one clause to have a valid output that is satisfiable.
        if (clauses.size() == 0) clauses.push_back({1});

        // Output the clauses for SAT solver
        int n_clauses = clauses.size();
        std::cout << n_clauses << ' ' << n_vars << endl;
        for (auto& clause : clauses) {
            for (auto& var : clause) {
                std::cout << var << ' ';
            }
            std::cout << " 0" << endl;
        }
    }
};

int main() {
    ios::sync_with_stdio(false);

    int n, m;
    cin >> n >> m;

    ConvertILPToSat converter(n, m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cin >> converter.A[i][j];
        }
    }
    for (int i = 0; i < n; i++) {
        cin >> converter.b[i];
    }

    converter.printEquisatisfiableSatFormula();

    return 0;
}
