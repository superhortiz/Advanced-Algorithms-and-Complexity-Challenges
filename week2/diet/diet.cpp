#include <algorithm>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cmath>
#include <limits>
using namespace std;

typedef std::vector<double> Column;
typedef std::vector<double> Row;
typedef std::vector<Row> Matrix;

const double precision = 0.001;
const double max_total_amount = 1000000000;
const double min_val = -std::numeric_limits<double>::max();

struct Position {
    Position(int column, int row): column(column), row(row) {}
    int column;
    int row;
};

Position SelectPivotElement(const Matrix &a, std::vector <bool> &used_rows, std::vector <bool> &used_columns) {
    // Select the pivot as the leftmost non-zero element in a non-pivot row
    Position pivot_element(0, 0);
    int size = a.size();

    for (int col = 0; col < size; ++col) {
        for (int row = 0; row < size; ++ row) {
            if (a[row][col] != 0 && !used_rows[row] && !used_columns[col]) {
                pivot_element.row = row;
                pivot_element.column = col;
                return pivot_element;
            }
        }
    }
    return pivot_element;
}

void SwapLines(Matrix &a, Column &b, std::vector<bool> &used_rows, Position &pivot_element) {
    std::swap(a[pivot_element.column], a[pivot_element.row]);
    std::swap(b[pivot_element.column], b[pivot_element.row]);

    bool temp = used_rows[pivot_element.column];
    used_rows[pivot_element.column] = used_rows[pivot_element.row];
    used_rows[pivot_element.row] = temp;

    pivot_element.row = pivot_element.column;
}

void ProcessPivotElement(Matrix &a, Column &b, const Position &pivot_element) {
    int size = a.size();
    int p = pivot_element.row;
    int q = pivot_element.column;

    // Scale coefficients
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i != p && j != q) {
                a[i][j] -= a[p][j] * a[i][q] / a[p][q];
            }
        }
        if (i != p)
            b[i] -= b[p] * a[i][q] / a[p][q];
    }

    // Make all other elements in the Pivot Column Equal to 0
    for (int i = 0; i < size; ++i) {
        if (i != p)
            a[i][q] = 0; 
    }

    // Scale row p by dividing by the coefficient a[p][q]
    for (int j = 0; j < size; ++j) {
        if (j != q)
            a[p][j] /= a[p][q];
    }
    b[p] /= a[p][q];

    // Set the coefficient a[p][q] to 1
    a[p][q] = 1.0;
}

void MarkPivotElementUsed(const Position &pivot_element, std::vector <bool> &used_rows, std::vector <bool> &used_columns) {
    used_rows[pivot_element.row] = true;
    used_columns[pivot_element.column] = true;
}

Column SolveEquation(Matrix &a, Column &b) {
    int size = a.size();

    std::vector <bool> used_columns(size, false);
    std::vector <bool> used_rows(size, false);
    for (int step = 0; step < size; ++step) {
        Position pivot_element = SelectPivotElement(a, used_rows, used_columns);
        SwapLines(a, b, used_rows, pivot_element);
        ProcessPivotElement(a, b, pivot_element);
        MarkPivotElementUsed(pivot_element, used_rows, used_columns);
    }

    for (int i = 0; i < size; ++i) {
        if (a[i][i] == 0) return std::vector<double>(size, std::nan("1"));
    }

    return b;
}

void get_subsets(vector<vector<int>>& subsets, vector<int> curr_subset, const int curr, const int size_set, const int size_subset) {
    // Function to recursively get all the possible subsets of m (size_subset)
    // inequalities from the total set of n (size_set) inequalities
    if (curr_subset.size() == size_subset) {
        subsets.push_back(curr_subset);
        return;
    }

    if (curr == size_set) return;

    get_subsets(subsets, curr_subset, curr + 1, size_set, size_subset);
    curr_subset.push_back(curr);
    get_subsets(subsets, curr_subset, curr + 1, size_set, size_subset);
}

bool check_inequalities(const int n, const int m, const Matrix& A, const vector<double>& b, const Column& x, const vector<bool>& in_equation) {
    // Check that all the variables are equal or greater than zero
    for (int i = 0; i < m; ++i)
        if (x[i] + precision < 0)
            return false;

    // Check all the inequalities that are not in the current equation (Ax <= b)
    for (int i = 0; i < n; ++i) {
        if (!in_equation[i]) {
            double lhs = 0;
            for (int j = 0; j < m; ++j) {
                lhs += A[i][j] * x[j];
            }
            if (lhs > b[i] + precision) {
                return false;
            }
        }
    }
    return true;
}

double get_objective(const int m, const vector<double>& c, const Column& x) {
    // Get the current value of the objective function c*v
    double result = 0;
    for (int i = 0; i < m; ++i) {
        result += c[i] * x[i];
    }
    return result;
}

double get_total_amount(const Column& x) {
    // Get the current value of sum(x)
    int size = x.size();
    double total_amount = 0;

    for (int i = 0; i < size; ++i) {
        total_amount += x[i];
    }
    return total_amount;
}

pair<int, vector<double>> solve_diet_problem(const int n, const int m, const Matrix& A, const vector<double>& b, const vector<double>& c) {
    // Initialize data structures
    Matrix a_equation(m, std::vector<double>(m));
    Column b_equation(m);
    double best_objective = min_val;
    Column best_solution(m);
    vector<bool> in_equation(m + n, false);

    // Get the subsets from all the inequalities
    vector<vector<int>> subsets;
    get_subsets(subsets, {}, 0, m + n + 1, m);

    // Iterate over the subsets
    for (const auto& subset : subsets) {
        // Restart in_equation vector, to know which inequalities we need to check
        fill(in_equation.begin(), in_equation.end(), false);

        // Populate vectors a_equation and b_equation according to the current subset,
        // to solve the linear equation: a_equation * x = b_equation
        for (int i = 0; i < m; ++i) {
            int index = subset[i];
            in_equation[index] = true;

            if (index < n) {
                a_equation[i] = A[index];
                b_equation[i] = b[index];
            } else if (index < m + n) {
                Row row(m, 0);
                row[index - n] = 1;
                a_equation[i] = row;
                b_equation[i] = 0;
            } else {
                Row row(m, 1);
                a_equation[i] = row;
                b_equation[i] = max_total_amount;
            }
        }

        // Get the solution of the linear equation
        Column x = SolveEquation(a_equation, b_equation);

        // Get the current objective function
        double curr_objective = get_objective(m, c, x);

        // Check if we can improve the best solution so far
        if (check_inequalities(n, m, A, b, x, in_equation) && 
            curr_objective > best_objective || (curr_objective == best_objective && get_total_amount(best_solution) == max_total_amount)) {
            best_solution = x;
            best_objective = curr_objective;
        }
    }

    // Check if there is solution
    if (best_objective == min_val) {
        return {-1, vector<double>(m, 0)};
    }

    // Check if the solution is bounded
    if (get_total_amount(best_solution) + precision >= max_total_amount) {
        return {1, vector<double>(m, 0)};
    }

    // Return the optimal solution
    return {0, best_solution};
}

int main(){
    int n, m;
    cin >> n >> m;
    Matrix A(n, vector<double>(m));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cin >> A[i][j];
        }
    }
    vector<double> b(n);
    for (int i = 0; i < n; i++) {
        cin >> b[i];
    }
    vector<double> c(m);
    for (int i = 0; i < m; i++) {
        cin >> c[i];
    }

    pair<int, vector<double>> ans = solve_diet_problem(n, m, A, b, c);

    switch (ans.first) {
        case -1: 
            printf("No solution\n");
            break;
        case 0: 
            printf("Bounded solution\n");
            for (int i = 0; i < m; i++) {
                printf("%.18f%c", ans.second[i], " \n"[i + 1 == m]);
            }
            break;
        case 1:
            printf("Infinity\n");
            break;
    }
    return 0;
}
