#include <cmath>
#include <iostream>
#include <vector>

const double EPS = 1e-6;
const int PRECISION = 20;

typedef std::vector<double> Column;
typedef std::vector<double> Row;
typedef std::vector<Row> Matrix;

struct Equation {
    Equation(const Matrix &a, const Column &b): a(a), b(b) {}
    Matrix a;
    Column b;
};

struct Position {
    Position(int column, int row): column(column), row(row) {}
    int column;
    int row;
};

Equation ReadEquation() {
    int size;
    std::cin >> size;
    Matrix a(size, std::vector <double> (size, 0.0));
    Column b(size, 0.0);
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < size; ++column)
            std::cin >> a[row][column];
        std::cin >> b[row];
    }
    return Equation(a, b);
}

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

Column SolveEquation(Equation equation) {
    Matrix &a = equation.a;
    Column &b = equation.b;
    int size = a.size();

    std::vector <bool> used_columns(size, false);
    std::vector <bool> used_rows(size, false);
    for (int step = 0; step < size; ++step) {
        Position pivot_element = SelectPivotElement(a, used_rows, used_columns);
        SwapLines(a, b, used_rows, pivot_element);
        ProcessPivotElement(a, b, pivot_element);
        MarkPivotElementUsed(pivot_element, used_rows, used_columns);
    }

    return b;
}

void PrintColumn(const Column &column) {
    int size = column.size();
    std::cout.precision(PRECISION);
    for (int row = 0; row < size; ++row)
        std::cout << column[row] << std::endl;
}

int main() {
    Equation equation = ReadEquation();
    Column solution = SolveEquation(equation);
    PrintColumn(solution);
    return 0;
}
