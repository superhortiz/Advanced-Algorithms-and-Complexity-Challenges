#include <ios>
#include <iostream>
#include <vector>

using namespace std;
const int n_colors = 3;

int var_number(const int vertex, const int color) {
    return n_colors * (vertex - 1) + color;
}

struct Edge {
    int from;
    int to;
};

struct ConvertGSMNetworkProblemToSat {
    int numVertices;
    vector<Edge> edges;

    ConvertGSMNetworkProblemToSat(int n, int m) :
        numVertices(n),
        edges(m)
    {}

    void printEquisatisfiableSatFormula() {
        // Compute variables
        int n_vars = n_colors * numVertices;
        int n_clauses = n_colors * edges.size() + (n_colors + 1) * numVertices;

        // Print the first line, number of clauses and number of variables
        cout << n_clauses << ' ' << n_vars << endl;

        // Include clauses for vertices to have at least one color
        for (int vertex = 1; vertex <= numVertices; ++vertex) {
            for (int color = 1; color <= n_colors; ++color) {
                cout << var_number(vertex, color) << ' ';
            }
            cout << 0 << endl;
        }

        // Include clauses for vertices to have just one color
        for (int vertex = 1; vertex <= numVertices; ++vertex) {
            for (int color1 = 1; color1 < n_colors; ++color1) {
                for (int color2 = color1 + 1; color2 <= n_colors; ++ color2) {
                    cout << -var_number(vertex, color1) << ' ' << -var_number(vertex, color2) << ' ' << 0 << endl;
                }
            }
        }

        // Include clauses for adjacent vertices not to have the same color
        for (auto& edge : edges) {
            for (int color = 1; color <= n_colors; ++color) {
                cout << -var_number(edge.from, color) << ' ' << -var_number(edge.to, color) << ' ' << 0 << endl;
            }
        }
    }
};

int main() {
    ios::sync_with_stdio(false);

    int n, m;
    cin >> n >> m;
    ConvertGSMNetworkProblemToSat converter(n, m);
    for (int i = 0; i < m; ++i) {
        cin >> converter.edges[i].from >> converter.edges[i].to;
    }

    converter.printEquisatisfiableSatFormula();

    return 0;
}
