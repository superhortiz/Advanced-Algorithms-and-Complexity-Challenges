#include <bits/stdc++.h>
using namespace std;

struct Edge {
    int from;
    int to;
};

struct ConvertHampathToSat {
    int numVertices;
    vector<Edge> edges;

    ConvertHampathToSat(int n, int m) :
        numVertices(n),
        edges(m)
    {  }

    int num_var(int vertex, int pos_path) {
        return (vertex - 1) * numVertices + pos_path;
    }

    void printEquisatisfiableSatFormula() {
        // Preprocess edges to get the adjacency list
        vector<unordered_set<int>> adjacency_list(numVertices + 1);
        int n_edges = 0;

        for (auto& edge : edges) {
            if (adjacency_list[edge.from].count(edge.to) == 0) {
                adjacency_list[edge.from].emplace(edge.to);
                adjacency_list[edge.to].emplace(edge.from);
                n_edges++;
            }
        }

        // Compute the total number of clauses and variables
        int n_clauses1 = numVertices;
        int n_clauses2 = numVertices * numVertices * (numVertices - 1) / 2;
        int n_clauses3 = numVertices;
        int n_clauses4 = numVertices * numVertices * (numVertices - 1) / 2;
        int n_clauses5 = (numVertices * (numVertices - 1) - 2 * n_edges) * (numVertices - 1);
        int n_clauses = n_clauses1 + n_clauses2 + n_clauses3 + n_clauses4 + n_clauses5;
        int n_vars = numVertices * numVertices;
        cout << n_clauses << ' ' << n_vars << endl;

        // Clauses:
        // Each vertex belongs to a path
        for (int vertex = 1; vertex <= numVertices; ++vertex) {
            for (int pos_path = 1; pos_path <= numVertices; ++pos_path) {
                cout << num_var(vertex, pos_path) << ' ';
            }
            cout << 0 << endl;
        }

        // Each vertex appears just once in a path
        for (int vertex = 1; vertex <= numVertices; ++vertex) {
            for (int pos_path1 = 1; pos_path1 < numVertices; ++pos_path1) {
                for (int pos_path2 = pos_path1 + 1; pos_path2 <= numVertices; ++ pos_path2) {
                    cout << -num_var(vertex, pos_path1) << ' ' << -num_var(vertex, pos_path2) << " 0" << endl;
                }
            }
        }

        // Each position in a path is occupied by some vertex
        for (int pos_path = 1; pos_path <= numVertices; ++pos_path) {
            for (int vertex = 1; vertex <= numVertices; ++vertex) {
                cout << num_var(vertex, pos_path) << ' ';
            }
            cout << 0 << endl;
        }

        // No two vertices occupy the same position of a path
        for (int pos_path = 1; pos_path <= numVertices; ++pos_path) {
            for (int vertex1 = 1; vertex1 < numVertices; ++vertex1) {
                for (int vertex2 = vertex1 + 1; vertex2 <= numVertices; ++vertex2) {
                    cout << -num_var(vertex1, pos_path) << ' ' << -num_var(vertex2, pos_path) << " 0" << endl;
                }
            }
        }

        // Two successive vertices on a path must be connected by an edge
        // (In other words, if there is no edge {𝑖, 𝑗} in 𝐸, then for any 𝑘, it cannot be the case that both 𝑥𝑖𝑘 and 𝑥𝑗(𝑘+1) are True.)
        for (int vertex1 = 1; vertex1 < numVertices; ++vertex1) {
            for (int vertex2 = vertex1 + 1; vertex2 <= numVertices; ++vertex2) {
                if (adjacency_list[vertex1].count(vertex2) == 0) {
                    for (int k = 1; k < numVertices; ++k) {
                        cout << -num_var(vertex1, k) << ' ' << -num_var(vertex2, k + 1) << " 0" << endl;
                        cout << -num_var(vertex1, k + 1) << ' ' << -num_var(vertex2, k) << " 0" << endl;
                    }
                }
            }
        }
    }
};

int main() {
    ios::sync_with_stdio(false);

    int n, m;
    cin >> n >> m;
    ConvertHampathToSat converter(n, m);
    for (int i = 0; i < m; ++i) {
        cin >> converter.edges[i].from >> converter.edges[i].to;
    }

    converter.printEquisatisfiableSatFormula();

    return 0;
}
