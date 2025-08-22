#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
using namespace std;

const int n_colors = 3;

int get_num_color(const char color) {
    switch (color) {
        case 'R': return 1; break;
        case 'G': return 2; break;
        case 'B': return 3; break;
        default: return -1; break;
    }
}

char get_color(const int num_color) {
    switch (num_color) {
        case 1: return 'R'; break;
        case 2: return 'G'; break;
        case 3: return 'B'; break;
        default: return '0'; break;
    }
}

// Function to get the variable number (1-based) given a vertice number (1-based)
// in the original graph and color
int vertex_to_var(const int vertex, const int color) {
    return n_colors * (vertex - 1) + color;
}

void make_clauses(const int n_vertices, const vector<pair<int, int>> &edges, const string colors, vector<vector<int>> &clauses) {
    // Compute total number 2-SAT variables
    int n_vars = n_colors * n_vertices;

    // Include clauses for adjacent vertices not to have the same color
    for (const auto& edge : edges) {
        int u = edge.first;
        int v = edge.second;
        for (int color = 1; color <= n_colors; ++color) {
            clauses.push_back({-vertex_to_var(u, color), -vertex_to_var(v, color)});
        }
    }

    // Include clauses not to have the same initial color
    for (int vertex = 1; vertex <= n_vertices; ++vertex) {
        clauses.push_back({-vertex_to_var(vertex, get_num_color(colors[vertex - 1]))});
    }

    // Include clauses for vertices to have at least one color
    vector<int> clause;
    for (int vertex = 1; vertex <= n_vertices; ++vertex) {
        clause.clear();

        for (int color = 1; color <= n_colors; ++color) {
            if (get_num_color(colors[vertex - 1]) != color) {
                clause.push_back(vertex_to_var(vertex, color));
            }
        }
        clauses.push_back(clause);

        // Include clauses for vertices to have just one color
        for (auto& literal : clause) {
            literal *= -1;
        }
        clauses.push_back(clause);
    }
}

void dfs_toposort(const vector<set<int>> &adj, vector<int> &used, vector<int> &order, int x) {
    used[x] = 1;
    for (const auto &adjacent : adj[x])
        if (used[adjacent] == 0)
            dfs_toposort(adj, used, order, adjacent);
    order.push_back(x);
}     

vector<int> toposort(const vector<set<int>> &adj) {
    vector<int> used(adj.size(), 0);
    vector<int> order;
    
    for (size_t vertex = 0; vertex < adj.size(); ++vertex) {
        if (used[vertex] == 0) {
            dfs_toposort(adj, used, order, vertex);
        }
    }

    reverse(order.begin(), order.end());
    return order;
}

vector<set<int>> reverse_graph(const vector<set<int>> &adj) {
    vector<set<int>> reversed_adj(adj.size(), set<int>());
    for (int v = 0; v < adj.size(); ++v) {
        for (const int w : adj[v]) {
            reversed_adj[w].emplace(v);
        }
    }
    return reversed_adj;
}

void dfs(const vector<set<int>> &adj, vector<int> &scc, int v, int num_scc) {
    scc[v] = num_scc;
    for (const auto &adjacent : adj[v])
        if (scc[adjacent] == -1)
            dfs(adj, scc, adjacent, num_scc);
}

vector<int> strongly_connected_components(const vector<set<int>> &adj, int &num_scc) {   
    // Reverse the graph
    vector<set<int>> reversed_graph = reverse_graph(adj);

    // Get a topological sort of the reversed graph,
    // which is equivalent to a reverse topological order of the original graph.
    vector<int> reverse_topo_order_original_graph = toposort(reversed_graph);

    // Use a vector to track visited nodes and store the corresponding scc
    vector<int> scc(adj.size(), -1);

    // Iterate over the original graph using the reverse postorder from the reversed graph
    for (const auto& vertex : reverse_topo_order_original_graph) {
        if (scc[vertex] == -1) {
            dfs(adj, scc, vertex, num_scc);
            // Each completed DFS traversal corresponds to a new strongly connected component
            num_scc++;
        }
    }
    return scc;
}

struct TwoSatisfiability {
    int numVars;
    TwoSatisfiability(int n): numVars(n) { }

    // Function to get the vertice number (0-based) in an Implication Graph given a variable number (1-based).
    int var_to_impl_vertex(const int var_num) {
        if (var_num < 0) {
            return 2 * (abs(var_num) - 1) + 1;
        }
        return 2 * (var_num - 1);
    }

    // Function to get the variable number (1-based) given a vertice number (0-based) in an Implication Graph.
    int impl_vertex_to_var(const int impl_vertex) {
        return impl_vertex / 2 + 1;
    }

    bool isSatisfiable(vector<int>& result, const vector<vector<int>> &clauses) {
        // Build the Implication Graph based on the given 2-CNF formula
        // Even vertices correspond to a variable, odd vertices to their negation.
        // 2 * numVars is the number of vertices in the Implication Graph
        vector<set<int>> adj_list(2 * numVars);

        // For each 1-clause (x1), introduce an edge -x1 → x1
        // For each 2-clause (x1 V x2), introduce two directed edges -x1 → x2 and -x2 → x1
        for (const auto& clause : clauses) {
            if (clause.size() == 1) {
                adj_list[var_to_impl_vertex(-clause[0])].emplace(var_to_impl_vertex(clause[0]));
            } else {
                adj_list[var_to_impl_vertex(-clause[0])].emplace(var_to_impl_vertex(clause[1]));
                adj_list[var_to_impl_vertex(-clause[1])].emplace(var_to_impl_vertex(clause[0]));
            }
        }

        // Get the strongly connected components of the Implication Graph
        int num_scc = 0;
        vector<int> scc = strongly_connected_components(adj_list, num_scc);

        // Check if the problem is unsatisfiable
        for (int var = 0; var < numVars; ++var) {
            // Check if there is any variable in the same SCC than its negation, then the problem is unsatisfiable
            if (scc[2 * var] == scc[2 * var + 1]) {
                return false;
            }
        }

        // Get the reverse topological order
        vector<vector<int>> reverse_topo_order(num_scc + 1);
        for (int vertex = 0; vertex < scc.size(); ++vertex) {
            int curr_scc = scc[vertex];
            reverse_topo_order[curr_scc].push_back(vertex);
        }

        // Iterate in reverse of the topological order
        for (const auto& curr_scc : reverse_topo_order) {
            for (const int vertex : curr_scc) {
                int var_num = impl_vertex_to_var(vertex);

                // If we already have a result for the corresponding variable, continue to the next vertex
                if (result[var_num - 1] == -1)
                    // For literals that dont have assignment, set them to 1 (transitivity) and set their negation to 0
                    result[var_num - 1] = (vertex % 2 == 0);
            }
        }

        return true;
    }
};

/*
  Arguments:
    * `n` - the number of vertices.
    * `edges` - list of edges, each edge is a pair (u, v), 1 <= u, v <= n.
    * `colors` - string consisting of `n` characters, each belonging to the set {'R', 'G', 'B'}.
  Return value:
    * If there exists a proper recoloring, return value is a string containing new colors, similar to the `colors` argument.
    * Otherwise, return value is an empty string.
*/
string assign_new_colors(const int n, const vector<pair<int, int>> &edges, const string &colors) {
    vector<vector<int>> clauses;
    make_clauses(n, edges, colors, clauses);
  
    // Compute total number of variables for 2-SAT
    int n_vars = n_colors * n;
    vector<int> result(n_vars, -1);

    TwoSatisfiability twoSat(n_vars);
    bool is_satisfiable = twoSat.isSatisfiable(result, clauses);

    if (is_satisfiable) {
        string answer;
        for (int i = 0; i < n; ++i) {
            for (int color = 1; color <= n_colors; ++color) {
                if (result[n_colors * i + color - 1] == 1) {
                    answer.push_back(get_color(color));
                    break;
                }
            }
        }
        return answer;
    }
    return "";
}

int main() {
    // Number of vertices (n) and number of edges (m) of the graph
    int n, m;
    cin >> n >> m;
    string colors;
    cin >> colors;

    vector<pair<int, int>> edges;
    for (int i = 0; i < m; i++) {
        int u, v;
        cin >> u >> v;
        edges.push_back(make_pair(u, v));
    }

    string new_colors = assign_new_colors(n, edges, colors);   

    if (new_colors.empty()) {
        cout << "Impossible" << endl;
    } else {
        cout << new_colors << endl;
    }
}
