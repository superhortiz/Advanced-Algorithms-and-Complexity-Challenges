#include <bits/stdc++.h>
using namespace std;

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

vector<int> strongly_connected_components(vector<set<int>> &adj, int &num_scc) {
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

struct Clause {
    int firstVar;
    int secondVar;
};

struct TwoSatisfiability {
    int numVars;
    vector<Clause> clauses;

    TwoSatisfiability(int n, int m) :
        numVars(n),
        clauses(m)
    {  }

    bool isSatisfiable_exponential(vector<int>& result) {
        for (int mask = 0; mask < (1 << numVars); ++mask) {
            // This iterates over all 2^numVars possible subsets.
            // Each 'mask' value uniquely represents a subset of items (0 to numVars-1).
            for (int i = 0; i < numVars; ++i) {
                // This inner loop iterates over each individual item/variable (from 0 to numVars-1).
                // For the current 'mask' (representing a specific subset), we check the status of each item 'i'.

                result[i] = (mask >> i) & 1;
                // This extracts the value of the i-th bit from 'mask'.
                // If the i-th bit is 1, it means item 'i' is IN the current subset represented by 'mask'.
                // If the i-th bit is 0, it means item 'i' is NOT IN the current subset.
                // So, 'result[i]' will be 1 if item 'i' is in the subset, and 0 otherwise.
            }

            bool formulaIsSatisfied = true;

            // Verify if every clause is satisfy by the current subset
            for (const Clause& clause : clauses) {
                bool clauseIsSatisfied = false;

                // Verify if the current clause is satisfied
                if (result[abs(clause.firstVar) - 1] == (clause.firstVar > 0)) {
                    clauseIsSatisfied = true;
                }

                if (result[abs(clause.secondVar) - 1] == (clause.secondVar > 0)) {
                    clauseIsSatisfied = true;
                }

                // We have found a clause that is not satisfied, we can break and continue with the next subset
                if (!clauseIsSatisfied) {
                    formulaIsSatisfied = false;
                    break;
                }
            }

            // We found one solution, we can return
            if (formulaIsSatisfied) {
                return true;
            }
        }

        // We iterated over 2^numVars possible solutions and we didnt find an assignment to satisfy the formula
        return false;
    }

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

    bool isSatisfiable(vector<int>& result) {
        // Build the Implication Graph based on the given 2-CNF formula
        // Even vertices correspond to a variable, odd vertices to their negation
        vector<set<int>> adj_list(2 * numVars);

        for (const Clause& clause : clauses) {
            adj_list[var_to_impl_vertex(-clause.firstVar)].emplace(var_to_impl_vertex(clause.secondVar));
            adj_list[var_to_impl_vertex(-clause.secondVar)].emplace(var_to_impl_vertex(clause.firstVar));
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

int main() {
    ios::sync_with_stdio(false);

    int n, m; // n variables, m clauses
    cin >> n >> m;
    TwoSatisfiability twoSat(n, m);
    for (int i = 0; i < m; ++i) {
        cin >> twoSat.clauses[i].firstVar >> twoSat.clauses[i].secondVar;
    }

    vector<int> result(n, -1);
    if (twoSat.isSatisfiable(result)) {
        cout << "SATISFIABLE" << endl;
        for (int i = 1; i <= n; ++i) {
            if (result[i-1]) {
                cout << i;
            } else {
                cout << -i;
            }
            if (i < n) {
                cout << " ";
            } else {
                cout << endl;
            }
        }
    } else {
        cout << "UNSATISFIABLE" << endl;
    }

    return 0;
}
