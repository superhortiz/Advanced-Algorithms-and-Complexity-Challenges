#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_set>

using std::vector;
typedef vector<vector<int> > Matrix;

const int INF = 1e9;

Matrix read_data() {
    int vertex_count, edge_count;
    std::cin >> vertex_count >> edge_count;
    Matrix graph(vertex_count, vector<int>(vertex_count, INF));
    for (int i = 0; i < edge_count; ++i) {
        int from, to, weight;
        std::cin >> from >> to >> weight;
        --from, --to;
        graph[from][to] = graph[to][from] = weight;
    }
    return graph;
}

std::pair<int, vector<int> > optimal_path(const Matrix& graph) {
    // Define source node equal to zero
    int s = 0;

    // Get the number of nodes in the graph
    int n = graph.size();

    // Initialize data structures
    std::unordered_set<int> subset, new_subset;
    std::vector<std::vector<int>> dp(n, std::vector<int>(1 << n, INF));
    std::vector<std::vector<int>> prev(n, std::vector<int>(1 << n, -1));

    // Populate the dynamic programming table (dp) with the base cases
    for (int i = 0; i < n; ++i) {
        if (i == s) continue;  // Skip the loop over the source node
        subset.insert(1 << s | 1 << i);  // Get the set of states to process
        dp[i][1 << s | 1 << i] = graph[s][i];  // Populate the table
        prev[i][1 << s | 1 << i] = s;  // Table to reconstruct the optimal path
    }

    // Iterate until visiting all the nodes in the graph
    for (int subset_size = 2; subset_size < n; ++subset_size) {
        // Iterate over the current processed nodes
        for (const auto& mask : subset) {
            // Iterate over all the possible nodes to add to the path
            for (int i = 0; i < n; ++i) {
                // Get the new state
                int new_mask = mask | (1 << i);

                // Skip in case the node i is already in the path
                if (new_mask == mask) continue;

                // Insert new_mask to new_subset to be processed in the next step
                new_subset.insert(new_mask);

                // Iterate over all predecessor of new_mask to find the optimal path
                for (int j = 0; j < n; ++j) {
                    // Prevents from trying to use j as a predecessor when j is not a member of the set of nodes represented by the mask
                    if (!(mask & (1 << j))) continue;
                    if (dp[i][new_mask] > dp[j][mask] + graph[j][i]) {
                        dp[i][new_mask] = dp[j][mask] + graph[j][i];
                        prev[i][new_mask] = j;  // Store predecessor
                    }
                }
            }
        }
        // Update subsets
        subset = new_subset;
        new_subset.clear();
    }

    // Iterate over the last state to get the optimal cost
    int result = INF;
    int last_node = -1;
    int last_mask = (1 << n) - 1;
    for (int i = 0; i < n; ++i) {
        if (result > dp[i][last_mask] + graph[i][s]) {
            result = dp[i][last_mask] + graph[i][s];
            last_node = i;
        }
    }

    // Reconstruct the optimal path
    std::vector<int> path;
    int curr_node = last_node;
    int curr_mask = last_mask;

    while (curr_node != -1) {
        path.push_back(curr_node);
        int next_node = prev[curr_node][curr_mask];
        curr_mask ^= (1 << curr_node);
        curr_node = next_node;
    }
    
    std::reverse(path.begin(), path.end());
    // path.push_back(s);

    if (result == INF)
        result = -1;
    for (size_t i = 0; i < path.size(); ++i)
        ++path[i];
    return std::make_pair(result, path);
}

void print_answer(const std::pair<int, vector<int>>& answer) {
    std::cout << answer.first << "\n";
    if (answer.first == -1)
        return;
    const vector<int>& path = answer.second;
    for (size_t i = 0; i < path.size(); ++i)
        std::cout << path[i] << " ";
    std::cout << "\n";
}

int main() {
    print_answer(optimal_path(read_data()));
    return 0;
}
