#include <iostream>
#include <vector>
#include <limits>

const int INF = std::numeric_limits<int>::max();

struct Vertex {
    int weight;
    std::vector<int> children;
};
typedef std::vector<Vertex> Graph;
typedef std::vector<int> Sum;

Graph ReadTree() {
    int vertices_count;
    std::cin >> vertices_count;

    Graph tree(vertices_count);

    for (int i = 0; i < vertices_count; ++i)
        std::cin >> tree[i].weight;

    for (int i = 1; i < vertices_count; ++i) {
        int from, to, weight;
        std::cin >> from >> to;
        tree[from - 1].children.push_back(to - 1);
        tree[to - 1].children.push_back(from - 1);
    }

    return tree;
}

void dfs(const Graph &tree, std::vector<int> &dp, int vertex, int parent) {
    // If the value for vertex is already processed, then return
    if (dp[vertex] != INF) return;

    int m1 = tree[vertex].weight;

    // If vertex doesnt have children then dp[vertex] is just its weight
    if (tree[vertex].children.size() == 0) {
        dp[vertex] = m1;
    } else {
        // Compute m1, which includes vertex and its grandchildren
        for (int child : tree[vertex].children) {
            if (child != parent) {
                for (int grandchild : tree[child].children) {
                    if (grandchild != vertex) {
                        dfs(tree, dp, grandchild, child);
                        m1 += dp[grandchild];
                    }
                }
            }
        }

        // Compute m0, which includes the children of vertex
        int m0 = 0;

        for (int child : tree[vertex].children) {
            if (child != parent) {
                dfs(tree, dp, child, vertex);
                m0 += dp[child];
            }
        }

        // Get the maximum and store it
        dp[vertex] = std::max(m0, m1);
    }
}

int MaxWeightIndependentTreeSubset(const Graph &tree) {
    size_t size = tree.size();
    if (size == 0)
        return 0;

    // Define the dynamic programming table to store results
    std::vector<int> dp(size, INF);

    // We start exploring from vertex 0
    dfs(tree, dp, 0, -1);

    // Return the maximum value
    return dp[0];
}

int main() {
    Graph tree = ReadTree();
    int weight = MaxWeightIndependentTreeSubset(tree);
    std::cout << weight << std::endl;
    return 0;
}
