#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <queue>
#include <limits>

using std::vector;
using std::cin;
using std::cout;

/* This class implements a bit unusual scheme for storing edges of the graph,
 * in order to retrieve the backward edge for a given edge quickly. */
class FlowGraph {
public:
    struct Edge {
        int from, to, capacity, flow;
    };

public:
    /* List of all - forward and backward - edges */
    vector<Edge> edges;

    /* These adjacency lists store only indices of edges in the edges list */
    vector<vector<size_t>> graph;

public:
    explicit FlowGraph(size_t n): graph(n) {}

    void add_edge(int from, int to, int capacity) {
        /* Note that we first append a forward edge and then a backward edge,
            * so all forward edges are stored at even indices (starting from 0),
            * whereas backward edges are stored at odd indices in the list edges */
        Edge forward_edge = {from, to, capacity, 0};
        // For backward edges, keep the capacity equals to zero, then
        // residual_capacity = capacity - flow in both cases.
        // Flow in backward edges is negative.
        Edge backward_edge = {to, from, 0, 0};

        // Add edges to the vectors graph and edges
        graph[from].push_back(edges.size());
        edges.push_back(forward_edge);
        graph[to].push_back(edges.size());
        edges.push_back(backward_edge);
    }

    size_t size() const {
        return graph.size();
    }

    const vector<size_t>& get_ids(int from) const {
        return graph[from];
    }

    const Edge& get_edge(size_t id) const {
        return edges[id];
    }

    void add_flow(size_t id, int flow) {
        /* To get a backward edge for a true forward edge (i.e id is even), we should get id + 1
            * due to the described above scheme. On the other hand, when we have to get a "backward"
            * edge for a backward edge (i.e. get a forward edge for backward - id is odd), id - 1
            * should be taken.
            *
            * It turns out that id ^ 1 works for both cases. Think this through! */
        edges[id].flow += flow;
        edges[id ^ 1].flow -= flow;
    }
};

class MaxMatching {
public:
    void Solve() {
        vector<vector<bool>> adj_matrix = ReadData();
        vector<int> matching = FindMatching(adj_matrix);
        WriteResponse(matching);
    }

private:
    vector<vector<bool>> ReadData() {
        int num_left, num_right;
        cin >> num_left >> num_right;
        vector<vector<bool>> adj_matrix(num_left, vector<bool>(num_right));
        for (int i = 0; i < num_left; ++i)
            for (int j = 0; j < num_right; ++j) {
                int bit;
                cin >> bit;
                adj_matrix[i][j] = (bit == 1);
            }
        return adj_matrix;
    }

    void WriteResponse(const vector<int>& matching) {
        for (int i = 0; i < matching.size(); ++i) {
            if (i > 0)
                cout << " ";
            if (matching[i] == -1)
                cout << "-1";
            else
                cout << (matching[i] + 1);
        }
        cout << "\n";
    }

    FlowGraph build_graph(const vector<vector<bool>>& adj_matrix) {
        // Build a graph based in the adj_matrix
        int num_left = adj_matrix.size();
        int num_right = adj_matrix[0].size();
        FlowGraph graph(num_left + num_right + 2);

        // Connect source with the left nodes
        for (int i = 1; i <= num_left; ++i) {
            graph.add_edge(0, i, 1);
        }

        // Connect right nodes with the target
        for (int j = num_left + 1; j <= num_left + num_right; ++j) {
            graph.add_edge(j, num_left + num_right + 1, 1);
        }

        // Connect left nodes with right nodes
        for (int i = 0; i < num_left; ++i) {
            for (int j = 0; j < num_right; ++j) {
                if (adj_matrix[i][j] == true) {
                    graph.add_edge(i + 1, j + num_left + 1, 1);
                }
            }
        }
        return graph;
    }

    void bfs(FlowGraph& graph, vector<int>& edge_to, int from) {
        // Perform BFS to find paths from source to target
        std::queue<int> queue;
        queue.push(from);
        edge_to[from] = 0;
    
        while (!queue.empty()) {
            int v = queue.front();
            queue.pop();
            const vector<size_t>& ids = graph.get_ids(v);
            for (const auto& id : ids) {
                const FlowGraph::Edge& edge = graph.get_edge(id);
                if (edge_to[edge.to] == -1 && edge.capacity - edge.flow > 0) {
                    queue.push(edge.to);
                    edge_to[edge.to] = id; // Get track of the edges
                }
            }
        }
    }

    vector<int> FindMatching(const vector<vector<bool>>& adj_matrix) {
        // Build the graph based on adj_matrix
        FlowGraph graph = build_graph(adj_matrix);
        
        // Initialize vector to keep track of the edges in the found path
        size_t n = graph.size();
        vector<int> edge_to(n);

        // Initialize flow
        int flow = 0;
        int from = 0;
        int to = n - 1;

        while (true) {
            // Reset vector for each iteration
            std::fill(edge_to.begin(), edge_to.end(), -1);
    
            // Use BFS to find the shortest path from source (from) to target (to)
            bfs(graph, edge_to, from);
    
            // If there is no path to target, break the while loop
            if (edge_to[to] == -1) break;
    
            // Find the minimum capacity in the found path
            int x = std::numeric_limits<int>::max();
            int curr = to;
            while (curr != from) {
                const FlowGraph::Edge& edge = graph.get_edge(edge_to[curr]);
                x = std::min(x, edge.capacity - edge.flow);
                curr = edge.from;
            }
    
            // Update the flow for edges in the path
            curr = to;
            while (curr != from) {
                const FlowGraph::Edge& edge = graph.get_edge(edge_to[curr]);
                graph.add_flow(edge_to[curr], x);
                curr = edge.from;
            }
    
            // Update total flow
            flow += x;
        }

        // Initialize vector to store the results
        int num_left = adj_matrix.size();
        vector<int> matching(num_left, -1);

        // Iterate over the left nodes
        for (int i = 1; i <= num_left; ++i) {
            const vector<size_t>& ids = graph.get_ids(i);

            // Iterate over adjacent nodes of each left node
            for (const auto& id : ids) {
                const FlowGraph::Edge edge = graph.get_edge(id);
                // Store when we find an edge with positive flow
                if (edge.flow == 1) {
                    matching[i - 1] = edge.to - num_left - 1;
                }
            }
        }
       
        return matching;
    }
};

int main() {
    std::ios_base::sync_with_stdio(false);
    MaxMatching max_matching;
    max_matching.Solve();
    return 0;
}
