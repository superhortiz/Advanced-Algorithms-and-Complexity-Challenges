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
    
class StockCharts {
public:
    void Solve() {
        vector<vector<int>> stock_data = ReadData();
        int result = MinCharts(stock_data);
        WriteResponse(result);
    }

private:
    vector<vector<int>> ReadData() {
        int num_stocks, num_points;
        cin >> num_stocks >> num_points;
        vector<vector<int>> stock_data(num_stocks, vector<int>(num_points));
        for (int i = 0; i < num_stocks; ++i)
            for (int j = 0; j < num_points; ++j) {
                cin >> stock_data[i][j];
            }
        return stock_data;
    }

    void WriteResponse(int result) {
        cout << result << "\n";
    }

    FlowGraph build_graph(const vector<vector<int>>& stock_data) {
        // Build a graph based in the stock_data
        int num_left = stock_data.size();
        int num_right = stock_data.size();
        FlowGraph graph(num_left + num_right + 2);

        // Connect source with the left nodes (source node = 0)
        for (int i = 1; i <= num_left; ++i) {
            graph.add_edge(0, i, 1);
        }

        // Connect right nodes with the target (target node = num_left + num_right + 1)
        for (int j = num_left + 1; j <= num_left + num_right; ++j) {
            graph.add_edge(j, num_left + num_right + 1, 1);
        }

        // Connect left nodes with right nodes
        for (int i = 0; i < num_left; ++i) {
            for (int j = 0; j < num_right; ++j) {
                // Add edges just when the stock i is strictly less than the stock j
                if (i != j && compare(stock_data[i], stock_data[j])) {
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

    int MinCharts(const vector<vector<int>>& stock_data) {
        // Build the graph based on stock_data
        FlowGraph graph = build_graph(stock_data);
        
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

        // The minimum number of charts is equal to the number of stock minus the maximum flow
        return stock_data.size() - flow;
    }

    // Returns true if stock1's price is strictly less than
    // stock2's price at every corresponding time point.
    bool compare(const vector<int>& stock1, const vector<int>& stock2) {
        for (int i = 0; i < stock1.size(); ++i)
            if (stock1[i] >= stock2[i])
                return false;
        return true;
    }
};

int main() {
    std::ios_base::sync_with_stdio(false);
    StockCharts stock_charts;
    stock_charts.Solve();
    return 0;
}
