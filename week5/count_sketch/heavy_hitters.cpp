#include <iostream>
#include <math.h>
#include <utility>
#include <stdio.h>
#include <time.h>
#include <map>
#include <cassert>
#include <random>  // For random number generation
#include <algorithm> // For std::nth_element
using namespace std;

// Seed for the random number engine (static for single initialization)
static random_device rd;
static default_random_engine eng;

// Define a big prime number
static const int big_prime = 1000000007;

class CountSketch {
private:
    const int num_hash_functions; // t
    const int num_buckets; // b

    vector<vector<long long>> counts;
    vector<int> a_coeffs;
    vector<int> b_coeffs;
    vector<int> c_coeffs;

    vector<int> get_random_vector(const int size, const int min_val, const int max_val) {
        // Define the range for random integers
        uniform_int_distribution<int> int_dist(min_val, max_val);

        // Initialize vector
        vector<int> random_ints_vector(size);

        // Populate the vector with random integers
        for (int i = 0; i < size; ++i) {
            random_ints_vector[i] = int_dist(eng);
        }

        return random_ints_vector;
    }

    // Hash Function: h(x) = ((ax + b) (mod p)) (mod m)
    int hash_function(const int a, const int b, const int x, const int p, const int m) {
        // Use long long for intermediate_result to prevent overflow
        long long intermediate_result = static_cast<long long>(a) * x + b;

        // First modulo operation: result_mod_p must be in [0, p-1]
        // Use the (value % p + p) % p pattern for correct positive modulo
        int result_mod_p = (intermediate_result % p + p) % p;

        // Second modulo operation: result must be in [0, m-1]
        return result_mod_p % m;
    }

    // Sign Hash Function: s(x) = cx (mod p)
    int sign_hash_function(const int c, const int x, const int p) {
        // Use long long for intermediate_result to prevent overflow
        long long intermediate_result = static_cast<long long>(c) * x;

        // First modulo operation: result_mod_p must be in [0, p-1]
        // Use the (value % p + p) % p pattern for correct positive modulo
        int result_mod_p = (intermediate_result % p + p) % p;

        // Second modulo operation: result must be in {-1, 1}
        return (result_mod_p % 2 == 0) ? 1 : -1;
    }

public:
    CountSketch(const int t, const int b) : num_hash_functions(t), num_buckets(b) {
        counts.resize(num_hash_functions, vector<long long>(num_buckets, 0));

        // Define vectors a, b and c to get the hash functions
        a_coeffs = get_random_vector(num_hash_functions, 1, big_prime - 1);
        b_coeffs = get_random_vector(num_hash_functions, 0, big_prime - 1);
        c_coeffs = get_random_vector(num_hash_functions, 1, big_prime - 1);
    }

    void add_data(const int item_id, const int frequency) {
        for (int row = 0; row < num_hash_functions; ++row) {
            int col = hash_function(a_coeffs[row], b_coeffs[row], item_id, big_prime, num_buckets);
            counts[row][col] += frequency * sign_hash_function(c_coeffs[row], item_id, big_prime);
        }
    }

    double estimate(const int item_id) {
        vector<long long> estimates;
        for (int row = 0; row < num_hash_functions; ++row) {
            int col = hash_function(a_coeffs[row], b_coeffs[row], item_id, big_prime, num_buckets);
            estimates.push_back(counts[row][col] * sign_hash_function(c_coeffs[row], item_id, big_prime));
        }
        
        // Take the median for all the rows
        sort(estimates.begin(), estimates.end());
        int middle_idx = num_hash_functions / 2;

        if (num_hash_functions % 2 == 1) {
            return estimates[middle_idx];
        }
        
        return (static_cast<double>(estimates[middle_idx - 1]) + static_cast<double>(estimates[middle_idx])) / 2;
    }
};

int main() {
    eng.seed(0);
    int n, t;
    scanf("%d %d", &n, &t);

    // num_hash_functions (t): Increased from typical theoretical bounds to improve confidence.
    // A higher 't' reduces the probability of a "bad" estimate due to adversarial inputs or bad luck.
    int number_hash_functions = 80;

    // num_buckets (b): This is the number of columns. A higher 'b' reduces the error magnitude
    // by minimizing hash collisions, which are the primary source of noise in sketch estimates.
    // The value 10000 was found through tuning to meet accuracy requirements for given constraints.
    const int number_buckets = 10000;

    // Initialize Count Sketch object
    CountSketch count_sketch(number_hash_functions, number_buckets);

    int id, value;
    // Read the first half of the stream, it corresponds to positive values
    for (int i = 0; i < n; ++i) {
        scanf("%d %d", &id, &value);
        count_sketch.add_data(id, value);
    }

    // Read the second half of the stream, it corresponds to negative values
    for (int i = 0; i < n; ++i) {
        scanf("%d %d", &id, &value);
        // Here we assign the negative sign
        count_sketch.add_data(id, -value);
    }

    // Read queries and output the estimate for each of them
    int num_queries = 0;
    scanf("%d", &num_queries);
    for (int q = 0; q < num_queries; ++q) {
        int query = 0;
        scanf("%d", &query);

        int estimate = count_sketch.estimate(query);
        if (estimate >= t)
            cout << "1 ";
        else
            cout << "0 ";
    }

    return 0;
}
