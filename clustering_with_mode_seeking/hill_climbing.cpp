#include <array>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <numeric>
#include <utility>
#include <set>
#include <unordered_set>

class Point {
private:
    std::array<double, 2> coords;

public:
    Point(double x_val, double y_val) : coords{ x_val, y_val } {}

    static double sqdist(const Point& a, const Point& b) {
        double dx = a.coords[0] - b.coords[0];
        double dy = a.coords[1] - b.coords[1];
        return dx*dx + dy*dy;
    }

    std::string str() const {
        std::ostringstream os;
        os << "(" << coords[0] << ", " << coords[1] << ")";
        return os.str();
    }
};

class HillClimbing {
private:
    std::vector<Point> cloud;
    std::vector<std::vector<size_t>> neighbors;
    std::vector<double> density;
    std::vector<size_t> parent;
    std::vector<size_t> label;

public:
    // accessors
    const std::vector<Point>& points()   const { return cloud; }
    const std::vector<double>& densities() const { return density; }
    const std::vector<size_t>& labels()    const { return label; }

    void read_data(const std::string& filename) {
        std::ifstream in(filename);
        if (!in) throw std::runtime_error("Could not open file: " + filename);

        cloud.clear();
        double x, y;
        while (in >> x >> y) cloud.emplace_back(x, y);

        if (cloud.empty())
            throw std::runtime_error("File was read, but no points found!");
    }

    void compute_neighbors(size_t k) {
        const size_t n = cloud.size();
        if (n == 0) return;
        // at most n-1 neighbors
        if (k > 0 && k > n - 1) k = n - 1;

        neighbors.assign(n, {});

        for (size_t i = 0; i < n; i++) {
            std::vector<std::pair<double, size_t>> dist_idx;
            dist_idx.reserve(n);

            for (size_t j = 0; j < n; j++) {
                double d = Point::sqdist(cloud[i], cloud[j]);
                dist_idx.emplace_back(d, j);
            }

            std::sort(dist_idx.begin(), dist_idx.end(),
                      [](const auto& a, const auto& b){ return a.first < b.first; });

            // skip self at position 0
            for (size_t t = 1; t <= k && t < dist_idx.size(); t++)
                neighbors[i].push_back(dist_idx[t].second);
        }
    }

    void compute_density(size_t k) {

        const size_t n = cloud.size();
        density.assign(n, 0.0);

        for (size_t i = 0; i < n; i++) {
            double acc = 0.0;
            for (size_t idx=0; idx<k; idx++) {
                acc += Point::sqdist(cloud[neighbors[i][idx]], cloud[i]);
            }
            if (!neighbors[i].empty()) {
                double mean_sq = acc / (double)neighbors[i].size();
                density[i] = (mean_sq > 0.0)
                    ? (1.0 / std::sqrt(mean_sq))
                    : std::numeric_limits<double>::infinity();
            } else {
                density[i] = 0.0;
            }
        }
    }

    void compute_forest(size_t k) {

        const size_t n = cloud.size();
        parent.assign(n, 0);

        for (size_t i = 0; i < n; i++) {
            double best = density[i];
            size_t best_id = i;
            for (size_t idx=0; idx<k; idx++) {
                if (density[neighbors[i][idx]] > best) {
                    best = density[idx];
                    best_id = idx;
                }
            }
            parent[i] = best_id;
        }
    }

    size_t find_root(size_t u) {
        if (parent[u] == u) return u;
        return parent[u] = find_root(parent[u]); // path compression
    }

    void compute_labels() {
        const size_t n = cloud.size();
        label.assign(n, 0);
        for (size_t i = 0; i < n; i++) label[i] = find_root(i);
    }

    int compute_nb_labels() {
        std::unordered_set<size_t> unique_labels;
        for (size_t lab : label) {
            unique_labels.insert(lab);
        }
        return (int)unique_labels.size();
    }

    std::vector<size_t> sort_by_density() const {
        const size_t n = cloud.size();
        std::vector<size_t> indices(n);
        
        for (size_t i = 0; i < n; i++) {
            indices[i] = i;
        }
        
        std::sort(indices.begin(), indices.end(), 
                  [this](size_t a, size_t b) {
                      return density[a] > density[b];
                  });
        
        return indices;
    }

    struct PeakInfo {
    size_t peak_id;
    double birth;
    double death;
    double prominence;
    };


    std::set<double> compute_persistence_tau(int k, double tau) {
        const size_t n = cloud.size();
        if (n == 0) return {};

        std::vector<size_t> P(n);
        std::iota(P.begin(), P.end(), 0);
        std::sort(P.begin(), P.end(), [&](size_t i, size_t j){
            if (density[i] != density[j]) return density[i] > density[j];
            return i < j;
        });

        parent.resize(n);
        for (size_t i = 0; i < n; ++i) parent[i] = i;

        std::vector<char> marked(n, 0);

        std::set<double> pers;

        for (size_t t = 0; t < P.size(); ++t) {
            size_t pi = P[t];

            size_t p = find_root(pi);

            const auto& knn = neighbors[pi];
            const size_t K = std::min<size_t>(k, knn.size());
            for (size_t j = 0; j < K; ++j) {
                size_t qj = find_root(knn[j]);
                p = find_root(p);

                if (qj == p) continue;

                size_t m, M;
                if (density[p] < density[qj]) {
                    m = p;  M = qj;
                } else {
                    m = qj; M = p;
                }

                if (density[m] < density[pi] + tau) {
                    parent[m] = M;
                    p = find_root(pi);
                } else {
                    if (!marked[m]) {
                        pers.insert(density[m] - density[pi]);
                        marked[m] = 1;
                    }
                }
            }
        }

        return pers;
    }
};

int main() {
        HillClimbing hc;
        hc.read_data("./test.xy");

        size_t k_density = 10;
        size_t k_graph   = 5;
        double tau       = 0.35;

        hc.compute_neighbors(k_density);
        hc.compute_density(k_density);
        hc.compute_forest(k_graph);
        hc.compute_labels();

        std::cout << "Number of labels : "
                  << hc.compute_nb_labels() << std::endl;

        auto pers = hc.compute_persistence_tau(k_graph, tau);
        std::cout << "Persistences of non-merged peaks :\n";
        for (double p : pers) {
            if (std::isinf(p))
                std::cout << "+inf\n";
            else
                std::cout << p << "\n";
    }
}