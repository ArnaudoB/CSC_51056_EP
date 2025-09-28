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

class Point {
private:
    std::array<double, 2> coords;

public:
    // constructor
    Point(double x_val, double y_val) : coords{ x_val, y_val } {}

    // squared Euclidean distance
    static double sqdist(const Point& point_a, const Point& point_b) {
        double diff_x = point_a.coords[0] - point_b.coords[0];
        double diff_y = point_a.coords[1] - point_b.coords[1];
        return diff_x * diff_x + diff_y * diff_y;
    }

    // string method
    std::string str() const {
        std::ostringstream os;
        os << "(" << coords[0] << ", " << coords[1] << ")";
        return os.str();
    }
};


class HillClimbing{

    private:
    std::vector<Point> cloud;
    std::vector<std::vector<size_t>> neighbors;
    std::vector<double> density;
    std::vector<size_t> parent;
    std::vector<size_t> label;

    public:

    std::vector<Point> cloud;
    std::vector<std::vector<size_t>> neighbors;
    std::vector<double> density;
    std::vector<size_t> parent;

    // function to read the data
    void read_data(const std::string& filename) {
        std::ifstream in(filename);
        if (!in) {
            throw std::runtime_error("Could not open file: " + filename);
        }

        double x, y;
        cloud.clear();
        while (in >> x >> y) {
            cloud.emplace_back(x, y);
        }

        if (cloud.empty()) {
            throw std::runtime_error("File was read, but no points found!");
        }
    };

    // naive function to compute the k-NN graph of the points
    void compute_neighbors(int k) {
        const size_t n = cloud.size();
        neighbors.assign(n, {});

        for (size_t i = 0; i < n; i++) {
            std::vector<std::pair<double, size_t>> dist_idx;
            dist_idx.reserve(n);

            for (size_t j = 0; j < n; j++) {
                double d = Point::sqdist(cloud[i], cloud[j]);
                dist_idx.emplace_back(d, j);
            }

            // sort by distance
            std::sort(dist_idx.begin(), dist_idx.end(),
                    [](auto& a, auto& b){ return a.first < b.first; });

            // take the k closest neighbors
            for (int t = 1; t <= k && t < (int)dist_idx.size(); t++) {
                neighbors[i].push_back(dist_idx[t].second);
            }
        }
    };

    // density function
    void compute_density(int k) {
        if (neighbors.empty()) {
            compute_neighbors(k);
        }

        size_t n = cloud.size();
        density.assign(n, 0.0);

        for (size_t i = 0; i < n; i++) {
            double acc = 0.0;
            for (int idx : neighbors[i]) {
                acc += Point::sqdist(cloud[idx], cloud[i]);
            }

            if (!neighbors[i].empty()) {
                double mean_sq = acc / neighbors[i].size();
                density[i] = (mean_sq > 0.0) ? (1.0 / std::sqrt(mean_sq))
                                            : std::numeric_limits<double>::infinity();
            } else {
                density[i] = 0.0;
            }
        }
    };

    // forest function
    void compute_forest(int k) {
        if (density.empty()) {
            compute_density(k);
        }

        size_t n = cloud.size();
        parent.assign(n, 0);

        for (size_t i = 0; i < n; i++) {
            double max_density = density[i];
            size_t id_max = i;
            for (size_t idx : neighbors[i]) {
                if (density[idx] > max_density) {
                    max_density = density[idx];
                    id_max = idx;
                }
            }
            parent[i] = id_max;
        }
    }

    // labels computing 
    int find_root(int u) {
        if (parent[u] == u) return u;
        return parent[u] = find_root(parent[u]); // path compression
    };
    void compute_labels() {
        size_t n = cloud.size();
        label.assign(n, 0);
        for (size_t i = 0; i < n; i++) {
            label[i] = find_root(i);
        }
    };
};

int main() {
    Point p1(1.0, 2.0);
    Point p2(4.0, 6.0);

    std::cout << "p1 = " << p1.str() << "\n";
    std::cout << "p2 = " << p2.str() << "\n";
    std::cout << "sqdist = " << Point::sqdist(p1, p2) << "\n";
}
