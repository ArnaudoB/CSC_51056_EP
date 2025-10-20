#include "read_filtration.hpp" 
#include <vector>
#include <array>
#include <map>
#include <set>
#include <algorithm>             // std::includes
#include <cstdint>               // uint8_t
#include <iostream>
#include <unordered_map>

static inline bool is_subset_eq(const std::set<int>& a, const std::set<int>& b) {
    return std::includes(b.begin(), b.end(), a.begin(), a.end());
}

std::vector<std::vector<uint8_t>> boundary_matrix_dense(const std::vector<simplex>& v_simplex) {
    const std::size_t n = v_simplex.size();
    std::vector<std::vector<uint8_t>> M(n, std::vector<uint8_t>(n, 0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            if (v_simplex[i].dim + 1 == v_simplex[j].dim &&
                is_subset_eq(v_simplex[i].vert, v_simplex[j].vert)) {
                M[i][j] = 1;
            }
        }
    }
    return M;
};

struct VecHash {
    size_t operator()(const std::vector<int>& v) const noexcept {
        size_t h = 1469598103934665603ull; // FNV-1a seed
        for (int x : v) {
            h ^= static_cast<size_t>(x) + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2);
        }
        return h;
    }
};

// TODO: make a better boundary_matrix function and make it handle sparse
// TODO: explore using bit masks for vertex encoding

std::vector<std::vector<uint8_t>> boundary_matrix_dense_fast(const std::vector<simplex>& F) {
    const size_t n = F.size();
    std::vector<std::vector<uint8_t>> M(n, std::vector<uint8_t>(n, 0));

    // Map each simplex's vertex set -> index
    std::unordered_map<std::vector<int>, int, VecHash> idx_of;
    idx_of.reserve(n * 2);
    for (int i = 0; i < (int)n; ++i) {
        // ensure sorted once if you can't guarantee it
        // auto v = S[i].vert; std::sort(v.begin(), v.end()); idx_of.emplace(std::move(v), i);
        std::vector<int> v(F[i].vert.begin(), F[i].vert.end());
        idx_of.emplace(std::move(v), i);
    }

    // For each column j (a k-simplex), generate all (k-1)-faces and set M[i][j]=1
    for (int j = 0; j < (int)n; ++j) {
        const auto& sj = F[j];
        if (sj.dim <= 0) continue;                              // 0-simplices have empty boundary
        std::vector<int> v(sj.vert.begin(), sj.vert.end());     // size = sj.dim + 1

        // Build each face = v without v[r]
        std::vector<int> face; face.reserve(v.size() - 1);
        for (size_t r = 0; r < v.size(); ++r) {
            face.clear();
            face.insert(face.end(), v.begin(), v.begin() + r); // vertices before v_r
            face.insert(face.end(), v.begin() + r + 1, v.end()); // vertices after v_r
            auto it = idx_of.find(face); // look-up the face in the dictionary
            if (it != idx_of.end()) {
                M[it->second][j] = 1;
            }
            // else: filtration may omit faces; ignore or handle as needed
        }
    }
    return M;
}

// sparse version

std::vector<std::array<std::size_t, 2>> boundary_matrix_sparse(const std::vector<simplex>& simplices) {
    std::vector<std::array<std::size_t, 2>> entries;
    const std::size_t n = simplices.size();

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            if (simplices[i].dim + 1 == simplices[j].dim &&
                is_subset_eq(simplices[i].vert, simplices[j].vert)) {
                entries.emplace_back(std::array<std::size_t, 2>{i, j});
            }
        }
    }
    return entries;

}

std::vector<std::vector<int>> boundary_matrix_sparse_fast(const std::vector<simplex>& simplices) {
    int n = simplices.size();

    // boundary_matrix initialization
    std::vector<std::vector<int>> boundary_matrix(n);

    // map each vertex-set to its index in the filtration
    std::map<std::set<int>, int> index;
    for (int i = 0; i < n; ++i)
        index[simplices[i].vert] = i;

    // fill: for each simplex, generate its (d-1)-faces by removing one vertex
    for (int j = 0; j < n; ++j) {
        const auto& vert = simplices[j].vert;
        const std::size_t d_plus_1 = vert.size();
        if (d_plus_1 <= 1) continue; // 0-simplex: empty column

        // build faces by skipping the r-th vertex in order
        for (std::size_t r = 0; r < d_plus_1; ++r) {
            std::set<int> face;
            std::size_t k = 0;
            for (int v : vert) {
                if (k != r) face.insert(v);
                ++k;
            }
            auto it = index.find(face);
            if (it != index.end()) {
                boundary_matrix[j].push_back(it->second);
            }
        }
        std::sort(boundary_matrix[j].begin(), boundary_matrix[j].end());
    }

    return boundary_matrix;
}


// int main() {
//     auto F = read_filtration("./src/filtrations/filtration_test.txt");
//     std::cout << "Read " << F.size() << " simplices.\n";

//     auto B = boundary_matrix_sparse_fast(F);
//     for (size_t j = 0; j < B.size(); ++j) {
//         std::cout << "Column " << j << " -> [ ";
//         for (int row : B[j])
//             std::cout << row << ' ';
//         std::cout << "]\n";
//     }

//     auto entries = boundary_matrix_sparse(F);
//     std::cout << "Nonzero entries: " << entries.size() << "\n";
//     for (auto [i, j] : entries) {
//         std::cout << "(" << i << ", " << j << ")\n";
//     };
// }
