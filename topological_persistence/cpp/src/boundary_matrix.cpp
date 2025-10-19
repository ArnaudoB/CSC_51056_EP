#include "read_filtration.hpp" 
#include <vector>
#include <array>
#include <set>
#include <algorithm>             // std::includes
#include <cstdint>               // uint8_t
#include <iostream>

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

// int main() {
//     auto F = read_filtration("./src/filtration.txt");
//     std::cout << "Read " << F.size() << " simplices.\n";

//     auto B = boundary_matrix_dense(F);
//     std::cout << "Boundary matrix size: " << B.size() << " x "
//               << (B.empty() ? 0 : B[0].size()) << "\n";

//     for (const auto& row : B) {
//         for (uint8_t v : row) std::cout << int(v);
//         std::cout << "\n";
//     }

//     auto entries = boundary_matrix_sparse(F);
//     std::cout << "Nonzero entries: " << entries.size() << "\n";
//     for (auto [i, j] : entries) {
//         std::cout << "(" << i << ", " << j << ")\n";
//     };
// }
