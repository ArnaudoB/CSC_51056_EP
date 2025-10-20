#include <vector>
#include <array>
#include <cstdint>
#include <cstddef>
#include <utility>
#include <algorithm> 
#include <limits>
#include "boundary_matrix.hpp"
#include <iostream>

static inline void col_xor(std::vector<std::vector<uint8_t>>& M,
                           std::size_t dst, std::size_t src) {
    const std::size_t m = M.size();
    for (std::size_t i = 0; i < m; ++i) M[i][dst] ^= M[i][src];
}

static inline std::size_t lowest_one(const std::vector<std::vector<uint8_t>>& M,
                                     std::size_t j) {
    const std::size_t m = M.size();
    for (std::size_t i = m; i-- > 0; ) if (M[i][j]) return i;
    return std::numeric_limits<std::size_t>::max();
}

void gaussian_elim_dense(std::vector<std::vector<uint8_t>>& M) {
    const std::size_t n = M.size();
    if (n == 0) return;

    const std::size_t NONE = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> pivot_of_row(n, NONE);

    for (std::size_t j = 0; j < n; ++j) {
        std::size_t low = lowest_one(M, j);

        while (low != NONE) {
            const std::size_t k = pivot_of_row[low];  // earlier pivot column at this row?
            if (k == NONE) break;                     // unique now
            col_xor(M, j, k);                         // eliminate that lowest 1
            low = lowest_one(M, j);                   // recompute
        }

        if (low != NONE) pivot_of_row[low] = j;       // register pivot for this row
    }
}


// sparse version


// symmetric difference of two sorted vectors
static inline std::vector<int>
symdiff(const std::vector<int>& a, const std::vector<int>& b) {
    std::vector<int> r;
    r.reserve(a.size() + b.size());
    std::size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] < b[j])       r.push_back(a[i++]);
        else if (b[j] < a[i])  r.push_back(b[j++]);
        else { ++i; ++j; }
    }
    while (i < a.size()) r.push_back(a[i++]);
    while (j < b.size()) r.push_back(b[j++]);
    return r;
}


std::vector<std::vector<int>>
gaussian_elim_sparse(std::vector<std::vector<int>> cols)
{
    const size_t n = cols.size();
    if (n == 0) return cols;

    const int NONE = -1;
    std::vector<int> pivot_of_row(n, NONE); // pivot_of_row[r] = column owning pivot at row r

    for (int j = 0; j < n; ++j) {
        auto low = [&]() -> int {
            return cols[j].empty() ? NONE : cols[j].back(); // "lowest" = largest row index
        };

        int r = low();
        while (r != NONE) {
            const int k = pivot_of_row[r];
            if (k == NONE) break;                 // free pivot row found
            cols[j] = symdiff(cols[j], cols[k]); // XOR columns
            r = low();                             // recompute lowest 1
        }
        if (r != NONE) pivot_of_row[r] = j; // record pivot
    }

    return cols; // reduced columns, same sparse format
}


// int main() {
//     auto F = read_filtration("./src/filtrations/filtration_test.txt");
//     auto B = boundary_matrix_sparse_fast(F);
//     std::cout << "Boundary matrix size: " << B.size() << " x "
//               << (B.empty() ? 0 : B[0].size()) << "\n";

//     for (const auto& row : B) {
//         for (uint8_t v : row) std::cout << int(v);
//         std::cout << "\n";
//     }
//     gaussian_elim_sparse(B);
//     std::cout << "Boundary matrix reduced size: " << B.size() << " x "
//               << (B.empty() ? 0 : B[0].size()) << "\n";

//     for (const auto& row : B) {
//         for (uint8_t v : row) std::cout << int(v);
//         std::cout << "\n";
//     }

//     auto B_sparse = boundary_matrix_dense(F);

//     gaussian_elim_dense(B_sparse);
//     for (const auto& row : B) {
//         for (uint8_t v : row) std::cout << int(v);
//         std::cout << "\n";
//     }

//     return 0;
// }