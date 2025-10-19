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
static inline std::vector<std::size_t>
symdiff(const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
    std::vector<std::size_t> r;
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


std::vector<std::array<std::size_t, 2>> gaussian_elim_sparse(const std::vector<std::array<std::size_t, 2>>& entries)
{
    if (entries.empty()) return {};

    // infer size
    std::size_t n = 0;
    for (const auto& e : entries)
        n = std::max(n, std::max(e[0], e[1]) + 1);

    // build columns: cols[j] = row indices with a 1 in column j
    std::vector<std::vector<std::size_t>> cols(n);
    cols.reserve(n);
    for (const auto& e : entries) cols[e[1]].push_back(e[0]);

    // sort each column so we can take lowest = back() and do linear-time xor
    for (auto& c : cols) std::sort(c.begin(), c.end());

    const std::size_t NONE = std::numeric_limits<std::size_t>::max();
    // which column has its pivot at row r (or NONE)
    std::vector<std::size_t> pivot_of_row(n, NONE);

    // left→right reduction: enforce unique lowest 1 per nonzero column
    for (std::size_t j = 0; j < n; ++j) {
        auto low = [&]() -> std::size_t {
            return cols[j].empty() ? NONE : cols[j].back(); // lowest = largest row
        };

        std::size_t r = low();
        while (r != NONE) {
            const std::size_t k = pivot_of_row[r];
            if (k == NONE) break;                // unique pivot row reached
            cols[j] = symdiff(cols[j], cols[k]);
            r = low();                           // recompute lowest 1
        }
        if (r != NONE) pivot_of_row[r] = j;      // register pivot for this row
    }

    // rebuild as sparse
    std::vector<std::array<std::size_t, 2>> out;
    std::size_t nnz = 0;
    for (const auto& c : cols) nnz += c.size();
    out.reserve(nnz);

    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i : cols[j]) out.push_back({i, j});

    return out;
}


// int main() {
//     auto F = read_filtration("./src/filtration.txt");
//     auto B = boundary_matrix_dense(F);
//     std::cout << "Boundary matrix size: " << B.size() << " x "
//               << (B.empty() ? 0 : B[0].size()) << "\n";

//     for (const auto& row : B) {
//         for (uint8_t v : row) std::cout << int(v);
//         std::cout << "\n";
//     }
//     gaussian_elim_dense(B);
//     std::cout << "Boundary matrix reduced size: " << B.size() << " x "
//               << (B.empty() ? 0 : B[0].size()) << "\n";

//     for (const auto& row : B) {
//         for (uint8_t v : row) std::cout << int(v);
//         std::cout << "\n";
//     }

//     auto B_sparse = boundary_matrix_sparse(F);

//     auto B_sparse_reduced = gaussian_elim_sparse(B_sparse);
//     std::cout << "Nonzero entries: " << B_sparse_reduced.size() << "\n";
//     for (auto [i, j] : B_sparse_reduced) {
//         std::cout << "(" << i << ", " << j << ")\n";
//     };

//     return 0;
// }