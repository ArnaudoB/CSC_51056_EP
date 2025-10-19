#include <vector>
#include <cstdint>


void gaussian_elim_dense(std::vector<std::vector<uint8_t>>& M);

std::vector<std::array<std::size_t, 2>> gaussian_elim_sparse(const std::vector<std::array<std::size_t, 2>>& entries);