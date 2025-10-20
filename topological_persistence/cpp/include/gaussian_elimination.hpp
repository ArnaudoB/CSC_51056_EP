#include <vector>
#include <cstdint>


void gaussian_elim_dense(std::vector<std::vector<uint8_t>>& M);

std::vector<std::vector<int>> gaussian_elim_sparse(std::vector<std::vector<int>> entries);