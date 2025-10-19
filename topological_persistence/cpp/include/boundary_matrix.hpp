#pragma once
#include <set>
#include <vector>
#include <cstdint>
#include "read_filtration.hpp"

std::vector<std::array<std::size_t, 2>> boundary_matrix_sparse(const std::vector<simplex>& simplices);

std::vector<std::vector<uint8_t>> boundary_matrix_dense(const std::vector<simplex>& v_simplex);
