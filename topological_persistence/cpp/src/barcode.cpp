#include <vector>
#include <cstdint>
#include <array>
#include <fstream>
#include <algorithm>
#include <limits>
#include <string>
#include <gaussian_elimination.hpp>
#include <read_filtration.hpp>
#include <boundary_matrix.hpp>

static inline std::size_t lowest_one(const std::vector<std::vector<uint8_t>>& R,
                                     std::size_t j) {
    const std::size_t m = R.size();
    for (std::size_t i = m; i-- > 0; )
        if (R[i][j]) return i;
    return std::numeric_limits<std::size_t>::max();
}

void write_barcode_from_reduced_dense(const std::vector<std::vector<uint8_t>>& R,
                                      const std::vector<int>& dim,
                                      const std::vector<float>& val,
                                      const std::string& out_path)
{
    const std::size_t n = R.size();
    if (n == 0) { std::ofstream(out_path).close(); return; }

    const std::size_t NONE = std::numeric_limits<std::size_t>::max();
    std::vector<bool> row_used_as_pivot(n, false);

    std::ofstream ofs(out_path);
    ofs.setf(std::ios::fixed, std::ios::floatfield);
    ofs.precision(6);

    // finite intervals from paired columns
    for (std::size_t j = 0; j < n; ++j) {
        std::size_t i = lowest_one(R, j);
        if (i != NONE) {
            row_used_as_pivot[i] = true;
            int k = dim[j] - 1;
            ofs << k << ' ' << val[i] << ' ' << val[j] << '\n';
        }
    }
    // infinite intervals from unpaired zero columns that are not pivot-rows
    for (std::size_t j = 0; j < n; ++j) {
        if (lowest_one(R, j) == NONE && !row_used_as_pivot[j]) {
            int k = dim[j];
            ofs << k << ' ' << val[j] << ' ' << "inf" << '\n';
        }
    }
    ofs.close();
}


// sparse version 


void write_barcode_from_reduced_sparse(const std::vector<std::array<std::size_t,2>>& entries,
                                    const std::vector<int>& dim,
                                    const std::vector<float>& val,
                                    const std::string& out_path)
{
    if (entries.empty()) { std::ofstream(out_path).close(); return; }

    // infer n 
    std::size_t n = 0;
    for (const auto& e : entries) n = std::max(n, std::max(e[0], e[1]) + 1);

    // build per-column row lists (sorted)
    std::vector<std::vector<std::size_t>> cols(n);
    for (const auto& e : entries) cols[e[1]].push_back(e[0]);
    for (auto& c : cols) std::sort(c.begin(), c.end());

    std::vector<bool> row_used_as_pivot(n, false);

    std::ofstream ofs(out_path);
    ofs.setf(std::ios::fixed, std::ios::floatfield);
    ofs.precision(6);

    // paired columns : finite intervals
    for (std::size_t j = 0; j < n; ++j) {
        if (!cols[j].empty()) {
            std::size_t i = cols[j].back();  // lowest = largest row
            row_used_as_pivot[i] = true;
            int k = dim[j] - 1;
            ofs << k << ' ' << val[i] << ' ' << val[j] << '\n';
        }
    }
    // zero columns not used as pivot-row : infinite intervals
    for (std::size_t j = 0; j < n; ++j) {
        if (cols[j].empty() && !row_used_as_pivot[j]) {
            int k = dim[j];
            ofs << k << ' ' << val[j] << ' ' << "inf" << '\n';
        }
    }
    ofs.close();
}

int main(){
    auto F = read_filtration("./src/filtration.txt");
    size_t n = F.size();
    std::vector<int> dims(n);
    std::vector<float> vals(n);
    for(size_t i=0; i<n; i++){
        dims[i] = F[i].dim;
        vals[i] = F[i].val;
    }
    auto B = boundary_matrix_dense(F);
    gaussian_elim_dense(B);
    write_barcode_from_reduced_dense(B, dims, vals, "./barcode_dense.txt");

    auto B_sparse = boundary_matrix_sparse(F);
    auto B_sparse_reduced = gaussian_elim_sparse(B_sparse);
    write_barcode_from_reduced_sparse(B_sparse_reduced, dims, vals, "./barcode_sparse.txt");
}