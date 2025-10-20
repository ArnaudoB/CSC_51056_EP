#include <vector>
#include <string>
#include <iostream>
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


void write_barcode_from_reduced_sparse(const std::vector<std::vector<int>>& cols,
                                       const std::vector<int>& dim,
                                       const std::vector<float>& val,
                                       const std::string& out_path)
{
    const int n = cols.size();

    if (n == 0) { std::ofstream(out_path).close(); return; }

    std::vector<bool> row_used_as_pivot(n, false);

    std::ofstream ofs(out_path);
    ofs.setf(std::ios::fixed, std::ios::floatfield);
    ofs.precision(6);

    // paired columns (finite intervals)
    for (int j = 0; j < n; ++j) {
        if (!cols[j].empty()) {
            const int i = cols[j].back();      // lowest = largest row index
            row_used_as_pivot[i] = true;

            const int k = dim[j] - 1;          // homology dimension for the pair
            ofs << k << ' ' << val[i]
                << ' ' << val[j] << '\n';
        }
    }

    // unpaired zero columns (infinite intervals)
    for (std::size_t j = 0; j < n; ++j) {
        if (cols[j].empty() && !row_used_as_pivot[j]) {
            const int k = dim[j];
            ofs << k << ' ' << val[j] << ' ' << "inf" << '\n';
        }
    }
    ofs.close();
}

int run_barcode(const std::string& filtration_path, const std::string& output_path) {
    std::cout << "Reading filtration..." << std::endl;
    auto F = read_filtration(filtration_path);
    std::cout << "Done." << std::endl;

    std::cout << "Parsing dimensions and values..." << std::endl;
    size_t n = F.size();
    std::vector<int> dims(n);
    std::vector<float> vals(n);
    for (size_t i = 0; i < n; i++) {
        dims[i] = F[i].dim;
        vals[i] = F[i].val;
    }
    std::cout << "Done." << std::endl;

    std::cout << "Computing the boundary matrix..." << std::endl;
    auto B_sparse = boundary_matrix_sparse_fast(F);
    std::cout << "Done." << std::endl;

    std::cout << "Reducing the boundary matrix..." << std::endl;
    auto B_sparse_reduced = gaussian_elim_sparse(B_sparse);
    std::cout << "Done." << std::endl;

    std::cout << "Computing the barcodes..." << std::endl;
    write_barcode_from_reduced_sparse(B_sparse_reduced, dims, vals, output_path);
    std::cout << "Done." << std::endl;
    return 0;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <filtration.txt> <output.txt>\n";
        return 1;
    }
    return run_barcode(argv[1], argv[2]);
}