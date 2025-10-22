// includes

#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <array>
#include <limits>
#include <chrono>

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

// struct definitions


typedef std::set<int> vertices;

struct simplex{
  int dim;
  float val;
  vertices vert;
};

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

// filtration reading

std::vector<simplex> read_filtration(std::string name){

  std::vector<simplex> F; F.clear();
  char input_file[100];
  sprintf(input_file, "%s", (char*) name.c_str());
  std::ifstream input(input_file);

  if (input){
    std::string line;
    while(getline(input,line)){
      std::stringstream stream(line);
      simplex s; s.vert.clear();
      s.dim = -1; stream >> s.val; stream >> s.dim; int i = 0;
      while(i <= s.dim){
        int f; stream >> f;
	s.vert.insert(f); i++;
      }
      if(s.dim != -1)
        F.push_back(s);
    }
  }
  else{std::cout << "Failed to read file " << name << std::endl;}
  std::sort(F.begin(), F.end(), [](const simplex& a, const simplex& b) {
    return std::tie(a.val, a.dim, a.vert) < std::tie(b.val, b.dim, b.vert);});
  return F;
};

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

// boundary matrix calculation

std::vector<std::vector<int>> boundary_matrix_dense_fast(const std::vector<simplex>& simplices) {
    const size_t n = simplices.size();

    // boundary_matrix initialization
    std::vector<std::vector<int>> boundary_matrix(n, std::vector<int>(n, 0));

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
                boundary_matrix[it->second][j] = 1;
            }
        }
    }

    return boundary_matrix;
};

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
};

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

// gaussian elimination

    // helper : symmetric difference of two sorted vectors

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

std::vector<std::vector<int>> gaussian_elim_sparse(std::vector<std::vector<int>> cols)
{
    const size_t n = cols.size();
    if (n == 0) return cols;

    const int NONE = -1;
    std::vector<int> pivot_of_row(n, NONE); // pivot_of_row[r] = column owning pivot at row r

    for (size_t j = 0; j < n; ++j) {
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

    return cols;
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

// barcode

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
    ofs.precision(17);

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

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

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

    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Computing the boundary matrix..." << std::endl;
    auto B_sparse = boundary_matrix_sparse_fast(F);
    std::cout << "Done." << std::endl;

    std::cout << "Reducing the boundary matrix..." << std::endl;
    auto B_sparse_reduced = gaussian_elim_sparse(B_sparse);
    std::cout << "Done." << std::endl;

    std::cout << "Computing the barcodes..." << std::endl;
    write_barcode_from_reduced_sparse(B_sparse_reduced, dims, vals, output_path);
    std::cout << "Done." << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
    return 0;
}

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

// d-sphere filtration

static inline bool less_lex_vertices(const vertices& A, const vertices& B) {
    return std::lexicographical_compare(A.begin(), A.end(), B.begin(), B.end());
}
struct FiltrationLess {
    bool operator()(const simplex& a, const simplex& b) const {
        if (a.val != b.val) return a.val < b.val;
        if (a.dim != b.dim) return a.dim < b.dim;
        return less_lex_vertices(a.vert, b.vert);
    }
};

std::vector<simplex> filtration_sphere(int d) {
    const int n = d + 1;                // number of vertices
    std::vector<int> t(n + 1, 0);
    for (int v = 3; v <= n; ++v) t[v] = v - 2; // vertex times: t[1]=0, t[2]=0, t[3]=1, t[4]=2, ..., t[n]=d-1

    std::vector<simplex> F;
    const int FULL = (1 << n) - 1;
    for (int mask = 1; mask < FULL; ++mask) {
        vertices V;
        int max_time = 0, cnt = 0;
        for (int j = 0; j < n; ++j) if (mask & (1 << j)) {
            int v = j + 1; V.insert(v); ++cnt;
            max_time = std::max(max_time, t[v]);
        }
        F.push_back(simplex{cnt - 1, static_cast<float>(max_time), std::move(V)});
    }
    std::sort(F.begin(), F.end(), FiltrationLess{});
    return F;
}

void save_sphere_filtration(int d, const std::string& path) {
    auto F = filtration_sphere(d);
    std::ofstream out(path);
    for (const auto& s : F) {
        out << s.val << ' ' << s.dim;
        for (int v : s.vert) out << ' ' << v;
        out << '\n';
    }
}

int run_barcodes_sphere(int max_d) {
    if (max_d < 1) {
        std::cerr << "max_d must be >= 1\n";
        return 1;
    }

    for (int i = 1; i <= max_d; ++i) {
        const std::string fil =
            std::string("./src/filtrations/filtration_") + std::to_string(i) + "_sphere.txt";
        const std::string out =
            std::string("./output/barcode_") + std::to_string(i) + "_sphere.txt";

        // build filtration for S^{i-1} : this is not optimal because we recompute every time the filtrations of S^{j}, j<i when computing the filtration of S^{i} but it's fine
        save_sphere_filtration(i, fil);

        // compute barcode from that filtration
        int rc = run_barcode(fil, out);
        if (rc != 0) {
            std::cerr << "Failed on d=" << i << " (rc=" << rc << ")\n";
            return rc;
        }
    }
    return 0;
};

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

int main(int argc, char** argv) {
    if (argc == 3) {
        return run_barcode(argv[1], argv[2]);
    }

    if (argc == 2) {
        try {
            int max_d = std::stoi(argv[1]);
            return run_barcodes_sphere(max_d);
        } catch (const std::exception& e) {
            std::cerr << "Invalid integer for max_d: " << argv[1] << "\n";
            return 1;
        }
    }

    std::cerr
        << "Usage:\n"
        << "  " << argv[0] << " <filtration.txt> <output.txt>\n"
        << "    Compute barcode for an existing filtration file.\n\n"
        << "  " << argv[0] << " <max_d>\n"
        << "    Generate filtrations and barcodes for spheres S^{0}..S^{max_d-1}\n"
        << "    (i.e., build boundary of a d-simplex for d=1..max_d).\n";
    return 1;
}