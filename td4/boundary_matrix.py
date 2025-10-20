from read_filtration import Simplex, read_filtration
from tqdm import tqdm
import sys

class SparseMatrix(set):
    def __init__(self, ens, m):
        super().__init__(ens)
        self.m = m

    def cols(self):
        return {j for (i,j) in self}
    def lines(self):
        return {i for (i,j) in self}

    def line(self,i):
        return super().intersection(set([(i,j) for j in range(len(self))]))
    
    def column(self,j):
        return super().intersection(set([(i,j) for i in range(len(self))]))

    def __str__(self):
        return super().__str__()
    def __repr__(self):
        return super().__repr__()
    def __iter__(self):
        return super().__iter__()
    def __contains__(self, item):
        return super().__contains__(item)
    def __len__(self):
        return self.m
    def add(self, element):
        super().add(element)
    def symmetric_difference(self, s):
        return SparseMatrix(super().symmetric_difference(s),len(self))
    
    def print(self):
        print_dense_boundary_matrix(sparse_to_dense(self))

def build_dense_boundary_matrix(filtration: list[Simplex]) -> list[list[int]]:
    m = len(filtration)
    boundary_matrix = [[0 for _ in range(m)] for _ in range(m)]
    for i in range(m):
        simplex_i = filtration[i]
        for j in range(i+1,m):
            simplex_j = filtration[j]
            if simplex_j.dim == simplex_i.dim + 1 and simplex_j.vert >= simplex_i.vert:
                boundary_matrix[i][j] = 1
    return boundary_matrix

def build_sparse_boundary_matrix(filtration: list[Simplex]) -> SparseMatrix:
    m = len(filtration)
    boundary_matrix = SparseMatrix(set(), m)
    for i in tqdm(range(m)):
        simplex_i = filtration[i]
        for j in range(i+1,m):
            simplex_j = filtration[j]
            if simplex_j.dim == simplex_i.dim + 1 and simplex_j.vert >= simplex_i.vert:
                boundary_matrix.add((i,j))
    return boundary_matrix


def print_dense_boundary_matrix(matrix: list[list[int]]):
    for row in matrix:
        print(" ".join(str(val) for val in row))

def sparse_to_dense(sparse_matrix:SparseMatrix) -> list[list[int]]:
    dense_matrix = [[0 for _ in range(len(sparse_matrix))] for _ in range(len(sparse_matrix))]
    for p in sparse_matrix:
        dense_matrix[p[0]][p[1]] = 1
    return dense_matrix

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Syntax: {sys.argv[0]} <filtration_file>")
        sys.exit(0)

    filtration_file = sys.argv[1]
    print("Reading filtration...")
    filtration = read_filtration(filtration_file)
    print("Done.")

    print("Boundary Matrix:")
    boundary_matrix = build_dense_boundary_matrix(filtration)
    print_dense_boundary_matrix(boundary_matrix)

    print("Sparse Boundary Matrix:")
    sparse_boundary_matrix = build_sparse_boundary_matrix(filtration)
    print(sparse_boundary_matrix)

    print_dense_boundary_matrix(sparse_to_dense(sparse_boundary_matrix))

