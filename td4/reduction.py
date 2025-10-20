from read_filtration import *
from boundary_matrix import *
import sys
from math import inf
from tqdm import tqdm

def dense_low(matrix : list[list[int]], j : int) -> int:
    low = 0
    
    for i in range(len(matrix)):
        if matrix[i][j] != 0:
            low = i
    return low


def dense_addMod2(tab1, tab2):
    assert len(tab1) == len(tab2)
    n = len(tab1)
    res = [0]*n
    for i in range(n):
        res[i] = (tab1[i] + tab2[i]) % 2

def sparse_addMod2(tab1 : SparseMatrix,tab2 : set) -> SparseMatrix:
    diff = tab1.symmetric_difference(tab2)
    return diff



def sparse_low(matrix: SparseMatrix, j: int) -> int:
    col = matrix.column(j)
    if len(col) == 0:
        return -1
    else:
        return max(col)[0]
    

def value_that_is_lower(tab,val):
    for el in tab:
        if el < val:
            return el
    return None


def sparse_gaussian_reduction(matrix: SparseMatrix, verbose = False) -> SparseMatrix:
    m  = len(matrix)
    low_to_col = [set() for _ in range(m)]
    for j in range(m):
        low_to_col[sparse_low(matrix,j)].add(j)

    for j in tqdm(matrix.cols()):
        while True :
            l = sparse_low(matrix,j)
            if l >= 0 and len(low_to_col[l]) > 1: 
                i = value_that_is_lower(low_to_col[l],j)
                if i is None:
                    break
                else:
                    matrix = sparse_addMod2(matrix,{(k,j) for (k,_) in matrix.column(i)})
                    print(f"{j} <- {j} + {i} ({i} < {j})")
            else:
                break
            if verbose:
                print("--------------------")
                matrix.print()
                print("--------------------")
            low_to_col[l].remove(j)
            low_to_col[sparse_low(matrix,j)].add(j)

    return matrix

def dense_gaussian_reduction(matrix : list[list[int]]) ->list[list[int]]:
    m  = len(matrix)
    low_to_col = {}
    for j in range(m):
        low_to_col[dense_low(matrix,j)] = j

    for j in range(m):
        while True :
            l = dense_low(matrix,j)
            if l > 0 and l in low_to_col:
                i = low_to_col[l]
                matrix[:][j] = dense_addMod2(matrix[:][j],matrix[:][i])
            else:
                break
        for i in range(m):
            low_to_col[dense_low(matrix,i)] = i
    return matrix


def sparse_barcode_builder(matrix : SparseMatrix, filtration:list[Simplex]) -> None:
    m = len(matrix)
    sides = set(range(m))
    barcode = []
    for j in tqdm(matrix.cols()):
        i = sparse_low(matrix,j)
        sides = sides.difference({i,j})
        barcode.append((filtration[i].dim,filtration[i].val,filtration[j].val))
    
    for i in sides:
        barcode.append((filtration[i].dim,filtration[i].val,inf))

    return sorted(barcode)

def print_text_barcode(barcode):
    for bar in barcode:
        print(f"{bar[0]} {bar[1]} {bar[2]}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Syntax: {sys.argv[0]} <filtration_file>")
        sys.exit(0)

    filtration_file = sys.argv[1]
    print("Reading filtration...")
    filtration = read_filtration(filtration_file)
    #print(filtration)
    print("Done.")

    print("Building boundary matrix...")
    boundary_matrix = build_sparse_boundary_matrix(filtration)
    #print(boundary_matrix)
    print("Done.")

    print("Boundary Matrix:")
    print_dense_boundary_matrix(sparse_to_dense(boundary_matrix))

    print("Reducing Boundary Matrix:")
    boundary_matrix = sparse_gaussian_reduction(boundary_matrix)
    print("Reduced Boundary Matrix:")
    boundary_matrix.print()

    print("computing barcode")
    barcode = sparse_barcode_builder(boundary_matrix,filtration)

    print(barcode)