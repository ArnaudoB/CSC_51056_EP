def boundary_matrix(l_simplex):
    """
    Dense representation (TODO: sparse) over Z/2Z of the boundary matrix of the
    simplices in l_simplex.
    """
    n_simplices = len(l_simplex)
    matrix = [[0 for _ in range(n_simplices)] for _ in range(n_simplices)]
    for i in range(n_simplices):
        for j in range(i+1, n_simplices):
            if l_simplex[i].vert <= l_simplex[j].vert and l_simplex[i].dim + 1 == l_simplex[j].dim:
                matrix[i][j] = 1
    return matrix

if __name__ == "__main__":
    from read_filtration import read_filtration
    F = read_filtration("./filtration.txt")
    BM = boundary_matrix(F)
    for row in BM:
        print(row)
