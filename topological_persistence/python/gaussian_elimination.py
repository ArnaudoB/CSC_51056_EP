def gaussian_elimination(boundary_matrix):

    """
    Perform Gaussian elimination over Z/2Z on the given boundary matrix.
    Returns the reduced matrix and a list of pivot columns.

    Complexity is O(n^4)
    """
    n = len(boundary_matrix)
    matrix = [row[:] for row in boundary_matrix]

    for j in range(n):
        # find the lowest value in the jth column
        low_j = -1
        for i in reversed(range(n)) :
            if matrix[i][j] != 0:
                low_j = i
                break
        if low_j == -1:
            continue
        flag = True
        while flag:
            flag = False
            # find the lowest value in the ith column
            low_i = -1
            for i in range(j):
                for k in reversed(range(n)) :
                    if matrix[k][i] != 0:
                        low_i = k
                        break
            if low_i != -1 and matrix[i][low_i] == matrix[j][low_j]:
                # add the two columns mod 2
                for k in range(n):
                    matrix[k][j] = (matrix[k][j] + matrix[k][i])%2
                flag = True

if __name__ == '__main__':
    from read_filtration import read_filtration
    from boundary_matrix import boundary_matrix
    F = read_filtration("./filtration.txt")
    BM = boundary_matrix(F)
    for row in BM:
        print(row)
    elimination = gaussian_elimination(BM)
    for row in elimination:
        print(row)

