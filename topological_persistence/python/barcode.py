from boundary_matrix import boundary_matrix
from gaussian_elimination import gaussian_elimination

def barcode(l_simplex):
    b_matrix = boundary_matrix(l_simplex)
    e_matrix = gaussian_elimination(b_matrix)
    n = len(e_matrix)
    bars = []
    low = [-1 for _ in range(n)]
    for j in range(n):
        for i in reversed(range(n)):
            if e_matrix[i][j] != 0:
                low[j] = i
                break
    paired = set()
    for j in range(n):
        if low[j] != -1:
            bars.append((low[j], j))
            paired.add(low[j])
    for i in range(n):
        if i not in paired:
            bars.append((i, None))
    return bars