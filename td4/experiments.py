from read_filtration import *
from boundary_matrix import *
from reduction import *

def aux_simplex_d_sphere(d):
    if d == 2:
        cercle = [[Simplex(0,0,{1}), Simplex(0,0,{2}),Simplex(0,1,{3})],
                  [Simplex(1,1,{1,2}), Simplex(1,1,{2,3}),Simplex(1,1,{1,3})]]
        return cercle
    elif d > 2:
        sphere = aux_simplex_d_sphere(d-1)
        lenghts = [len(dim) for dim in sphere]
        for k in range(d-1):
            if k == 0: 
                sphere[0].append(Simplex(0,0,{1+len(sphere[0])}))
            else:
                to_add = [Simplex(k,d-1,sphere[k-1][i].vert.union(sphere[k-1][j].vert)) 
                                    for i in range(lenghts[k-1],len(sphere[k-1])) for j in range(i)]
                sphere[k].extend(to_add)
        
        sphere.append([Simplex(k,d-1,{j for j in range(d+1) if j!=i}) for i in range(d+1)])

        return sphere
    
def simplex_d_sphere(d):
    sphere = aux_simplex_d_sphere(d)
    return [s for dim_simplices in sphere for s in dim_simplices]

def filtration_d_sphere(d):
    F = simplex_d_sphere(d)
    F = sorted(F, key=lambda s: (s.val, list(s.vert), s.dim))

if __name__ == "__main__":

    for d in aux_simplex_d_sphere(4):
        print(d)
