#!/usr/bin/env python3

import sys

class Simplex:

    def __init__(self, dim: int, val: float, vert: set[int]):
        self.dim = dim
        self.val = val
        self.vert = vert

    def __repr__(self):
        return str((self.dim, self.val, self.vert))
    
    def __eq__(self, value):
        return self.vert == value.vert

def read_filtration(name: str):
    F = []
    try:
        with open(name, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                # structure attendue: val dim v0 v1 ... v_dim
                if len(parts) < 2:
                    continue
                val = float(parts[0])
                dim = int(parts[1])
                if dim == -1:
                    continue
                expected = 2 + (dim + 1)
                if len(parts) < expected:
                    # ligne incomplète ignorée
                    continue
                verts = {int(x) for x in parts[2:2 + dim + 1]}
                F.append(Simplex(dim=dim, val=val, vert=verts))
    except OSError:
        print(f"Failed to read file {name}")
    F = sorted(F, key=lambda s: (s.val, list(s.vert), s.dim))
    return F


def main(argv) -> int:
    if len(argv) != 2:
        print(f"Syntax: {argv[0]} <filtration_file>")
        return 0

    name = argv[1]
    print("Reading filtration...")
    F = read_filtration(name)
    print("Done.")

    for s in F:
        verts_sorted = sorted(s.vert)
        verts_str = ", ".join(str(v) for v in verts_sorted)
        print(f"{{val={s.val}; dim={s.dim}; [{verts_str}]}}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
