from math import inf
import matplotlib.pyplot as plt

class bar():
    def __init__(self, dim, birth, death):
        self.dim = dim
        self.birth = birth
        self.death = death

    def __repr__(self):
        return f"({self.dim}, {self.birth}, {self.death})"
    
def plot_barcode(barcode:list[bar],save_path: str = None):

    max_death = max([bar.death for bar in barcode if bar.death != inf] + [bar.birth for bar in barcode])
    y = 0
    # Group bars by dimension
    bars_by_dim = {}
    for b in barcode:
        if b.dim not in bars_by_dim:
            bars_by_dim[b.dim] = []
        bars_by_dim[b.dim].append(b)
    
    # Create subplots
    n_dims = len(bars_by_dim)
    fig, axes = plt.subplots(n_dims, 1, figsize=(10, 3*n_dims), sharex=True)
    if n_dims == 1:
        axes = [axes]
    
    # Plot each dimension
    for idx, (dim, bars) in enumerate(sorted(bars_by_dim.items(),reverse=True)):
        ax = axes[idx]
        y_local = 0
        for b in bars:
            if b.death == inf:
                ax.hlines(y_local, b.birth, max_death*1.3, colors='b')
                ax.plot([b.birth], [y_local], 'bo')
                ax.arrow(max_death*1.3, y_local, 0.5, 0, head_width=0.3, head_length=0.2, fc='b', ec='b')
            else:
                ax.hlines(y_local, b.birth, b.death, colors='b')
                ax.plot([b.birth, b.death], [y_local, y_local], 'bo')
            y_local += 1
        ax.set_ylim(-1, y_local)
        ax.set_ylabel(f"Dim {dim}")
        ax.set_yticks([])
    
    axes[-1].set_xlim(0, max_death*1.3+2)
    axes[-1].set_xlabel("Filtration Value")
    fig.suptitle("Barcode", fontsize=16)
    plt.subplots_adjust(hspace=0)
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()

def read_barcode_from_file(filename: str) -> list[bar]:
    barcode = []
    with open(filename, 'r') as f:
        for line in f:
            dim, birth, death = line.strip().split()
            dim = int(dim)
            birth = float(birth)
            if death == 'inf':
                death = inf
            else:
                death = float(death)
            barcode.append(bar(dim, birth, death))
    return barcode

if __name__ == "__main__":
    import pathlib
    barcode_dir = pathlib.Path("barcode")
    for filename in barcode_dir.glob("*.txt"):
        print(f"\nProcessing {filename}")
        barcode = read_barcode_from_file(filename)
        print(barcode)
        plot_barcode(barcode, save_path=f"barcode/{filename.stem}.png")
