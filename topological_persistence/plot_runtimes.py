import numpy as np
import matplotlib.pyplot as plt

d = np.arange(1, 21)  # dimensions
t_ms = np.array([
    0.1, 0.1, 0.1, 0.1, 0.2, 0.6, 0.9, 2, 7, 12,
    24, 54, 128, 332, 716, 1774, 3760, 8576, 21017, 46664
], dtype=float)

f_best = (d + 1) * (2.0 ** d)
f_worst = 8.0 ** d + d * 2.0 ** d
f_med = 4.0 ** d + d * 2.0 ** d

def fit_scale(f, y):
    num = np.dot(f, y)
    den = np.dot(f, f)
    return num / den

c1 = fit_scale(f_best, t_ms)
c2 = fit_scale(f_worst, t_ms)
c3 = fit_scale(f_med, t_ms)

t_best = c1 * f_best
t_worst = c2 * f_worst
t_med = c3 * f_med

# Plot
plt.figure(figsize=(7, 4.5))
plt.plot(d, t_ms, 'o-', linewidth=2, markersize=5, label='Measured runtime')
plt.plot(d, t_best, '--', linewidth=2, label=r'Theory (best case) $\; \propto (d+1)\,2^{d}$ (fitted)', c="green")
plt.plot(d, t_worst, '--', linewidth=2, label=r'Theory (worst case)$\; \propto d\,2^{d} + 8^{d}$ (fitted)', c="red")
plt.plot(d, t_med, '--', linewidth=2, label=r'Theory (medium case)$\; \propto d\,2^{d} + 4^{d}$ (fitted)', c="orange")

plt.title(r'Runtime vs. Dimension $d$ (Barcode computation of the $d$-sphere)')
plt.xlabel('Dimension $d$')
plt.ylabel('Time (ms)')
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.xticks(d)
plt.legend()
plt.tight_layout()
plt.savefig("./cpp/output/runtimes_d_sphere.png")
