# Truss: 1D Finite Element Analysis of a Tapered Bar

A small Python script that solves an axially loaded, linearly tapered bar with the finite element method (FEM), compares the result against the closed-form solution, and studies how the error shrinks as the mesh is refined.

## Problem description

A bar of length `L` has a cross-section that grows linearly along its length:

```
A(x) = A0 * (1 + x/L),     0 <= x <= L
```

- **x = 0** is the free end, where the axial load `P` is applied.
- **x = L** is the fixed end (displacement = 0).
- Since `P` is negative, the load points in the -x direction, pulling the bar away from the support, so the bar is in tension.

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Young's modulus | `E` | 70 GPa |
| Area at free end | `A0` | 12.5 cm² (12.5e-4 m²) |
| Length | `L` | 0.5 m |
| Applied load | `P` | -5 kN |

### Exact solution

Because the axial force is constant along the bar, the exact fields are:

```
u(x)     = -(|P| L / (E A0)) * ln( 2 / (1 + x/L) )
sigma(x) =   |P| / ( A0 (1 + x/L) )
```

These are used as the reference for all error measurements.

## Method

1. **Mesh**: the bar is divided into `n` equal 2-node bar elements (`n + 1` nodes).
2. **Element stiffness**: each element uses the cross-sectional area evaluated at its midpoint, so

   ```
   k_e = (E * A_e / L_e) * [[ 1, -1],
                            [-1,  1]]
   ```

3. **Assembly**: element matrices are summed into the global stiffness matrix.
4. **Loads and boundary conditions**: `P` is applied at node 0, and the last node (x = L) is fixed.
5. **Solve**: the reduced system `K_ff u_f = F_f` is solved with `numpy.linalg.solve`.
6. **Post-processing**: stress is constant within each element, `sigma_e = E (u2 - u1) / L_e`.

The script repeats this for meshes of **1, 2, 4, 8, 16, 32 and 64 elements**.

## Requirements

- Python 3.8+
- [NumPy](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/)

```bash
pip install numpy matplotlib
```

## Usage

Run the script from this folder (replace `truss.py` with the actual file name):

```bash
python truss.py
```

Two figures are displayed; nothing is written to disk by default.

## Output

### Figure 1: Three-panel comparison (meshes of 1, 2, 4, 8, 16 elements)

![Displacement error, stress and stress error for several meshes](images/mesh_comparison.png)

| Panel | Content |
|-------|---------|
| Left | Displacement error (FE - exact), in µm. It is zero at the fixed end and shrinks quickly as the mesh is refined. |
| Centre | Stress as a staircase per element, with the exact curve (dashed black). The steps converge onto the exact curve. |
| Right | Stress error (FE - exact), in MPa. The error is largest at element ends and is sawtooth-shaped, shrinking with each refinement. |

### Figure 2: Convergence plot

![Convergence of the free-end displacement error](images/convergence.png)

Log-log plot of the free-end displacement error `|u_FE(0) - u_exact(0)|` against the number of elements (1 to 64).

| Elements | Free-end error (m) |
|----------|--------------------|
| 1 | 7.57e-07 |
| 2 | 2.12e-07 |
| 4 | 5.51e-08 |
| 8 | 1.39e-08 |
| 16 | 3.48e-09 |
| 32 | 8.72e-10 |
| 64 | 2.18e-10 |

The fitted slope is about -2, so the free-end displacement converges at second order: doubling the number of elements reduces the error by roughly a factor of four.

### Saving the figures

To regenerate the images, add these lines before each figure is closed or shown (the first after `plt.tight_layout()`, the second after the convergence plot's `plt.title(...)`):

```python
plt.savefig('images/mesh_comparison.png', dpi=150)   # Figure 1
plt.savefig('images/convergence.png', dpi=150)       # Figure 2
```

## Code structure

| Section | Purpose |
|---------|---------|
| Parameters | Material, geometry and load definitions |
| `uExact`, `sExact` | Closed-form displacement and stress |
| Mesh loop | Builds, assembles and solves the FE model for each mesh size |
| `results` | Dictionary `results[n] = (x, u, stress)` holding each mesh's solution |
| `freeEndError` | Free-end displacement error for each mesh, used in the convergence plot |
| Plots | Error and stress comparison, then convergence study |

## Customising

- **Geometry, material, load**: edit `E`, `A`, `L`, `P` at the top of the script.
- **Mesh sizes**: change the `meshes` list (the plotted subset is `plotMeshes`).
- **Different taper**: update `Ae`, `uExact` and `sExact` consistently.
- **Different boundary conditions**: modify the `bc` list and the load vector `F`.

## Notes

- Element area is taken at the element midpoint, which is a simple and accurate approximation for a linear taper.
- Stress is piecewise constant, so stress error is largest at element ends and decreases as the mesh is refined.
- The stiffness matrix is stored as a dense array, which is fine for this small problem but would need sparse storage for large meshes.
