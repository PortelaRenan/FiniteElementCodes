# Finite Element Method in Python

Python implementations and worked solutions to selected problems from **Daryl L. Logan, *A First Course in the Finite Element Method*, 6th Edition** (Cengage Learning).

Every problem is solved with a small, readable FE library written from scratch (no black-box solvers), then **verified** against an analytical solution, the textbook's answer, or an independent solver. The goal is to show the full path from theory to working code: derivation, implementation, verification, and interpretation.

> **Disclaimer:** This is an independent study project. It is not affiliated with or endorsed by the author or publisher. Problem statements are paraphrased, and all solutions and code are my own work. Please buy the book to follow along.

---

## Highlights

- From-scratch FE library: bar, 2D/3D truss, beam, frame, CST and Q4 elements
- Hand derivations (scanned notes) alongside the code
- Weighted residual methods (collocation, Galerkin, least squares) compared against exact solutions
- Symmetry and boundary-condition handling demonstrated on real textbook problems
- Automated tests that check every solution against known results
- Plots of deformed shapes, stress distributions, and convergence studies

---

## Repository structure

```
.
├── fem/                    # Reusable FE library
│   ├── elements/           # bar.py, truss2d.py, truss3d.py, beam.py, frame.py, cst.py, q4.py
│   ├── assembly.py         # Global stiffness assembly
│   ├── boundary.py         # Constraints (elimination / penalty), symmetry
│   ├── loads.py            # Nodal forces, energy-equivalent distributed loads
│   ├── solver.py           # Linear solver and post-processing
│   └── plotting.py         # Deformed shape and stress plots
├── ch02_stiffness_method/
├── ch03_truss/
│   ├── p3_46/              # Each problem: README.md, solve.py, notes/
│   ├── p3_58/
│   └── p3_60_weighted_residuals/
├── ch04_beams/
├── ch05_frames_grids/
├── ch06_plane_stress_strain/
├── ...                     # Further chapters
├── tests/                  # pytest suite (analytical and textbook checks)
├── docs/                   # Theory notes and derivations
├── requirements.txt
└── README.md
```

Each problem folder contains:

| File | Purpose |
|---|---|
| `README.md` | Problem summary, method, results table, verification |
| `solve.py` | Script that reproduces the results |
| `notes/` | Hand derivations (scans or LaTeX) |

---

## Progress

Mark each chapter as work is completed. Chapter titles follow the book's table of contents; adjust if your edition differs.

| Ch. | Topic | Status |
|---|---|---|
| 2 | Introduction to the stiffness (displacement) method | ☐ |
| 3 | Development of truss equations (incl. weighted residual methods) | ☐ |
| 4 | Development of beam equations | ☐ |
| 5 | Frame and grid equations | ☐ |
| 6 | Plane stress and plane strain (CST) | ☐ |
| 7 | Practical considerations in modeling | ☐ |
| 8 | Linear-strain triangle | ☐ |
| 9 | Axisymmetric elements | ☐ |
| 10 | Isoparametric formulation (Q4) | ☐ |
| 11 | Three-dimensional stress analysis | ☐ |
| 13 | Heat transfer and mass transport | ☐ |
| 15 | Thermal stress | ☐ |
| 16 | Structural dynamics | ☐ |

---

## Featured problems

### Problem 3.46: Truss solved using symmetry
Plane truss with a central load. Only half the structure is modeled, using roller conditions on the symmetry axis and a halved load and member area on the axis member. The results are checked against full-model equilibrium and member compatibility.

### Problem 3.58: Energy-equivalent nodal forces
Consistent nodal loads for linearly and quadratically varying distributed loads, computed by integrating the shape functions. Verified by confirming that the nodal forces sum to the total applied load.

### Problem 3.60: Tapered bar with weighted residuals
A variable-area bar solved with a cubic polynomial using **collocation**, **Galerkin**, and **least squares**. Tip displacements are compared against the exact solution (u′ = P/(EA)), with error and convergence discussed.

*(Add further featured problems here, each with a one-line description and a result image.)*

---

## Results gallery

Add images here once generated:

| Deformed truss | Stress contour | Convergence study |
|---|---|---|
| `docs/img/truss.png` | `docs/img/stress.png` | `docs/img/convergence.png` |

---

## Installation

```bash
git clone https://github.com/<your-username>/<repo-name>.git
cd <repo-name>
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

**Requirements:** Python 3.10+, NumPy, SciPy, Matplotlib, SymPy, pytest.

## Usage

Run any problem script directly:

```bash
python ch03_truss/p3_46/solve.py
```

Minimal example using the library:

```python
from fem.elements.truss2d import Truss2D
from fem.assembly import assemble
from fem.solver import solve

# Define nodes, elements, supports, and loads, then:
K = assemble(nodes, elements)
U, reactions = solve(K, loads, fixed_dofs)
```

## Testing

```bash
pytest -v
```

Tests compare computed results with analytical solutions and textbook answers, and run automatically on every push via GitHub Actions.

---

## Verification approach

1. **Analytical:** closed-form solutions where available (bars, beams, tapered bars)
2. **Equilibrium:** global and nodal force balance checked in code
3. **Textbook answers:** compared to the book's published results
4. **Independent solver:** selected problems cross-checked with an open-source solver such as CalculiX or FEniCS

Each problem README states which checks were applied and the observed error.

---

## Roadmap

- [ ] Complete Chapters 2 to 6 problem sets
- [ ] Add convergence studies (h-refinement) for beam and plane stress elements
- [ ] Add modal analysis and time-dependent heat transfer
- [ ] Capstone: bracket with stress concentration and mesh sensitivity study
- [ ] Jupyter notebooks with narrated walkthroughs

---

## Reference

Logan, D. L. *A First Course in the Finite Element Method*, 6th ed., Cengage Learning.

## License

Code is released under the MIT License (see `LICENSE`). Book content, figures, and problem text remain the property of the publisher and are not reproduced here.

## Contact

**Your Name**: [LinkedIn](https://linkedin.com/in/your-profile) | [Email](mailto:you@example.com)

Feedback, corrections, and suggestions are welcome through issues or pull requests.
