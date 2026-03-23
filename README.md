# Finite Element Codes
Codes in Matlab and Python based on Finite Element Methods

---

# Buckling Analysis of a Two-Element Beam

This Python script performs a **linear buckling analysis** of a 2D beam composed of two Euler-Bernoulli elements with circular hollow cross-sections. It computes the **critical buckling load** (\(P_{cr}\)) using the **eigenvalue problem** formulation:

\[
\mathbf{K_f} \mathbf{v} = \lambda \mathbf{K_g} \mathbf{v}
\]

where:  
- \(\mathbf{K_f}\) is the **global stiffness matrix**,  
- \(\mathbf{K_g}\) is the **global geometric stiffness matrix**,  
- \(\lambda\) are the **eigenvalues**, the smallest of which corresponds to the **critical buckling load**.  

The code allows for varying material and geometric properties along elements, and applies boundary conditions by masking fixed degrees of freedom.

---

## Features

- Supports multiple elements with variable diameters and hollow sections.  
- Computes local and global stiffness matrices.  
- Solves the generalized eigenvalue problem to find the critical load.  
- Efficient use of `numpy` and `scipy.linalg.eig`.  

---

## Problem Illustration

![Buckling Beam Illustration](./a_technical_illustration_illustration_of_a_bucklin.png)

*Figure: Buckling mode of a two-element cantilever beam with different hollow circular cross-sections. Node numbers, cross-sections, and axial load \(P\) are shown.*

---

## Usage

Run the script using Python:

bash
python beam_buckling.py
