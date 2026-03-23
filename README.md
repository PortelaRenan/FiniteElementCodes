# Finite Element Codes
Codes in Matlab and Python based on Finite Element Methods


Code: Euler Critical Load  - (Buckling Analysis of a Two-Element Beam Using Euler-Bernoulli Theory)

Description:
This Python script performs a linear buckling analysis of a 2D beam composed of two Euler-Bernoulli elements with circular hollow cross-sections. It computes the critical buckling load (𝑃_cr) using the eigenvalue problem formulation:

𝐾_𝑓 𝑣 = 𝜆 𝐾_𝑔

where:

* 𝐾_𝑓  is the global stiffness matrix,
𝐾
𝑔
K
g
	​

 is the global geometric stiffness matrix,
𝜆
λ are the eigenvalues, the smallest of which corresponds to the critical buckling load.

The code allows for varying material and geometric properties along elements, and applies boundary conditions by masking fixed degrees of freedom.

Features:

Supports multiple elements with variable diameters and hollow sections.
Computes local and global stiffness matrices.
Solves the generalized eigenvalue problem to find the critical load.
Efficient use of numpy and scipy.linalg.eig.

Usage Example:

python beam_buckling.py

Output:

Pcr = 12345.67 N

Dependencies:

numpy
scipy
