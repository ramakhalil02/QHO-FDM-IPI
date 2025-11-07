# Eigenvalue Problem for the Quantum Harmonic Oscillator (QHO)

This project is from independent research and development conducted during my **Master of Science program in Computational Physics**. The entire methodology, analysis, and implementation of all numerical solvers were executed solely by me, and the full work is documented in the accompanying report and code.

---

## Overview

This work focuses on numerically solving the **time-independent Schrödinger equation** for the one-dimensional Quantum Harmonic Oscillator (QHO). It serves as a comprehensive comparison of two key numerical methods: the **Finite Difference Method (FDM)** diagonalization and the **Inverse Power Iteration (IPI)**, evaluating their accuracy, stability, and computational efficiency.

The system's behavior is analyzed in both its unperturbed state (compared to analytical solutions) and a scenario where a perturbation is introduced.

## Core Objectives

1.  **Solve the QHO Eigenvalue Problem:** Compute the energy eigenvalues ($E_n$) and corresponding eigenfunctions ($\psi_n(x)$) using numerical methods.
2.  **Accuracy Assessment:** Compare numerical eigenvalues and eigenvectors against the exact **algebraic solutions**.
3.  **Efficiency Comparison:** Benchmark the runtime performance of **FDM diagonalization** versus the **IPI method with a shift** across varying matrix/grid sizes.
4.  **Perturbation Analysis:** Investigate the effect of introducing a time-dependent **Symmetric Gaussian potential** on the system's energy levels and wave functions.

## Methods Implemented

The computational core of the project, found in the `QHO_FDM_IPI.py` script, implements the following techniques:

| Method | Role in Project | Key Implementation Detail |
| :--- | :--- | :--- |
| **Finite Difference Method (FDM)** | Constructs the pentadiagonal Hamiltonian matrix. | Utilizes a **five-point stencil** to approximate the second derivative, achieving $\mathcal{O}(h^4)$ accuracy. Full matrix diagonalization is used to find all states. |
| **Inverse Power Iteration (IPI)** | Targets specific eigenvalues closest to a given shift $\xi$. | Solves the generalized eigenvalue problem using the implicit **shift-invert method**, providing high efficiency for targeted excited states. |
| **Perturbation Model** | Dynamic system analysis. | A **Symmetric Gaussian potential** is added to the QHO potential to model a localized time-dependent interaction. |

## Key Findings

* **High Accuracy:** Both FDM and IPI achieved excellent agreement with the analytical QHO solutions, with small deviations typically limited to higher-energy states.
* **Runtime Efficiency:** The **Inverse Power Iteration (IPI) method scaled significantly better** than the full FDM diagonalization, proving far more efficient for large systems when only a subset of eigenvalues is required (as detailed in Section 4.3 of the report).
* **Perturbation Impact:** The introduction of the symmetric Gaussian potential primarily affected the **symmetric states** (even quantum numbers $n = 0, 2, 4, \dots$) due to their high probability density at the perturbation's center (refer to Figure 4 in the report).

## Project Files

* **`QHO_FDM_IPI.py`**: Python script containing the numerical implementation of the FDM and IPI solvers using `numpy` and `scipy`.
* **`Eigenvalue_FDM_IPI.pdf`**: The final comprehensive report containing the theoretical background, detailed derivations, numerical results (tables and figures), and discussion.
