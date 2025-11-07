# Eigenvalue Problem for the Quantum Harmonic Oscillator (QHO)

This project is from independent research and development conducted during my **Master of Science program in Computational Physics**. The entire methodology, analysis, and implementation of all numerical solvers were executed solely by me, and the full work is documented in the accompanying report and code.

---

## Table of Contents

* [About The Project](#about-the-project)
* [Core Objectives](#core-objectives)
* [Languages and Libraries](#languages-and-libraries)
* [Methods Implemented](#methods-implemented)
* [Key Findings](#key-findings)
* [Getting Started](#getting-started)
* [Full Project Report](#full-project-report)
* [Contact](#contact)

---

## About The Project

This work focuses on numerically solving the **time-independent Schrödinger equation** for the one-dimensional Quantum Harmonic Oscillator (QHO). It serves as a comprehensive comparison of two key numerical methods: the **Finite Difference Method (FDM)** diagonalization and the **Inverse Power Iteration (IPI)**, evaluating their accuracy, stability, and computational efficiency.

The system's behavior is analyzed in both its unperturbed state (compared to analytical solutions) and a scenario where a perturbation is introduced.

## Core Objectives

1.  **Solve the QHO Eigenvalue Problem:** Compute the energy eigenvalues ($E_n$) and corresponding eigenfunctions ($\psi_n(x)$) using numerical methods.
2.  **Accuracy Assessment:** Compare numerical eigenvalues and eigenvectors against the exact **algebraic solutions**.
3.  **Efficiency Comparison:** Benchmark the runtime performance of **FDM diagonalization** versus the **IPI method with a shift** across varying matrix/grid sizes.
4.  **Perturbation Analysis:** Investigate the effect of introducing a time-dependent **Symmetric Gaussian potential** on the system's energy levels and wave functions.

---

## Languages and Libraries

| Category | Tools & Libraries | Competency Demonstrated |
| :--- | :--- | :--- |
| **Language** | Python | Efficient development and handling of complex number mathematics. |
| **Numerical** | NumPy, SciPy | Advanced array manipulation and solving systems of linear equations for the TISE. |
| **Visualization** | Matplotlib | Generating high-quality physics focused data plots. |

---

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

---

## Getting Started

### Execution

To run the simulation and generate the results and visualizations, execute the core solver script:

```bash
python QHO_FDM_IPI.py 
```
---

## Full Project Report

For a complete breakdown of the theoretical derivations and full results, please see the final project report:

[**Full QHO Project Report (PDF)**](Eigenvalue_FDM_IPI.pdf)

---

## Contact

I'm happy to hear your feedback or answer any questions about this project!

**Author** Rama Khalil 

**Email**  rama.khalil.990@gmail.com
