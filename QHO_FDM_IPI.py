import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse import linalg as sla
import time
from matplotlib.lines import Line2D

# Constants
hbar = 1
m = 1
omega = 1
xmin = -10
xmax = 10 
Nx = 500
Neigs = 10
Nround = 8
A = 10.0
alpha = 5.0
tol = 1e-8
max_iter=10

def harmonic_potential(x, A, alpha, pert=False):
    V = x**2
    if pert:
        V += A * np.exp(-alpha * x**2)
    return V

def build_hamiltonian(xmin, xmax, Nx, A, alpha, pert=False):
    x = np.linspace(xmin, xmax, Nx)
    dx = x[1] - x[0]
    V = harmonic_potential(x, A, alpha, pert)
    H = sparse.lil_matrix((Nx, Nx))

    for i in range(2, Nx - 2):
        H[i, i - 2] = -1 / 12
        H[i, i - 1] = 4 / 3
        H[i, i]     = -2.5
        H[i, i + 1] = 4 / 3
        H[i, i + 2] = -1 / 12
    H = -H / dx**2
    for i in range(Nx):
        H[i, i] += V[i]

    return H.tocsc(), x, V

def inverse_power_iteration(H, tol, shift, max_iter):
    N = H.shape[0]
    I = sparse.identity(N, format='csc')
    H_shifted = H - shift * I
    lu = sla.splu(H_shifted)

    v = np.random.rand(N)
    v /= np.linalg.norm(v)

    for _ in range(max_iter):
        w = lu.solve(v)
        w /= np.linalg.norm(w)
        if np.linalg.norm(w - v) < tol:
            break
        v = w

    eigval = w @ (H @ w)
    return eigval, w

def find_multiple_eigenpairs_inverse(H, tol, shifts, max_iter):
    eigenvalues = []
    eigenvectors = []

    for shift in shifts:
        eigval, eigvec = inverse_power_iteration(H, tol , shift, max_iter=max_iter)
        if any(np.isclose(eigval, ev, tol) for ev in eigenvalues):
            continue
        eigenvalues.append(eigval)
        eigenvectors.append(eigvec)

    return np.array(eigenvalues), np.column_stack(eigenvectors)

def analytical_eigenvalues(n, omega=1):
    return np.array([hbar * omega * 2 * (i + 0.5) for i in range(n)])

def check_orthonormality(psi_matrix, method="FD", tol=1e-6):
    norms = np.linalg.norm(psi_matrix, axis=0)
    all_normalized = np.all(np.abs(norms - 1) < tol)

    overlap = psi_matrix.T @ psi_matrix
    identity = np.eye(psi_matrix.shape[1])
    is_orthogonal = np.all(np.abs(overlap - identity) < tol)
    return all_normalized and is_orthogonal


def plot_last_eigenvectors(x, evt_fd, evt_ipi, Neigs, count=3):
    plt.figure(figsize=(8, 6))
    for i in range(Neigs - count, Neigs):
        plt.plot(x, evt_fd[:, i], color='blue', linestyle='--', linewidth=1)
        plt.plot(x, evt_ipi[:, i], color='red', linestyle='-', linewidth=1)

    custom_lines = [
        Line2D([0], [0], color='blue', linestyle='--', lw=2),
        Line2D([0], [0], color='red', linestyle='-', lw=2)
    ]
    plt.legend(custom_lines, ['FD (dashed)', 'IPI (solid)'], fontsize='small')
    plt.xlabel("x")
    plt.ylabel(r"$\psi_n(x)$")
    plt.title(f"Last {count} eigenvectors of the Harmonic Oscillator")
    plt.grid()
    plt.tight_layout()
    plt.show()

def run(Neigs, Nround, pert=False, shifts=None, title="Harmonic Oscillator"):
    H, x, V = build_hamiltonian(xmin, xmax, Nx, A, alpha, pert=pert)

    # FD
    t0 = time.time()
    evl_fd, evt_fd = sla.eigs(H, k=Neigs, which='SM')
    t_fd = time.time() - t0
    evl_fd = np.real(evl_fd)
    for i in range(Neigs):
        evt_fd[:, i] /= np.linalg.norm(evt_fd[:, i])

    # IPI
    if shifts is None:
        shifts = [2 * (i + 0.5) for i in range(Neigs)]
    t0 = time.time()
    evl_ipi, evt_ipi = find_multiple_eigenpairs_inverse(H, tol, shifts, max_iter)
    t_ipi = time.time() - t0
    idx = np.argsort(evl_ipi)
    evl_ipi = evl_ipi[idx]
    evt_ipi = evt_ipi[:, idx]

    print(f"\nFD method time       : {t_fd:.4f} s")
    print(f"IPI method time      : {t_ipi:.4f} s")

    print("\nLast few FD eigenvalues:")
    print(np.round(np.sort(evl_fd[-5:-1]), Nround))

    print("\nLast few IPI eigenvalues:")
    print(np.round(evl_ipi[-5:-1], Nround))

    if not pert:
        evl_analytical = analytical_eigenvalues(Neigs)
        print("\nLast few analytical eigenvalues:")
        print(np.round(evl_analytical[-5:-1], Nround))

    fd_ok = check_orthonormality(evt_fd, "FD")
    ipi_ok = check_orthonormality(evt_ipi, "IPI")
    
    print(f"\nFD eigenvectors orthonormal: {fd_ok}")
    print(f"IPI eigenvectors orthonormal: {ipi_ok}")
    
    plot_last_eigenvectors(x, evt_fd, evt_ipi, Neigs)

if __name__ == "__main__":
    run(Neigs, Nround, pert=False, title="Unperturbed Harmonic Oscillator")
    # run(Neigs, Nround, pert=True, title="Perturbed Harmonic Oscillator")
