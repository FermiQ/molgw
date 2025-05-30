# `m_linear_algebra.f90`

## Overview

The `m_linear_algebra` Fortran module in MOLGW provides a collection of essential linear algebra routines. These routines are primarily wrappers around standard Basic Linear Algebra Subprograms (BLAS) and Linear Algebra Package (LAPACK) routines, offering a consistent interface within the MOLGW environment. The module includes functionalities such as matrix diagonalization (eigenvalue and eigenvector computation), matrix inversion, matrix trace, symmetry checks, vector cross product, and more specialized algorithms like Davidson diagonalization and joint diagonalization. It supports double precision real (`dp`), single precision real (`sp`), and double precision complex (`cdp`) data types for many operations through generic interfaces.

## Key Components

*   **Generic Interfaces**:
    *   **`invert`**: Wraps `invert_dp`, `invert_inplace_dp` (for real matrices), `invert_cdp`, `invert_inplace_cdp` (for complex matrices).
    *   **`invert_symmetric`**: Wraps `invert_symmetric_inplace_dp` for real symmetric matrices.
    *   **`diagonalize_wo_vectors`**: Wraps `diagonalize_wo_vectors_dp` (computes eigenvalues only for real symmetric matrices).
    *   **`diagonalize`**: Wraps routines for real symmetric (`_dp`, `_sp`) and complex Hermitian (`_cdp`) matrices, including `_inplace` versions.
    *   **`matrix_lower_to_full`**: Wraps `matrix_lower_to_full_dp` (symmetric) and `matrix_lower_to_full_cdp` (Hermitian).

*   **Matrix Manipulation and Query**:
    *   **`matrix_lower_to_full_dp(matrix)`**: Fills the upper triangle of a real matrix `matrix` assuming it is symmetric (upper(i,j) = lower(j,i)).
    *   **`matrix_lower_to_full_cdp(matrix)`**: Fills the upper triangle of a complex matrix `matrix` assuming it is Hermitian (upper(i,j) = CONJG(lower(j,i))).
    *   **`matrix_trace(matrix)` (Function)**: Returns the trace of a real square matrix.
    *   **`matrix_trace_cmplx(matrix)` (Function)**: Returns the trace of a complex square matrix.
    *   **`matrix_is_symmetric(matrix)` (Function)**: Checks if a real square matrix is symmetric within a tolerance (1.0e-5).

*   **Matrix Inversion**:
    *   **`invert_dp(matrix, matrix_inv)` / `invert_cdp(matrix, matrix_inv)`**: Computes the inverse of `matrix` and stores it in `matrix_inv`.
    *   **`invert_inplace_dp(matrix)` / `invert_inplace_cdp(matrix)`**: Computes the inverse of a general real/complex matrix `matrix` in place using LU factorization (LAPACK's `DGETRF`/`ZGETRF` followed by `DGETRI`/`ZGETRI`).
    *   **`invert_symmetric_inplace_dp(matrix)`**: Computes the inverse of a real symmetric matrix `matrix` in place using Bunch-Kaufman factorization (LAPACK's `DSYTRF` followed by `DSYTRI`).

*   **Matrix Diagonalization (Eigen Solvers)**:
    *   **`diagonalize_wo_vectors_dp(flavor, matrix, eigval)`**: Computes only the eigenvalues (`eigval`) of a real symmetric matrix. `matrix` is overwritten.
    *   **`diagonalize_dp(flavor, matrix, eigval, eigvec)` / `diagonalize_sp(...)` / `diagonalize_cdp(...)`**: Computes all eigenvalues (`eigval`) and eigenvectors (`eigvec`) of a matrix (real symmetric for `_dp`/`_sp`, complex Hermitian for `_cdp`). The input `matrix` is copied to `eigvec` first.
        *   `flavor` (Character): Selects the LAPACK driver:
            *   'R': `DSYEVR` / `ZHEEVR` (eigenvalues/vectors using relatively robust representations).
            *   'D': `DSYEVD` / `ZHEEVD` (divide and conquer).
            *   'S': `SSYEV` (single precision standard).
            *   'E': `SSYEVD` (single precision divide and conquer).
            *   Default: `DSYEV` / `ZHEEV` (standard).
    *   **`diagonalize_inplace_dp(...)` / `diagonalize_inplace_sp(...)` / `diagonalize_inplace_cdp(...)`**: Similar to `diagonalize`, but the input `matrix` is directly overwritten with the eigenvectors.

*   **Iterative Eigensolver**:
    *   **`diagonalize_davidson(tolerance, nstep, ham, neig, eigval, eigvec)`**: Implements the Davidson iterative algorithm to find a specified number (`neig`) of eigenvalues and eigenvectors of a large real symmetric Hamiltonian matrix `ham`. Useful when only a few eigensolutions are needed.

*   **Vector Operations**:
    *   **`orthogonalize(vec)`**: Orthonormalizes a set of column vectors `vec(:, :)` using the Gram-Schmidt procedure.
    *   **`check_unitarity(cmat)`**: Verifies if a given complex matrix `cmat` is unitary (C<sup>&dagger;</sup>C = I and CC<sup>&dagger;</sup> = I) within a numerical tolerance.
    *   **`cross_product(u1, u2, u3)`**: Computes the vector cross product **u3** = **u1** &times; **u2**.
    *   **`determinant_3x3_matrix(mat)` (Function)**: Calculates the determinant of a 3x3 real matrix.

*   **`joint_diagonalization(A, tol, V, converged)`**:
    *   Performs a Jacobi-like joint diagonalization of a series of `n` real symmetric matrices A(:,:,i) of size `m`x`m`. It finds a single orthogonal transformation matrix `V` such that all V<sup>T</sup>A(:,:,i)V are made as diagonal as possible. Iterates until convergence based on `tol`.

## Important Variables/Constants

*   **`flavor` (Character, Input for diagonalization routines)**: Allows selection of different LAPACK eigensolver algorithms (e.g., 'R', 'D', 'S', 'E', or default). This choice can impact performance and numerical accuracy depending on the matrix properties.
*   Numerical tolerances are hardcoded in routines like `matrix_is_symmetric` (1.0e-5) and `check_unitarity` (1.0e-9).

## Usage Examples

This module provides general-purpose linear algebra functionalities that are widely used across MOLGW.

Diagonalizing a Fock matrix in an SCF procedure:
```fortran
USE m_linear_algebra
USE m_definitions
IMPLICIT NONE

REAL(DP) :: fock_matrix(100, 100), mo_coeffs(100, 100), mo_energies(100)
CHARACTER(LEN=1) :: diag_flavor = 'D' ! Use DSYEVD (divide and conquer)

! ... Fock_matrix is constructed ...

CALL diagonalize_inplace_dp(diag_flavor, fock_matrix, mo_energies)
! fock_matrix now contains the eigenvectors (mo_coeffs)
mo_coeffs = fock_matrix 
! mo_energies contains the eigenvalues
```

Inverting an overlap matrix (if needed, though S<sup>-1/2</sup> is more common):
```fortran
USE m_linear_algebra
USE m_definitions
IMPLICIT NONE

REAL(DP) :: S_matrix(50,50), S_inverse(50,50)
! ... S_matrix is constructed ...

CALL invert_dp(S_matrix, S_inverse)
! S_inverse now holds the inverse of S_matrix
```

## Dependencies and Interactions

*   **`m_definitions`**: For precision kinds `dp` and `sp`.
*   **`m_warning`**: For error handling (`die`, `issue_warning`).
*   **LAPACK and BLAS Libraries**: This module is fundamentally a set of wrappers around LAPACK routines (e.g., `DGETRF`, `DGETRI`, `DSYTRF`, `DSYTRI`, `DSYEV`, `DSYEVD`, `DSYEVR`, `ZHEEV`, `ZHEEVD`, `ZHEEVR`, `SSYEV`, `SSYEVD`) and BLAS routines (e.g., `DOT_PRODUCT`, `MATMUL` implicitly, `DSYMV`, `DGEMM`). The correct linking to these numerical libraries is essential for MOLGW to function.
*   **Other MOLGW Modules**:
    *   SCF routines (`m_scf`) heavily rely on `diagonalize` or `diagonalize_hamiltonian_scalapack` (from `m_hamiltonian_tools`, which might use routines from here or ScaLAPACK directly).
    *   `m_hamiltonian_tools` uses `setup_x_matrix` which involves diagonalization.
    *   Any module performing transformations or requiring matrix factorizations might use routines from `m_linear_algebra`.
```
