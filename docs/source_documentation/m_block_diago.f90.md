# `m_block_diago.f90`

## Overview

The `m_block_diago` Fortran module in MOLGW provides routines for solving the specific block matrix eigenvalue problem that arises in linear response theories such as Time-Dependent Density Functional Theory (TDDFT), Random Phase Approximation (RPA), and the Bethe-Salpeter Equation (BSE). This problem is typically formulated as:

```
(  A  B  ) ( X )    ( X )
( -B -A  ) ( Y )  = ( Y ) * Omega
```
or a similar variant, where `A` and `B` are matrices derived from orbital energy differences and electron repulsion integrals, `X` and `Y` are the solution vectors, and `Omega` are the excitation energies. This module offers different algorithms to perform this diagonalization, including methods based on Cholesky decomposition and the iterative Davidson algorithm, with support for ScaLAPACK for parallel execution.

## Key Components

*   **`diago_4blocks_chol(amb_matrix, apb_matrix, desc_apb, bigomega, xpy_matrix, xmy_matrix, desc_x)`**:
    *   **Purpose**: Solves the block matrix eigenvalue problem using a Cholesky decomposition-based approach. This method typically transforms the problem into a standard symmetric eigenvalue problem of the form `(A-B)^(1/2) (A+B) (A-B)^(1/2) Z = Z Omega^2` or `(A+B)^(-1/2) (A-B) (A+B)^(-1/2) Z = Z Omega^(-2)`.
    *   If ScaLAPACK is available, it uses the `PDBSSOLVER1` routine. Otherwise, it employs standard LAPACK routines (`DPOTRF`, `DSYGST`, `DSYEV`, `DTRMM`, `DTRSM`).
    *   **Output**: `bigomega` (excitation energies), `xpy_matrix` (X+Y vectors), `xmy_matrix` (X-Y vectors).

*   **`diago_4blocks_rpa_sca(amb_diag_rpa, apb_matrix, desc_apb, bigomega, xpy_matrix, desc_x)`**:
    *   **Purpose**: A specialized ScaLAPACK-based routine designed for cases where the (A-B) matrix is diagonal. This is common in standard RPA (direct RPA) calculations where (A-B) consists only of orbital energy differences.
    *   It solves the symmetric eigenproblem `(A-B)^(1/2) * (A+B) * (A-B)^(1/2) * Z' = Z' * Omega^2`.
    *   Can optionally use the ELPA library for diagonalization if available.
    *   **Output**: `bigomega` (excitation energies), `xpy_matrix` (transformed eigenvectors, from which X+Y can be reconstructed).

*   **`diago_4blocks_davidson(toldav, nstep, amb_diag_rpa, amb_matrix, apb_matrix, desc_apb, bigomega, xpy_matrix, xmy_matrix, desc_x)`**:
    *   **Purpose**: Implements the Davidson iterative algorithm to find a subset of the lowest excitation energies and corresponding eigenvectors. This method is particularly useful for large systems where forming or diagonalizing the full A and B matrices is computationally prohibitive.
    *   It iteratively builds a subspace using preconditioned residual vectors and solves the eigenvalue problem within this smaller subspace.
    *   **Input**: `toldav` (convergence tolerance for residuals), `nstep` (maximum Davidson iterations), `amb_diag_rpa` (diagonal preconditioner, typically orbital energy differences).
    *   **Output**: `bigomega` (excitation energies), `xpy_matrix` (X+Y vectors), `xmy_matrix` (X-Y vectors) for the converged eigenvalues.

## Important Variables/Constants

*   **`amb_matrix` (Real(dp), Intent(inout))**: Represents the (A-B) block of the Casida matrix.
*   **`apb_matrix` (Real(dp), Intent(inout))**: Represents the (A+B) block of the Casida matrix.
*   **`desc_apb`, `desc_x` (Integer, Intent(in))**: ScaLAPACK array descriptors for `apb_matrix` (and `amb_matrix`) and the output eigenvector matrices respectively.
*   **`bigomega(:)` (Real(dp), Intent(out))**: Array storing the computed excitation energies (eigenvalues &Omega;).
*   **`xpy_matrix(:, :)` (Real(dp), Intent(out))**: Array storing the (X+Y) solution vectors.
*   **`xmy_matrix(:, :)` (Real(dp), Intent(out))**: Array storing the (X-Y) solution vectors.
*   **`amb_diag_rpa(:)` (Real(dp), Intent(in))**: Used in `diago_4blocks_rpa_sca` and as a preconditioner in `diago_4blocks_davidson`, typically containing orbital energy differences (&epsilon;<sub>a</sub> - &epsilon;<sub>i</sub>).
*   **`postscf_diago_flavor` (Character string)**: Input parameter (from `m_inputparam`) that can select specific LAPACK/ScaLAPACK diagonalization routines (e.g., 'STD' for standard, 'DC' for divide-and-conquer).

## Usage Examples

These subroutines are typically called by other modules within MOLGW, such as `m_linear_response` or `m_spectra`, after the (A-B) and (A+B) matrices (or their action on vectors) have been constructed.

```fortran
! Conceptual call within a linear response routine
USE m_block_diago
! ... A_minus_B_matrix and A_plus_B_matrix are prepared ...
! ... ScaLAPACK descriptors desc_AB and desc_eigenvecs are set up ...
REAL(DP), ALLOCATABLE :: eigenvalues(:), eigenvectors_XplusY(:,:), eigenvectors_XminusY(:,:)
INTEGER :: n_excitations_to_find, matrix_size

n_excitations_to_find = ...
matrix_size = ... ! Size of the A/B matrices
ALLOCATE(eigenvalues(n_excitations_to_find))
ALLOCATE(eigenvectors_XplusY(matrix_size, n_excitations_to_find)) ! Or ScaLAPACK local sizes
ALLOCATE(eigenvectors_XminusY(matrix_size, n_excitations_to_find)) ! Or ScaLAPACK local sizes

IF (use_davidson_solver) THEN
  CALL diago_4blocks_davidson(davidson_tolerance, max_davidson_iter, &
                              A_minus_B_diagonal_elements, A_minus_B_matrix, A_plus_B_matrix, desc_AB, &
                              eigenvalues, eigenvectors_XplusY, eigenvectors_XminusY, desc_eigenvecs)
ELSE
  CALL diago_4blocks_chol(A_minus_B_matrix, A_plus_B_matrix, desc_AB, &
                          eigenvalues, eigenvectors_XplusY, eigenvectors_XminusY, desc_eigenvecs)
END IF
```

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_timing`, `m_warning`, `m_memory`, `m_mpi`, `m_scalapack`.
*   **Linear Algebra**: `m_linear_algebra` (for general matrix operations or if ScaLAPACK is not used).
*   **Input Parameters**: `m_inputparam` (for `postscf_diago_flavor`).
*   **External Libraries**:
    *   **LAPACK/BLAS**: Required for serial matrix operations (e.g., `DPOTRF`, `DSYGST`, `DSYEV`, `DTRMM`, `DTRSM`, `DSYMM`, `DGEMM`).
    *   **ScaLAPACK**: (Optional, if `HAVE_SCALAPACK` is defined) Used for parallel, distributed-memory diagonalization and matrix operations (e.g., `PDBSSOLVER1`, `PDLACPY`, `PDSYMM`, `PDGEMM`).
    *   **ELPA (Eigenvalue Solvers for Petaflop Applications)**: (Optional, if `HAVE_ELPA` is defined) Can be used as an alternative, high-performance eigensolver within `diago_4blocks_rpa_sca`.
*   **Calling Modules**: This module is primarily called by modules that set up and require the solution of the Casida equations or similar linear response eigenvalue problems, such as `m_linear_response` or `m_spectra`.
```
