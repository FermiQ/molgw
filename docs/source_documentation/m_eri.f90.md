# `m_eri.f90`

## Overview

The `m_eri` Fortran module in MOLGW is central to handling two-electron repulsion integrals (ERIs). It manages the storage, indexing, and screening of 2-center, 3-center, and 4-center Coulomb integrals. Key functionalities include identifying and neglecting numerically small shell pairs via Cauchy-Schwarz screening to reduce computational cost, setting up efficient indexing schemes for accessing stored integrals, and managing the distribution of 3-center ERIs in parallel computing environments (MPI/ScaLAPACK). It also provides routines for saving and loading 4-center ERIs to/from disk.

## Key Components

*   **Module Variables for Storing ERIs**:
    *   `eri_4center(:)` (Real(dp), Allocatable, Public): A 1D array storing all significant 4-center, 2-electron integrals (ij|kl). Accessed via `index_eri` function.
    *   `eri_4center_lr(:)` (Real(dp), Allocatable, Public): Similar to `eri_4center`, but for the long-range component of ERIs if range-separated Coulomb interaction is used.
    *   `eri_3center(:, :)` (Real(dp), Allocatable, Public): A 2D array storing 3-center ERIs, typically of the form (AuxiliaryBasisFunction | BasisFunctionPair). Dimensions are (`npair`, `nauxil_local`), where `nauxil_local` refers to the number of auxiliary basis functions on the current MPI process.
    *   `eri_3center_lr(:, :)` (Real(dp), Allocatable, Public): Similar to `eri_3center`, but for the long-range component.

*   **`prepare_eri(basis)`**:
    *   **Purpose**: Initializes the ERI handling system.
    *   Sets the integral screening tolerance (`TOL_INT`) based on `integral_level` from input.
    *   Calls `identify_negligible_shellpair` to screen out shell pairs that would lead to negligible ERIs.
    *   Calls `setup_shell_index`, `setup_shellpair`, and `setup_basispair` to create efficient indexing schemes for non-negligible integrals.
    *   Calculates `nint_4center`, the total number of 4-center ERIs to be stored (using 8-byte integers for potentially large counts).

*   **Screening and Indexing Subroutines**:
    *   `identify_negligible_shellpair(basis)`: Uses Cauchy-Schwarz inequality on diagonal shell-pair integrals (e.g., (ij|ij)) to flag shell pairs whose contributions to 4-center ERIs will likely be below `TOL_INT`.
    *   `setup_shell_index(basis)`: Associates each basis function with its shell index.
    *   `setup_shellpair(basis)`: Creates a list of non-negligible shell pairs, ordered to optimize calls to integral libraries.
    *   `setup_basispair(basis)`: Creates a unique index for each non-negligible pair of basis functions (ij) and stores the actual (i,j) indices. Defines `npair`.
    *   `index_eri(ibf, jbf, kbf, lbf)` (Function): Returns a 1D packed index for the 4-center integral (ibf jbf|kbf lbf).
    *   `index_pair(ibf, jbf)` (Function): Returns a 1D packed index for a pair of basis functions (ibf, jbf).
    *   `negligible_basispair(ibf, jbf)` (Function): Returns `.TRUE.` if the shell pair corresponding to basis functions `ibf` and `jbf` was marked as negligible.

*   **ERI Access Functions**:
    *   `eri(ibf, jbf, kbf, lbf)` (Elemental Function): Provides access to the stored 4-center integral (ibf jbf|kbf lbf). Returns 0.0 if either pair (ibf,jbf) or (kbf,lbf) is negligible.
    *   `eri_lr(ibf, jbf, kbf, lbf)` (Function): Similar to `eri` but accesses `eri_4center_lr`.

*   **Parallel Distribution for 3-center ERIs**:
    *   `distribute_auxil_basis(nbf_auxil_basis)`: Distributes auxiliary basis functions across MPI processes using a ScaLAPACK-like block-cyclic distribution. Sets up mappings between global and local auxiliary basis indices.
    *   `distribute_auxil_basis_lr(...)`: Same for long-range auxiliary basis.
    *   `reshuffle_distribution_3center()`: If necessary, redistributes the `eri_3center` array from an initial distribution (e.g., suitable for calculation) to a final distribution (e.g., suitable for specific ScaLAPACK operations).

*   **I/O Routines**:
    *   `dump_out_eri(rcut)`: Writes the `eri_4center` array to a binary file (`molgw_eri.data` or `molgw_eri_lr.data`).
    *   `read_eri(rcut)` (Logical Function): Reads `eri_4center` from a binary file. Returns `.TRUE.` if successful.

*   **Deallocation Routines**:
    *   `deallocate_eri_4center()`, `deallocate_eri_4center_lr()`, `deallocate_eri()`, `destroy_eri_3center_lowerlevel()`, `destroy_eri_3center_lr()`: Manage deallocation of the large ERI arrays and associated indexing arrays.

## Important Variables/Constants

*   **`TOO_LOW_EIGENVAL` (Real(dp), Parameter)**: A threshold for small eigenvalues, likely used in the context of forming the (V<sup>-1/2</sup>) matrix for RI, value 1.0e-6.
*   **`TOL_INT` (Real(dp), Private)**: Tolerance used for screening ERIs based on Cauchy-Schwarz estimates. Its value is set based on the `integral_level` input parameter.
*   **`nint_4center` (Integer(kind=int8), Protected)**: The total number of 4-center ERIs to be stored. Uses an 8-byte integer to accommodate large numbers.
*   **`npair` (Integer, Protected)**: The number of unique, non-negligible pairs of basis functions (ij).
*   **`nauxil_global`, `nauxil_local` (Integer, Public/Protected)**: Global and local (per MPI process) counts of auxiliary basis functions, crucial for managing distributed 3-center ERI arrays.
*   **`desc_eri3`, `desc_eri3_lr` (Integer, Public)**: ScaLAPACK array descriptors for the distributed 3-center ERI arrays.

## Usage Examples

This module is primarily used internally by MOLGW.
1.  **Initialization**: `prepare_eri` is called early in a calculation after the basis set is defined.
2.  **ERI Calculation**: Other modules (e.g., `m_eri_calculate`, or interfaces to Libcint/Libint) compute the ERIs and populate the arrays managed by `m_eri` (e.g., `eri_4center`, `eri_3center`).
3.  **ERI Access**: Modules performing post-HF calculations (like `m_mp2`, `m_gw_selfenergy`, `m_ci`) then use the `eri(i,j,k,l)` function or directly access `eri_3center` (with appropriate parallel considerations) to get the required integrals.

```fortran
! Conceptual: How another module might use m_eri after ERIs are calculated

USE m_eri
USE m_definitions
IMPLICIT NONE

INTEGER :: i, j, k, l ! Basis function indices
REAL(DP) :: J_term, K_term

! ... (Assume i, j, k, l are valid indices for occupied/virtual orbitals) ...

! Get a Coulomb integral (ii|jj)
J_term = eri(i, i, j, j)

! Get an exchange integral (ij|ji)
K_term = eri(i, j, j, i)

! ... use J_term and K_term in a calculation ...
```

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_mpi`, `m_memory`, `m_warning`, `m_scalapack`, `m_timing`, `m_inputparam`.
*   **Basis Set**: `m_basis_set` (provides `basis_set` type, shell and basis function details).
*   **Transformations**: `m_cart_to_pure` (used by `transform_libint_to_molgw` to convert integrals if basis is 'PURE').
*   **Integral Libraries**: `m_libint_tools` and `m_libcint_tools` are used in `identify_negligible_shellpair` for computing the diagonal (ii|jj) type integrals needed for Cauchy-Schwarz screening.
*   **ERI Calculation Modules (e.g., `m_eri_calculate`)**: These modules are responsible for actually computing the ERI values and storing them into the arrays (`eri_4center`, `eri_3center`) defined and managed by `m_eri`.
*   **Post-HF Method Modules**: Modules implementing methods like MP2, GW, BSE, CI, etc., will heavily rely on this module to access the computed two-electron integrals.
*   The efficiency of ERI handling, particularly the screening and indexing, is critical for the performance of the entire MOLGW package for correlated calculations.
```
