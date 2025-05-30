# `m_eri_ao_mo.f90`

## Overview

The `m_eri_ao_mo` Fortran module in MOLGW is responsible for transforming two-electron repulsion integrals (ERIs) from the Atomic Orbital (AO) basis to the Molecular Orbital (MO) basis. This transformation is a critical step for most post-Hartree-Fock quantum chemistry calculations, as these methods are typically formulated in the MO basis. The module provides routines to perform these transformations for 3-center (Resolution of Identity - RI) and 4-center ERIs, handling both real and complex-valued molecular orbital coefficients (e.g., for exact two-component relativistic calculations, x2c).

## Key Components

*   **Module Variables for Storing MO ERIs**:
    *   `eri_3center_eigen(:,:,:,:)` (Real(dp), Protected, Allocatable): Stores 3-center ERIs transformed to the MO basis, having the form (Q|&phi;<sub>i</sub>&phi;<sub>j</sub>), where Q is an auxiliary basis function and &phi;<sub>i</sub>, &phi;<sub>j</sub> are molecular orbitals. Dimensions: (nauxil_local, n_mo_i, n_mo_j, nspin).
    *   `eri_3center_eigen_lr(:,:,:,:)` (Real(dp), Protected, Allocatable): Similar to `eri_3center_eigen`, but for the long-range part of the ERIs if range separation is used.
    *   `eri_4center_eigen_uks(:,:,:,:)` (Real(dp), Protected, Allocatable): Stores 4-center ERIs fully transformed to the MO basis: (&phi;<sub>i</sub>&phi;<sub>j</sub>|&phi;<sub>k</sub>&phi;<sub>l</sub>). Primarily for unrestricted Kohn-Sham (UKS) or spin-restricted calculations where `nspin=1` effectively.
    *   `eri_3center_eigen_cmplx(:,:,:,:)` (Complex(dp), Protected, Allocatable): Complex version of `eri_3center_eigen` for use with complex MO coefficients.
    *   `eri_3center_eigen_lr_cmplx(:,:,:,:)` (Complex(dp), Protected, Allocatable): Complex version for long-range 3-center ERIs.

*   **ERI Access Functions (MO basis)**:
    *   **`eri_eigen(istate, jstate, ijspin, kstate, lstate, klspin)` (Function)**: Returns the 4-index MO integral (&phi;<sub>istate</sub>&phi;<sub>jstate</sub>|&phi;<sub>kstate</sub>&phi;<sub>lstate</sub>).
        *   If `has_auxil_basis` is true, it reconstructs the 4-center MO integral from 3-center MO integrals: &sum;<sub>Q</sub> (istate jstate|Q)(Q|kstate lstate).
        *   Otherwise, it directly returns a pre-transformed 4-center integral from `eri_4center_eigen_uks`.
    *   **`eri_eigen_ri(...)`, `eri_eigen_ri_lr(...)`, `eri_eigen_ri_cmplx(...)`, `eri_eigen_ri_lr_cmplx(...)` (Functions)**: These are specialized versions that *always* reconstruct the 4-index MO integrals using the RI approximation from the respective 3-center MO integral arrays (real, long-range real, complex, long-range complex).
    *   **`eri_eigen_ri_x2c(...)` (Function)**: Computes 4-index MO integrals for x2c calculations using complex 3-center MO integrals, correctly summing over spin components.
    *   **`eri_eigen_ri_paral(...)` (Pure Function)**: A `pure` version of `eri_eigen_ri`, intended for use in contexts requiring function purity. It performs the dot product over auxiliary functions.

*   **AO to MO Transformation Subroutines**:
    *   **`calculate_eri_4center_eigen(c_matrix, istate, ijspin, eri_eigenstate_i)`**: Performs a partial AO-to-MO transformation for 4-center integrals, transforming two AO indices of (AO<sub>p</sub> AO<sub>q</sub> | AO<sub>r</sub> AO<sub>s</sub>) to get (MO<sub>istate</sub> AO<sub>q</sub> | AO<sub>r</sub> AO<sub>s</sub>) and subsequently further transforming to get (MO<sub>istate</sub> MO<sub>jstate</sub> | MO<sub>kstate</sub> MO<sub>lstate</sub>). It is optimized using the 8-fold permutation symmetry of ERIs.
    *   **`calculate_eri_4center_eigen_uks(c_matrix, nstate_min, nstate_max)`**: Transforms all 4-center AO ERIs to the MO basis for a specified range of MOs, specifically for spin-restricted cases (nspin=1). Stores results in `eri_4center_eigen_uks`.
    *   **`calculate_eri_3center_eigen(...)`**: Transforms 3-center AO ERIs (Q|&mu;&nu;) to the MO basis (Q|&phi;<sub>i</sub>&phi;<sub>j</sub>). This is effectively a two-index transformation on the &mu;&nu; pair for each auxiliary function Q. Handles real coefficients.
    *   **`calculate_eri_3center_eigen_lr(...)`**: Same as above, but for long-range 3-center ERIs.
    *   **`calculate_eri_3center_eigen_cmplx(...)`**: Same as `calculate_eri_3center_eigen` but for complex MO coefficients.
    *   **`calculate_eri_x2c(...)`**: A driver for x2c calculations. It adapts the complex relativistic MO coefficients (which have alpha and beta spin components) and calls `calculate_eri_3center_eigen_cmplx` to perform the transformation.

*   **`form_erimol(...)`**:
    *   Performs a full four-index transformation of AO ERIs to MO ERIs: (&phi;<sub>p</sub>&phi;<sub>i</sub>|&phi;<sub>j</sub>&phi;<sub>k</sub>) = &sum;<sub>&mu;&nu;&lambda;&sigma;</sub> C<sub>&mu;p</sub>C<sub>&nu;i</sub>C<sub>&lambda;j</sub>C<sub>&sigma;k</sub> (&mu;&nu;|&lambda;&sigma;).
    *   Handles both real (`c_matrix`, `ERImol`) and complex (`c_matrix_cmplx`, `ERImol_cmplx`) MO coefficients.

*   **Deallocation Subroutines**:
    *   `destroy_eri_4center_eigen_uks()`, `destroy_eri_3center_eigen()`, `destroy_eri_3center_eigen_cmplx()`, `destroy_eri_3center_eigen_x2c()`: Deallocate the respective MO ERI arrays.

## Important Variables/Constants

*   **`has_auxil_basis` (Logical from `m_inputparam`)**: A crucial flag that determines the strategy for obtaining 4-index MO integrals. If true, the RI approximation is used with 3-center integrals.
*   **`eri_4center(:)` (Real(dp) from `m_eri`)**: Array of 4-center ERIs in the AO basis.
*   **`eri_3center(:, :)` (Real(dp) from `m_eri`)**: Array of 3-center ERIs in the AO basis (Aux|AO AO).
*   **`c_matrix(:, :, :)` (Real(dp) or Complex(dp))**: Input array of MO coefficients in the AO basis.

## Usage Examples

This module is primarily used by other modules within MOLGW that perform correlated calculations in the MO basis.
1.  After AO integrals are computed (by `m_eri_calculate` or Libcint/Libint wrappers and stored in `m_eri`), and MO coefficients `c_matrix` are available (e.g., from an SCF calculation via `m_scf`).
2.  One of the `calculate_eri_..._eigen` subroutines is called to transform the necessary ERIs into the MO basis and store them in the module's public arrays (e.g., `eri_3center_eigen` or `eri_4center_eigen_uks`).
3.  Subsequent routines (e.g., in `m_mp2`, `m_gw_selfenergy`, `m_ci`, `m_build_bse`) then use functions like `eri_eigen(...)` or `eri_eigen_ri(...)` to retrieve these MO integrals for their computations.

```fortran
! Conceptual: Within an MP2 energy calculation module
USE m_eri_ao_mo
USE m_definitions
IMPLICIT NONE

REAL(DP) :: c_matrix_scf(n_ao, n_mo, n_spin_scf)
! ... (c_matrix_scf is obtained) ...

! Transform 3-center ERIs to MO basis (if RI-MP2)
IF (has_auxil_basis) THEN
  CALL calculate_eri_3center_eigen(c_matrix_scf)
END IF

! ... later, to get an integral (ia|jb) needed for MP2 ...
REAL(DP) :: mo_integral_iajb
INTEGER :: i_occ, a_virt, j_occ, b_virt, spin_channel

mo_integral_iajb = eri_eigen(i_occ, a_virt, spin_channel, j_occ, b_virt, spin_channel) 
! This will use eri_3center_eigen if has_auxil_basis, or eri_4center_eigen_uks otherwise.
```

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_mpi`, `m_scalapack`, `m_memory`, `m_warning`, `m_timing`, `m_inputparam`.
*   **Basis Set**: `m_basis_set` (not directly used in many functions here, but `c_matrix` is in this basis).
*   **AO Integrals**: `m_eri` is a primary dependency, providing the raw AO-basis ERIs (`eri_4center`, `eri_3center`) and their indexing/screening information.
*   **Linear Algebra**: Uses BLAS routines like `DGEMM` (implicitly via `MATMUL` or directly) and `DOT_PRODUCT` for the transformation steps.
*   **Parallelism**: OpenMP pragmas (`!$OMP`) are used for loop parallelization in `calculate_eri_4center_eigen`. MPI communication (`auxil%sum`) is used in the `eri_eigen_ri...` functions to sum contributions when 3-center integrals are distributed.
*   This module acts as a crucial bridge: it takes AO integrals (from `m_eri` and underlying integral engines) and MO coefficients (from `m_scf` or other sources) and produces MO integrals needed by various electron correlation methods.
```
