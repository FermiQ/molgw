# `m_eri_calculate.f90`

## Overview

The `m_eri_calculate` Fortran module in MOLGW is responsible for orchestrating the computation of two-electron repulsion integrals (ERIs). It handles 2-center, 3-center, and 4-center ERIs, which are fundamental for many quantum chemistry methods. The module interfaces with underlying integral libraries (Libcint or a native Libint implementation wrapped in `m_libint_tools`) to perform the actual mathematical calculations over primitive Gaussian basis functions. It works in conjunction with the `m_eri` module, which manages the storage, indexing, and screening of these integrals. This module includes logic for Resolution of Identity (RI) approximations, which involve 2-center and 3-center integrals with an auxiliary basis set.

## Key Components

*   **`calculate_eri(print_eri_, basis, rcut)`**:
    *   **Purpose**: Main driver subroutine for calculating 4-center ERIs.
    *   If `incore_` (input parameter) is true, it allocates memory for `eri_4center` (or `eri_4center_lr` if `rcut` is significant for range separation).
    *   It first attempts to read the integrals from a pre-existing file (`molgw_eri.data` or `molgw_eri_lr.data`) using the `read_eri` function (from `m_eri`).
    *   If reading fails or is not requested, it calls `calculate_eri_4center` to compute them.
    *   If `print_eri_` is true, it calls `dump_out_eri` (from `m_eri`) to save the computed integrals.

*   **`calculate_eri_4center(basis, rcut)`**:
    *   **Purpose**: Computes all significant 4-center ERIs (ij|kl) by iterating over pairs of non-negligible shell pairs (determined by `m_eri`).
    *   For each quartet of shells, it calls either Libcint (`cint2e_cart`) or a native Libint wrapper (`libint_4center`) to get integrals over primitive Cartesian Gaussians.
    *   These are then transformed to the final basis representation (Cartesian or Pure Spherical Harmonics) using `transform_libint_to_molgw`.
    *   Stores the results in the `eri_4center` or `eri_4center_lr` array (from `m_eri`).

*   **`calculate_eri_4center_shell(...)`**:
    *   **Purpose**: Calculates 4-center ERIs for a single, specific quartet of shells.

*   **`calculate_eri_4center_shell_grad(...)`**:
    *   **Purpose**: Calculates the gradients of 4-center ERIs for a specific quartet of shells with respect to nuclear coordinates. Calls `libint_4center_grad`.

*   **`calculate_eri_ri(basis, auxil_basis, rcut)`**:
    *   **Purpose**: Main driver for calculating ERIs using the Resolution of Identity (RI) approximation.
    *   Calls `calculate_integrals_eri_2center_scalapack` to compute (Q|P) integrals.
    *   Calls either `calculate_inverse_eri_2center_scalapack` (for V<sup>-1</sup>) or `calculate_inverse_sqrt_eri_2center_scalapack` (for V<sup>-1/2</sup>) depending on `eri3_genuine_`.
    *   Calls `calculate_integrals_eri_3center_scalapack` (for genuine (Q|ij)) or `calculate_eri_3center_scalapack` (for RI-transformed (Q|P)<sup>-1/2</sup>(ij|P)).

*   **RI-related Subroutines (ScaLAPACK enabled)**:
    *   **`calculate_integrals_eri_2center_scalapack(auxil_basis, rcut, mask_auxil)`**: Computes the 2-center ERI matrix V<sub>QP</sub> = (Q|1/r<sub>12</sub>|P) over auxiliary basis functions. Handles optional masking for recalculations.
    *   **`calculate_inverse_eri_2center_scalapack(auxil_basis, rcut)`**: Computes the inverse of the 2-center ERI matrix V<sup>-1</sup>. Uses ScaLAPACK's `PDPOTRF`/`PDPOTRI` (if positive definite) or a general inversion.
    *   **`calculate_inverse_sqrt_eri_2center_scalapack(auxil_basis, rcut)`**: Computes V<sup>-1/2</sup> by diagonalizing V, taking `1/sqrt(eigenvalue)` for non-negligible eigenvalues, and transforming back. This is key for standard RI.
    *   **`calculate_integrals_eri_3center_scalapack(basis, auxil_basis, rcut, mask, mask_auxil)`**: Computes "genuine" 3-center ERIs (Q|&mu;&nu;) where Q is an auxiliary basis function and &mu;,&nu; are AO basis functions. Stores results in `eri_3center`.
    *   **`calculate_eri_3center_scalapack(basis, auxil_basis, rcut)`**: Computes RI-transformed 3-center integrals M<sub>Q,(&mu;&nu;)</sub> = &sum;<sub>P</sub> (Q|P)<sup>-1/2</sup> (&mu;&nu;|P). Stores results in `eri_3center`.

*   **`calculate_eri_approximate_hartree(...)`**:
    *   **Purpose**: Calculates an approximate Hartree potential matrix V<sub>&mu;&nu;</sub> by contracting 3-center ERIs (&mu;&nu;|Q) with coefficients of a density fitted to Gaussian functions (Q).

*   **`destroy_eri_3center()`**:
    *   Deallocates module-level arrays associated with 2-center and 3-center ERIs (e.g., `eri_2center`, `eri_2center_inv`, `eri_2center_sqrtinv`, `eri_3center`).

## Important Variables/Constants

*   **`eri_2center(:, :)` (Real(dp), Protected, Allocatable)**: Stores the matrix of 2-center integrals (Q|P) over auxiliary basis functions.
*   **`eri_2center_inv(:, :)` (Real(dp), Protected, Allocatable)**: Stores V<sup>-1</sup>.
*   **`eri_2center_sqrtinv(:, :)` (Real(dp), Protected, Allocatable)**: Stores V<sup>-1/2</sup>.
*   **`eri_2center_sqrt(:, :)` (Real(dp), Protected, Allocatable)**: Stores V<sup>1/2</sup> (experimental/conditional).
*   **`incore_` (Logical, from `m_inputparam`)**: If true, 4-center ERIs are stored in memory. Otherwise, they might be recomputed on the fly or an out-of-core algorithm would be needed (though the latter is not explicitly shown in this module for 4c).
*   **`eri3_genuine_` (Logical, from `m_inputparam`)**: If true, "genuine" 3-center integrals (Q|&mu;&nu;) are used. Otherwise, the RI approximation M<sub>Q,(&mu;&nu;)</sub> = &sum;<sub>P</sub> (Q|P)<sup>-1/2</sup> (&mu;&nu;|P) is formed.
*   **`rcut` (Real(dp), Intent(in))**: Cutoff parameter for range-separated Coulomb interactions. If zero, standard Coulomb interaction is used.
*   **Arrays from `m_eri`**: This module populates `eri_4center`, `eri_4center_lr`, and `eri_3center`, `eri_3center_lr` which are declared in `m_eri`.

## Usage Examples

This module is primarily used internally by MOLGW's main computational workflow.
1.  `m_eri::prepare_eri` is called first to set up screening and indexing.
2.  Then, `calculate_eri` (for 4-center) or `calculate_eri_ri` (for RI-based methods) is called.

```fortran
! Conceptual: In main MOLGW program or a high-level driver
USE m_eri_calculate
USE m_eri
USE m_basis_set
USE m_inputparam

TYPE(basis_set) :: orbital_basis, aux_basis
LOGICAL :: do_print_eri = .FALSE.
REAL(DP) :: rcut_coulomb = 0.0_dp ! Standard Coulomb

! ... (orbital_basis and aux_basis are initialized) ...
! ... (input parameters like incore_, integral_level, eri3_genuine_ are set) ...

CALL prepare_eri(orbital_basis) ! From m_eri

IF (has_auxil_basis) THEN
  CALL calculate_eri_ri(orbital_basis, aux_basis, rcut_coulomb)
ELSE
  CALL calculate_eri(do_print_eri, orbital_basis, rcut_coulomb)
END IF

! Now ERIs are available in arrays from m_eri for subsequent calculations.
```

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_mpi`, `m_scalapack`, `m_memory`, `m_warning`, `m_timing`, `m_inputparam`.
*   **Basis Set & ERI Storage**: `m_basis_set` (for basis/shell information), `m_eri` (for storage arrays like `eri_4center`, `eri_3center`, and indexing schemes like `index_shellpair`, `negligible_shellpair`).
*   **Transformations**: `m_cart_to_pure` (for `transform_libint_to_molgw` and `transform_libcint_to_molgw` which convert integrals from the Cartesian form produced by Libint/Libcint to the desired final representation).
*   **Integral Engines**:
    *   `m_libint_tools`: Provides Fortran wrappers (`libint_4center`, `libint_3center`, `libint_2center`, `libint_4center_grad`) to a native Libint library or similar for performing the actual ERI calculations over primitive Gaussians.
    *   `m_libcint_tools`: (If `HAVE_LIBCINT` is defined) Provides Fortran wrappers (`cint2e_cart`, `cint3c2e_cart`, `cint2c2e_cart`, etc.) to the Libcint library.
*   **ScaLAPACK/BLAS**: Used for matrix operations in `calculate_inverse_sqrt_eri_2center_scalapack` (`diagonalize_sca`, `PDGEMR2D`, `PDSCAL`) and `calculate_eri_3center_scalapack` (`PDGEMM`).
*   This module is responsible for the computationally intensive step of ERI calculation and relies heavily on external (Libint/Libcint) or internal optimized libraries for this. The choice of RI vs. 4-center methods, and the use of screening, significantly impacts performance and memory usage.
```
