# `m_libcint_tools.f90`

## Overview

The `m_libcint_tools` Fortran module serves as the primary interface layer between MOLGW and the external Libcint library. Libcint is a C library specialized in the efficient computation of molecular integrals over Gaussian basis functions. This module defines the necessary Fortran data structures that mirror Libcint's expectations, provides interfaces to call various C functions within Libcint for computing one- and two-electron integrals (and their derivatives), and includes routines to transform the raw integral data from Libcint's native ordering and normalization conventions to those used internally by MOLGW. The functionality of this module is conditional on MOLGW being compiled and linked against a compatible version of the Libcint library (indicated by the `HAVE_LIBCINT` preprocessor macro).

## Key Components

*   **Libcint Data Structure Definitions (Parameters)**:
    *   Defines integer parameters (e.g., `LIBCINT_CHARGE_OF`, `LIBCINT_PTR_COORD`, `LIBCINT_ATOM_OF`, `LIBCINT_ANG_OF`, `LIBCINT_PTR_EXP`, `LIBCINT_PTR_ENV_START`) that map Fortran array indices to the specific data layout required by Libcint's `atm` (atom information), `bas` (basis shell information), and `env` (environment containing coordinates, exponents, coefficients) arrays.

*   **External C Function Interfaces**:
    *   A comprehensive set of `INTERFACE` blocks declaring external C functions from the Libcint library. These cover a wide range of integrals:
        *   One-electron: `cint1e_ovlp_cart`, `cint1e_ovlp_sph` (overlap), `cint1e_kin_cart` (kinetic), `cint1e_nuc_cart` (nuclear attraction), `cint1e_rinv_cart` (1/r operator), `cint1e_r_cart` (dipole), `cint1e_rr_cart` (quadrupole), `cint1e_ipovlp_cart` (overlap derivatives), `cint1e_ipkin_cart` (kinetic derivatives), `cint1e_iprinv_cart` (nuclear attraction derivatives), `cint1e_cg_irxp_cart` (angular momentum), `cint1e_giao_irjxp_cart` (GIAO).
        *   Two-electron: `cint2e_cart` (4-center ERIs), `cint2c2e_sph`/`_cart` (2-center ERIs), `cint3c2e_cart`/`_sph` (3-center ERIs), `cint3c1e_cart` (3-center overlap-like).
    *   These interfaces specify the binding name in C and the expected Fortran types corresponding to C types.

*   **`check_capability_libcint(lmax)` (Subroutine)**:
    *   **Purpose**: Probes the linked Libcint library to determine its capabilities.
    *   Sets `lmax` to the minimum of the input `lmax` and Libcint's actual maximum supported angular momentum (`LMAX_LIBCINT`, hardcoded as 8 in this module but ideally should come from Libcint itself).
    *   Determines if Libcint supports range-separated integrals (by testing a 2-center ERI with a known screening parameter). Sets `libcint_has_range_separation`.
    *   Determines the ordering of pure p-type spherical harmonics (px,py,pz vs. py,pz,px) used by Libcint by computing a known p-p overlap. Sets `pypzpx_order`.
    *   Calculates and stores normalization factors (`libcint_pure_norm`) for pure spherical harmonics as returned by Libcint by computing s-s, p-p, d-d, etc., overlaps with unit coefficients.

*   **Libcint Environment Setup**:
    *   **`set_rinv_origin_libcint(x0, env_local)`**: Modifies a local copy of the Libcint environment array `env_local` to set the origin `x0` for r<sup>-1</sup> type integrals (e.g., nuclear attraction).
    *   **`set_erf_screening_length_libcint(basis, rcut)`**: Sets the range-separation parameter &omega; (omega = 1/rcut) in the global Libcint environment array stored within the `basis` object. Handles cases for standard Coulomb (rcut=0) and screened Coulomb.
    *   **`init_libcint(basis1, basis2)`**: Populates the `LIBCINT_atm`, `LIBCINT_bas`, and `LIBCINT_env` arrays within the `basis1` object. These arrays are structured precisely as Libcint expects. It can prepare an environment for one basis set (`basis1`) or two combined basis sets (e.g., `basis1` as orbital basis and `basis2` as auxiliary basis, where shells from `basis2` are appended after `basis1`).
    *   **`destroy_libcint(basis)`**: Deallocates the `LIBCINT_atm`, `LIBCINT_bas`, and `LIBCINT_env` arrays stored in the `basis` object.

*   **Fortran Wrappers for Libcint Calls**:
    *   **`libcint_3center(...)`, `libcint_overlap(...)`, `libcint_kinetic(...)`, `libcint_elecpot(...)`, `libcint_elecpot_grad(...)`, `libcint_gth_projector(...)`**: These subroutines act as Fortran-friendly wrappers. They take shell and primitive information in MOLGW's native format, construct temporary `tmp_atm`, `tmp_bas`, `tmp_env` arrays formatted for Libcint for the specific shells involved in an integral, and then call the corresponding external C Libcint function.

*   **Integral Transformation Routines**:
    *   **`transform_libcint_to_molgw_2d(gaussian_type, am1, am2, array_in, matrix_out)`**: Transforms a 1D array of 2D integrals (e.g., S<sub>&mu;&nu;</sub>) from Libcint's Cartesian output order and normalization to MOLGW's desired order and normalization (which can be Cartesian or Pure Spherical Harmonics).
    *   **`transform_libcint_to_molgw_3d(gaussian_type_left, am1, gaussian_type_right, am2, am3, array_in, matrix_out)`**: Similar for 3D integrals (e.g., 3-center integrals (Q|&mu;&nu;), transforming &mu;&nu; indices).

## Important Variables/Constants

*   **`pypzpx_order` (Logical, Protected)**: Indicates if Libcint uses a (py,pz,px)-like ordering for p-orbitals instead of (px,py,pz). Determined at runtime by `check_capability_libcint`.
*   **`libcint_pure_norm(0:LMAX_LIBCINT)` (Real(dp), Private)**: Stores normalization correction factors for pure spherical harmonics from Libcint.
*   **`LIBCINT_PTR_ENV_START` (Integer, Parameter)**: Defines the starting offset in the `env` array where basis set coefficients and exponents are stored.
*   **`libcint_has_range_separation` (Logical, Protected)**: True if the linked Libcint library was compiled with support for range-separated integrals.

## Usage Examples

This module is not intended for direct use in MOLGW input files. It's a backend component used by other Fortran modules that require molecular integrals.

Conceptual usage by another module (e.g., `m_hamiltonian_onebody`):
```fortran
! In another module, e.g., m_hamiltonian_onebody.f90
USE m_libcint_tools
USE m_basis_set
IMPLICIT NONE

TYPE(basis_set) :: my_basis
REAL(DP) :: overlap_matrix(my_basis%nbf, my_basis%nbf)
! ... my_basis is initialized ...

! Initialize Libcint environment for my_basis
CALL init_libcint(my_basis) 
! (This prepares my_basis%LIBCINT_atm, my_basis%LIBCINT_bas, my_basis%LIBCINT_env)

! To calculate overlap for a pair of shells (ishell, jshell from my_basis):
INTEGER(C_INT) :: am_i, contr_i, am_j, contr_j
REAL(C_DOUBLE) :: center_i(3), center_j(3)
REAL(C_DOUBLE), ALLOCATABLE :: alpha_i(:), coeffs_i(:), alpha_j(:), coeffs_j(:)
REAL(C_DOUBLE), ALLOCATABLE :: s_ij_cart(:) ! Integrals in Cartesian representation
REAL(DP), ALLOCATABLE :: s_ij_transformed(:,:) ! Transformed integrals

! 1. Extract shell data for ishell and jshell into Libcint-compatible format
CALL set_libint_shell(my_basis%shell(ishell), am_i, contr_i, center_i, alpha_i, coeffs_i) ! (A helper not in m_libcint_tools)
CALL set_libint_shell(my_basis%shell(jshell), am_j, contr_j, center_j, alpha_j, coeffs_j) ! (A helper not in m_libcint_tools)

ALLOCATE(s_ij_cart(num_cart_i * num_cart_j))

! 2. Call the Fortran wrapper which calls the C Libcint function
CALL libcint_overlap(am_i, contr_i, center_i, alpha_i, coeffs_i, &
                     am_j, contr_j, center_j, alpha_j, coeffs_j, &
                     s_ij_cart)
! Note: The actual libcint_overlap in this file sets up its own temporary env.
! The global basis%LIBCINT_env is used by direct cint1e_ovlp_cart calls.

! 3. Transform to MOLGW's convention
CALL transform_libcint_to_molgw_2d(my_basis%gaussian_type, am_i, am_j, s_ij_cart, s_ij_transformed)
! ... store s_ij_transformed into the global overlap_matrix ...

DEALLOCATE(alpha_i, coeffs_i, alpha_j, coeffs_j, s_ij_cart, s_ij_transformed)
! ...
CALL destroy_libcint(my_basis) ! At the very end
```

## Dependencies and Interactions

*   **Libcint C Library**: This module is critically dependent on the external Libcint library. All `cint*_cart`/`cint*_sph` functions are external C functions from this library.
*   **`m_definitions`**: For `dp` and C- interoperability kinds (`C_INT`, `C_DOUBLE`).
*   **`m_cart_to_pure`**: For `double_factorial` (used in normalization constants) and potentially for transformation matrices if Libcint output needs further processing not covered by `transform_libcint_to_molgw_*`.
*   **`m_basis_set`**: For the `basis_set` derived type, which is augmented by this module to include `LIBCINT_atm`, `LIBCINT_bas`, and `LIBCINT_env` members. Also uses `number_basis_function_am`.
*   **`m_atoms`**: For atomic coordinates and charges used in `init_libcint`.
*   **Modules Calculating Integrals**: Modules like `m_hamiltonian_onebody` and `m_eri_calculate` will call the Fortran wrapper subroutines provided by `m_libcint_tools` (e.g., `libcint_overlap`, `libcint_3center`) or directly use the C interfaces after setting up the global Libcint environment via `init_libcint`.
```
