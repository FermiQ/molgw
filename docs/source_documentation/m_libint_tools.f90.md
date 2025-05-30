# `m_libint_tools.f90`

## Overview

The `m_libint_tools` Fortran module serves as an interface layer to the Libint library (specifically, it appears to target Libint2 based on include paths and parameter names like `LIBINT2_SUPPORT_ONEBODY`). Libint is a high-performance C++ library for evaluating molecular electron repulsion integrals and other quantum chemical integrals over Gaussian basis functions. This module provides Fortran-callable C-bindings to wrapped Libint functions and includes utility routines to transform the output from Libint (which is typically in a specific Cartesian ordering and normalization) into the conventions used by MOLGW (which may involve pure spherical harmonics and different normalization).

## Key Components

*   **C Function Interfaces (`INTERFACE` blocks)**:
    *   This module declares a series of Fortran interfaces to C functions that, in turn, call the Libint library. These interfaces are conditionally compiled based on Libint's capabilities (e.g., `LIBINT2_SUPPORT_ONEBODY`, `LIBINT2_DERIV_ONEBODY_ORDER > 0`).
    *   **`libint_init(ammax, has_onebody, has_gradient)`**: Initializes the Libint library and returns its capabilities, such as the maximum angular momentum (`ammax`) supported.
    *   **One-Electron Integrals**:
        *   `libint_overlap(...)`: Computes overlap integrals.
        *   `libint_kinetic(...)`: Computes kinetic energy integrals.
        *   `libint_elecpot(...)`: Computes nuclear attraction integrals.
    *   **Gradients of One-Electron Integrals**:
        *   `libint_overlap_grad(...)`: Gradients of overlap integrals.
        *   `libint_kinetic_grad(...)`: Gradients of kinetic energy integrals.
        *   `libint_elecpot_grad(...)`: Gradients of nuclear attraction integrals.
    *   **Two-Electron Integrals (ERIs)**:
        *   `libint_2center(...)`: Computes 2-center ERIs.
        *   `libint_3center(...)`: Computes 3-center ERIs.
        *   `libint_4center(...)`: Computes 4-center ERIs.
    *   **Gradients of Two-Electron Integrals**:
        *   `libint_4center_grad(...)`: Gradients of 4-center ERIs.

*   **Transformation Subroutines**:
    *   **`transform_molgw_to_molgw_2d(...)`**: This routine's name might be slightly misleading in context. It appears to take data already transformed from Libint's Cartesian output (or directly computed in MOLGW's Cartesian representation) and applies further transformations, likely from Cartesian to pure spherical harmonics using matrices from `m_cart_to_pure`.
    *   **`transform_libint_to_molgw` (Generic Interface)**:
        *   `transform_libint_to_molgw_2d(...)`: Transforms 2D arrays of integrals (e.g., one-electron integrals over shells <i|op|j>) from Libint's Cartesian output to MOLGW's convention.
        *   `transform_libint_to_molgw_3d(...)`: Transforms 3D arrays of integrals (e.g., 3-center ERIs (ij|k)) from Libint's output.
        *   `transform_libint_to_molgw_4d(...)`: Transforms 4D arrays of integrals (e.g., 4-center ERIs (ij|kl)) from Libint's output.
        These routines use normalization factors and transformation matrices from `m_cart_to_pure` to handle the conversion to MOLGW's desired basis type (Cartesian or Pure Spherical Harmonics) and normalization.
    *   **`transform_libint_to_molgw_gth_projector(...)`**: A specialized transformation for GTH pseudopotential projector integrals, where one side of the integral (the projector from Libint) is treated as pure spherical harmonic and normalized accordingly.

*   **`set_libint_shell(shell, amA, contrdepthA, A, alphaA, cA)` (Subroutine)**:
    *   A utility to convert shell data from MOLGW's `shell_type` (defined in `m_basis_set`) into separate arrays (`amA`, `contrdepthA`, `A`, `alphaA`, `cA`) as expected by the C-wrapped Libint functions.

## Important Variables/Constants

*   The module uses C interoperability kinds from `ISO_C_BINDING` (e.g., `C_INT`, `C_DOUBLE`) in its interface declarations.
*   Preprocessor macros like `LIBINT2_SUPPORT_ONEBODY`, `LIBINT2_DERIV_ONEBODY_ORDER`, `LIBINT2_DERIV_ERI_ORDER` (from `libint2/libint2_params.h`) control the conditional compilation of specific interfaces, ensuring that MOLGW only attempts to use Libint features that were enabled during Libint's compilation.

## Usage Examples

This module serves as a backend and is not directly called by users. Other MOLGW modules that compute integrals (like `m_hamiltonian_onebody` for one-electron integrals or `m_eri_calculate` for two-electron integrals) would use these routines when Libcint is not available or when Libint is preferred for certain integral types.

Conceptual flow for calculating overlap integrals for two shells (`shell_i`, `shell_j`) from a `basis_set` object `bs`:
```fortran
USE m_libint_tools
USE m_basis_set
IMPLICIT NONE

TYPE(basis_set) :: bs
TYPE(shell_type) :: shell_i, shell_j
INTEGER :: i_shell_idx, j_shell_idx
! ... (bs, i_shell_idx, j_shell_idx are initialized) ...
shell_i = bs%shell(i_shell_idx)
shell_j = bs%shell(j_shell_idx)

INTEGER(C_INT) :: am_i_c, cd_i_c, am_j_c, cd_j_c
REAL(C_DOUBLE) :: center_i_c(3), center_j_c(3)
REAL(C_DOUBLE), ALLOCATABLE :: alpha_i_c(:), coeffs_i_c(:), alpha_j_c(:), coeffs_j_c(:)
REAL(C_DOUBLE), ALLOCATABLE :: s_ij_libint_cart(:) ! 1D array from libint (Cartesian)
REAL(DP), ALLOCATABLE :: s_ij_molgw(:,:)          ! 2D array for MOLGW (final representation)
INTEGER :: num_cart_i, num_cart_j, num_final_i, num_final_j

! 1. Prepare shell data for Libint C-wrapper
CALL set_libint_shell(shell_i, am_i_c, cd_i_c, center_i_c, alpha_i_c, coeffs_i_c)
CALL set_libint_shell(shell_j, am_j_c, cd_j_c, center_j_c, alpha_j_c, coeffs_j_c)

! 2. Allocate space for Libint output (Cartesian integrals)
num_cart_i = number_basis_function_am('CART', shell_i%am)
num_cart_j = number_basis_function_am('CART', shell_j%am)
ALLOCATE(s_ij_libint_cart(num_cart_i * num_cart_j))

! 3. Call the Fortran interface to the C-wrapped Libint function
CALL libint_overlap(am_i_c, cd_i_c, center_i_c, alpha_i_c, coeffs_i_c, &
                    am_j_c, cd_j_c, center_j_c, alpha_j_c, coeffs_j_c, &
                    s_ij_libint_cart)

! 4. Transform integrals to MOLGW's convention (e.g., pure spherical, different normalization)
num_final_i = number_basis_function_am(bs%gaussian_type, shell_i%am)
num_final_j = number_basis_function_am(bs%gaussian_type, shell_j%am)
ALLOCATE(s_ij_molgw(num_final_i, num_final_j))

CALL transform_libint_to_molgw(bs%gaussian_type, shell_i%am, shell_j%am, &
                               s_ij_libint_cart, s_ij_molgw)
                               
! ... s_ij_molgw now contains the transformed integrals ...

DEALLOCATE(alpha_i_c, coeffs_i_c, alpha_j_c, coeffs_j_c, s_ij_libint_cart, s_ij_molgw)
```

## Dependencies and Interactions

*   **Libint (Libint2) C++ Library**: This module is fundamentally a bridge to an external Libint library. The C functions declared in the interfaces are wrappers around actual Libint calls (likely implemented in separate C++ wrapper files like `libint_onebody.cpp`, `libint_twobody.cpp`, etc., within MOLGW's source). This dependency is managed by the `NO_LIBINT` preprocessor flag.
*   **`m_definitions`**: For `dp` (double precision kind) and C interoperability kinds (`C_INT`, `C_DOUBLE`).
*   **`m_cart_to_pure`**: Crucial for the `transform_libint_to_molgw_*` routines, providing transformation matrices (`cart_to_pure_norm`) and functions like `get_gaussian_type_tag` and `number_basis_function_am`.
*   **`m_basis_set`**: For the `shell_type` definition used in `set_libint_shell` and for `basis%gaussian_type`.
*   **Integral Calculation Modules**: Modules like `m_hamiltonian_onebody` (for one-electron integrals) and `m_eri_calculate` (for two-electron integrals) will invoke the routines in `m_libint_tools` if Libcint is not used or if Libint is specifically chosen for a certain integral type. These modules then use the transformed integrals from `m_libint_tools`.
```
