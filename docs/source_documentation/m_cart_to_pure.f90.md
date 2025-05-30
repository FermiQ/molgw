# `m_cart_to_pure.f90`

## Overview

The `m_cart_to_pure` Fortran module provides the necessary data structures and routines to facilitate the transformation of Gaussian basis functions (or related quantities like integrals or coefficients) between their Cartesian representation and their pure spherical harmonic representation. This is a common requirement in quantum chemistry software as different integral libraries or computational methods may prefer one representation over the other. The module pre-calculates and stores these transformation matrices up to a defined maximum angular momentum (`MOLGW_LMAX`).

## Key Components

*   **`transform` (Type)**:
    *   A derived data type that encapsulates a transformation matrix.
    *   `matrix(:, :)` (Real(dp), Allocatable): A 2D array holding the coefficients for the transformation.

*   **`cart_to_pure(0:MOLGW_LMAX, 2)` (Type(transform), Allocatable, Protected)**:
    *   A 2D array where each element `cart_to_pure(l, type)` stores a `transform` object.
    *   `l`: Angular momentum.
    *   `type`: Either `CARTG` (for Cartesian-to-Cartesian, which is an identity matrix) or `PUREG` (for Cartesian-to-Pure Spherical Harmonics).
    *   `cart_to_pure(l, PUREG)%matrix` stores the coefficients C<sub>pure,m</sub> = &sum;<sub>xyz</sub> T<sub>xyz,m</sub> C<sub>cart,xyz</sub>. The matrix dimensions are (number_of_cartesian_functions, number_of_pure_functions).

*   **`cart_to_pure_norm(0:MOLGW_LMAX, 2)` (Type(transform), Allocatable, Protected)**:
    *   Similar to `cart_to_pure`, but these matrices are intermediate forms used during the setup, likely related to unnormalized or differently normalized transformation coefficients before final scaling.

*   **`CARTG`, `PUREG` (Integer, Parameters)**:
    *   Constants (1 and 2 respectively) used as symbolic indices for the second dimension of `cart_to_pure` and `cart_to_pure_norm` to specify the target Gaussian type.

*   **`get_gaussian_type_tag(gaussian_type)` (Function)**:
    *   Converts a character string representation of Gaussian type ('CART' or 'PURE') into its corresponding integer tag (`CARTG` or `PUREG`).

*   **`number_basis_function_am(gaussian_type, am)` (Pure Function)**:
    *   Calculates the number of basis functions for a given angular momentum `am` and `gaussian_type`.
    *   For 'CART': `(am+1)*(am+2)/2`.
    *   For 'PURE': `2*am+1`.

*   **`double_factorial(intin)` (Pure Function)**:
    *   Computes the double factorial (n!!) for a given integer `intin`. Uses a hardcoded `SELECT CASE` for values up to 31. Returns 1.0 for -1 and 0.

*   **`setup_cart_to_pure_transforms(pypzpx_order_in)` (Subroutine)**:
    *   The main initialization routine. It computes the transformation coefficients based on standard mathematical formulas involving binomial coefficients (`cnk`) and permutation-related terms (`ank`).
    *   Populates the `cart_to_pure` and `cart_to_pure_norm` arrays for all angular momenta from 0 to `MOLGW_LMAX`.
    *   `pypzpx_order_in` (Logical): An input flag that adjusts the p-orbital transformation to match a specific ordering convention (e.g., if the pure p-orbitals are ordered as p<sub>y</sub>, p<sub>z</sub>, p<sub>x</sub> instead of the more standard p<sub>x</sub>, p<sub>y</sub>, p_<sub>z</sub> relative to the Cartesian input order).

*   **`destroy_cart_to_pure_transforms()` (Subroutine)**:
    *   Deallocates the memory occupied by the `cart_to_pure` and `cart_to_pure_norm` arrays.

*   **`cnk(n, k)` (Pure Function)**:
    *   Calculates the binomial coefficient C(n,k) = n! / (k! * (n-k)!).

*   **`ank(n, k)` (Pure Function)**:
    *   Calculates the product n * (n-1) * ... * (k+1). This is related to permutations P(n, n-k) = n! / k!.

## Important Variables/Constants

*   **`MOLGW_LMAX` (Integer, Preprocessor Constant)**: Defined in `molgw.h`, this specifies the maximum angular momentum for which transformation matrices are pre-calculated and stored.
*   **`pypzpx_order_in` (Logical)**: Input to `setup_cart_to_pure_transforms`. If true, implies a non-standard ordering for p-orbitals (y, z, x) relative to the input Cartesian order when constructing the pure spherical versions. The default seems to be x,y,z for pure p-orbitals if this is false.

## Usage Examples

The transformation matrices stored in this module are used internally by MOLGW when converting quantities like molecular orbital coefficients or integrals from a Cartesian basis representation to a pure spherical harmonic basis representation, or vice-versa (though the inverse transform is not explicitly stored, it can be derived).

```fortran
USE m_cart_to_pure
IMPLICIT NONE

REAL(DP), ALLOCATABLE :: cart_coeffs(:), pure_coeffs(:)
INTEGER :: l, n_cart_l, n_pure_l
LOGICAL :: p_orbital_custom_order = .FALSE. ! Example value

! Initialize transformations (usually done once at the start of a program)
CALL setup_cart_to_pure_transforms(p_orbital_custom_order)

l = 2 ! Example: d-functions
n_cart_l = number_basis_function_am('CART', l) ! Should be 6
n_pure_l = number_basis_function_am('PURE', l) ! Should be 5

ALLOCATE(cart_coeffs(n_cart_l))
ALLOCATE(pure_coeffs(n_pure_l))

! Assume cart_coeffs is filled with coefficients of Cartesian d-functions
! (e.g., d_xx, d_xy, d_xz, d_yy, d_yz, d_zz in a specific order)
! ... fill cart_coeffs ...

! To transform Cartesian coefficients to Pure spherical harmonic coefficients:
! C_pure(k) = SUM_i T(i,k) * C_cart(i)
! where T is cart_to_pure(l, PUREG)%matrix, with dimensions (n_cart_l, n_pure_l)
! This means the operation is pure_coeffs = MATMUL(TRANSPOSE(cart_to_pure(l, PUREG)%matrix), cart_coeffs)
! if cart_coeffs is a column vector.
! If cart_coeffs is a row vector, then pure_coeffs = MATMUL(cart_coeffs, cart_to_pure(l, PUREG)%matrix)

! Example assuming column vectors for coefficients:
! Let T_matrix = cart_to_pure(l, PUREG)%matrix
! pure_coeffs(1:n_pure_l) = MATMUL(TRANSPOSE(T_matrix), cart_coeffs(1:n_cart_l))

! ... use pure_coeffs ...

DEALLOCATE(cart_coeffs, pure_coeffs)

! Clean up transformations at the end of the program
CALL destroy_cart_to_pure_transforms()
```

## Dependencies and Interactions

*   **`m_definitions`**: For the `dp` double precision kind parameter.
*   **`m_warning`**: For the `die` error handling subroutine.
*   **`molgw.h`**: Provides the `MOLGW_LMAX` preprocessor constant, which dictates the upper limit of angular momentum for which transformations are generated.
*   **Internal Logic**: The core of the module lies in the mathematical formulas used in `setup_cart_to_pure_transforms` to derive the transformation coefficients. These formulas are standard in quantum chemistry for relating Cartesian Gaussians to their pure spherical harmonic counterparts.
*   **Other MOLGW Modules**: Modules dealing with basis set representation (`m_basis_set`), integral calculations, or coefficient manipulation might use the transformation matrices provided by `m_cart_to_pure` to ensure consistency or to interface with libraries that expect a specific representation.
```
