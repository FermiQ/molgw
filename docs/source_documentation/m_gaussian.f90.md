# `m_gaussian.f90`

## Overview

The `m_gaussian` Fortran module provides the fundamental data structures and routines for working with individual primitive Cartesian Gaussian functions. A Cartesian Gaussian function is of the form (x-x0)<sup>nx</sup> * (y-y0)<sup>ny</sup> * (z-z0)<sup>nz</sup> * exp(-&alpha; * |**r**-**R<sub>0</sub>**|<sup>2</sup>). This module includes functionalities for initializing these Gaussians, evaluating their values at specific points in real or Fourier space, calculating their gradients and Laplacians, and computing various one-electron integrals (overlap, kinetic energy, nuclear attraction) between two such Gaussians using recurrence relations (Obara-Saika). It also provides routines for generalized oscillator strength (GOS) elements.

## Key Components

*   **`gaussian` (Type)**:
    *   A derived data type that encapsulates all necessary information for a single primitive Cartesian Gaussian function.
    *   `nx, ny, nz` (Integer): Cartesian exponents for x, y, and z directions.
    *   `am` (Integer): Total angular momentum (l = nx+ny+nz).
    *   `amc` (Character(len=1)): Single-letter representation of angular momentum (e.g., 's', 'p', 'd').
    *   `alpha` (Real(dp)): Gaussian exponent.
    *   `x0(3)` (Real(dp)): Coordinates of the Gaussian center.
    *   `norm_factor` (Real(dp)): The normalization constant applied to the Gaussian, defined as `common_norm_factor * (1.0 / sqrt(double_factorial(2*nx-1) * ...))`.
    *   `common_norm_factor` (Real(dp)): Part of the normalization factor: `(2/&pi;)^0.75 * 2^l * &alpha;^((2l+3)/4)`.

*   **`init_gaussian_cart(nx, ny, nz, alpha, x0, ga)` (Subroutine)**:
    *   Initializes a `gaussian` type variable `ga` with the given Cartesian exponents, exponent `alpha`, and center `x0`. It also computes and stores the normalization factors.

*   **`eval_gaussian(ga, x)` (Pure Function)**:
    *   Evaluates the value of the (normalized) Cartesian Gaussian `ga` at a given point `x(3)`.

*   **`eval_gaussian_fourier(ga, k)` (Function)**:
    *   Computes the Fourier transform of the (normalized) Cartesian Gaussian `ga` at a given momentum vector `k(3)`.
    *   The Fourier transform definition used is &phi;(k) = 1/(2&pi;)<sup>3</sup> &int; dr e<sup>-i**k**&middot;**r**</sup> &phi;(**r**).
    *   Contains hardcoded analytical formulas for angular momenta up to l=8 (for even powers of k) and l=7 (for odd powers of k).

*   **`eval_gaussian_grad(ga, x)` (Function)**:
    *   Computes the three Cartesian components of the gradient &nabla;&phi;(**r**) of the Gaussian `ga` at point `x(3)`.

*   **`eval_gaussian_lapl(ga, x)` (Function)**:
    *   Computes components related to second derivatives. The current implementation seems to calculate components of &nabla;(&nabla;&phi;) rather than the scalar Laplacian &nabla;<sup>2</sup>&phi;. Specifically, it returns a 3-component vector where each component `i` is &sum;<sub>j</sub> &part;<sup>2</sup>&phi;/(&part;x<sub>i</sub>&part;x<sub>j</sub>). *Correction: On closer inspection, the formula seems to be for diagonal elements of the Hessian, i.e., &part;<sup>2</sup>&phi;/&part;x<sub>i</sub><sup>2</sup> for each component i.*

*   **`overlap_recurrence(ga, gb, s_ab)` (Subroutine)**:
    *   Calculates the overlap integral S<sub>ab</sub> = <g<sub>a</sub>|g<sub>b</sub>> between two Gaussians `ga` and `gb` using Obara-Saika recurrence relations.

*   **`kinetic_recurrence(ga, gb, k_ab)` (Subroutine)**:
    *   Calculates the kinetic energy integral T<sub>ab</sub> = <g<sub>a</sub>|-1/2 &nabla;<sup>2</sup>|g<sub>b</sub>> using Obara-Saika recurrence relations.

*   **`nucleus_recurrence(zatom, c, ga, gb, v_ab)` (Subroutine)**:
    *   Calculates the nuclear attraction integral V<sub>ab</sub> = <g<sub>a</sub>|Z/|**r**-**C**||g<sub>b</sub>> using Obara-Saika recurrence relations and the Boys function (obtained via an interface to `boys_function_c`).
    *   `zatom`: Nuclear charge Z.
    *   `c(3)`: Coordinates of the nucleus C.

*   **`evaluate_gos(ga, gb, qvec, gos_ab)` (Subroutine)**:
    *   Calculates the generalized oscillator strength (GOS) matrix element <g<sub>a</sub>|e<sup>i**q**&middot;**r**</sup>|g<sub>b</sub>> for a given momentum transfer vector `qvec`.
    *   Uses helper functions `g` (for 1D integrals related to Hermite polynomials) and `f` (for binomial expansions).

*   **`boys_function_c(fnt, n, t)` (Interface)**:
    *   An interface to an external C-compatible function (likely implemented in Fortran elsewhere or in C) that computes Boys function values F<sub>m</sub>(T).

*   **Helper functions**: `norm_factor`, `compare_gaussian`, `print_gaussian`.

## Important Variables/Constants

*   **`ga%norm_factor`**: The full normalization constant for a primitive Cartesian Gaussian. All evaluation and integral routines in this module return values corresponding to *normalized* primitive Gaussians.
*   **`MOLGW_LMAX` (Preprocessor Constant from `molgw.h`)**: Used in `m_cart_to_pure` which `m_gaussian` uses for `double_factorial`.

## Usage Examples

This module is primarily a low-level utility used by `m_basis_set` to construct and operate on contracted basis functions.

```fortran
! Conceptual example of initializing and evaluating a Gaussian:
USE m_gaussian
USE m_definitions
IMPLICIT NONE

TYPE(gaussian) :: my_gaussian
REAL(DP) :: center(3), point_of_eval(3)
REAL(DP) :: alpha_exp, value

center(1) = 0.0_dp; center(2) = 0.0_dp; center(3) = 0.0_dp
alpha_exp = 1.0_dp

! Initialize a px Gaussian: (x-x0)^1 * exp(-alpha*r^2)
CALL init_gaussian_cart(1, 0, 0, alpha_exp, center, my_gaussian)

point_of_eval(1) = 0.5_dp; point_of_eval(2) = 0.0_dp; point_of_eval(3) = 0.0_dp
value = eval_gaussian(my_gaussian, point_of_eval)

WRITE(*,*) "Value of px Gaussian at (0.5, 0, 0): ", value
```

## Dependencies and Interactions

*   **`m_definitions`**: For `dp` (double precision kind), `pi`, and `im` constants.
*   **`m_mpi`**: Imported but not directly used by most functions in this module.
*   **`m_cart_to_pure`**: For the `double_factorial` function, which is used in calculating normalization factors.
*   **`m_string_tools`**: For `orbital_momentum_name`.
*   **External Boys Function**: The `nucleus_recurrence` subroutine depends on an external `boys_function_c` implementation for evaluating Boys functions, which are essential for nuclear attraction integrals. This is typically provided by `boys_function.f90` or a linked math library.
*   **`m_basis_set`**: This module is the primary user of `m_gaussian`. Contracted basis functions in `m_basis_set` are composed of `gaussian` objects defined here. Routines in `m_basis_set` for evaluating contracted functions or their integrals call the corresponding routines in `m_gaussian` for each primitive component.
```
