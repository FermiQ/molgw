# `libint_onebody_grad.cpp`

## Overview

This C++ source file is a component of MOLGW that provides wrapper functions to compute the **gradients** of one-body molecular integrals using the `libint2` library. Specifically, it calculates the derivatives of overlap, kinetic energy, and nuclear attraction integrals with respect to the Cartesian coordinates of the atomic centers. These gradients are crucial for calculating atomic forces in molecular simulations. The functions are designed with `extern "C"` linkage for compatibility with Fortran code and handle contracted Gaussian basis functions.

The compilation and availability of these functions are conditional upon `libint2` being present, configured to support one-body integrals, and compiled with support for at least first-order derivatives (`LIBINT2_DERIV_ONEBODY_ORDER > 0`).

## Key Components

The file defines three primary functions for calculating the gradients of different one-body integrals:

*   **`libint_overlap_grad(amA, contrdepthA, A, alphaA, cA, amB, contrdepthB, B, alphaB, cB, overlapABx, overlapABy, overlapABz)`**:
    *   **Purpose**: Computes the gradient of the overlap integrals &nabla;<sub>R</sub> S<sub>&mu;&nu;</sub> = &nabla;<sub>R</sub> <&phi;<sub>&mu;</sub> | &phi;<sub>&nu;</sub>> with respect to the Cartesian coordinates of the centers of the basis functions &phi;<sub>&mu;</sub> (at A) and &phi;<sub>&nu;</sub> (at B). `libint2` typically provides these as dS/dAx, dS/dAy, dS/dAz.
    *   **Parameters**:
        *   `amA`, `amB` (Input, `int`): Angular momenta of basis functions.
        *   `contrdepthA`, `contrdepthB` (Input, `int`): Contraction depths.
        *   `A[]`, `B[]` (Input, `double[3]`): Coordinates of centers A and B.
        *   `alphaA[]`, `alphaB[]` (Input, `double[]`): Primitive exponents.
        *   `cA[]`, `cB[]` (Input, `double[]`): Contraction coefficients.
        *   `overlapABx[]`, `overlapABy[]`, `overlapABz[]` (Output, `double[]`): Arrays to store the x, y, and z components of the gradient of the overlap integrals. `libint2` returns derivatives organized by d/d(center1_coord), d/d(center2_coord). This wrapper seems to extract the d/d(center_A) components.

*   **`libint_kinetic_grad(amA, contrdepthA, A, alphaA, cA, amB, contrdepthB, B, alphaB, cB, kineticABx, kineticABy, kineticABz)`**:
    *   **Purpose**: Computes the gradient of the kinetic energy integrals &nabla;<sub>R</sub> T<sub>&mu;&nu;</sub> = &nabla;<sub>R</sub> <&phi;<sub>&mu;</sub> | -1/2 &nabla;<sup>2</sup> | &phi;<sub>&nu;</sub>>.
    *   **Parameters**: Similar to `libint_overlap_grad`, with `kineticABx[]`, `kineticABy[]`, `kineticABz[]` as output arrays for the gradient components.

*   **`libint_elecpot_grad(amA, contrdepthA, A, alphaA, cA, amB, contrdepthB, B, alphaB, cB, C, elecpotAx, elecpotAy, elecpotAz, elecpotBx, elecpotBy, elecpotBz)`**:
    *   **Purpose**: Computes the gradient of the nuclear attraction integrals &nabla;<sub>R</sub> V<sub>&mu;&nu;</sub> = &nabla;<sub>R</sub> <&phi;<sub>&mu;</sub> | 1/|**r** - **R**<sub>C</sub>| | &phi;<sub>&nu;</sub>> with respect to the coordinates of centers A, B (and implicitly C, handled by `libint2`).
    *   **Parameters**: Similar to `libint_overlap_grad`, with the addition of:
        *   `C[]` (Input, `double[3]`): Coordinates of the potential center C.
        *   `elecpotAx[]`, `elecpotAy[]`, `elecpotAz[]` (Output, `double[]`): Arrays for gradient components with respect to center A's coordinates.
        *   `elecpotBx[]`, `elecpotBy[]`, `elecpotBz[]` (Output, `double[]`): Arrays for gradient components with respect to center B's coordinates.
        (Note: `libint2` provides dV/dCx, dV/dCy, dV/dCz as well, which are not explicitly exposed as separate output arrays by this wrapper but are part of the `inteval[0].targets[0]` data structure from which specific derivatives are extracted).

**Common Internal Workflow**:
The internal workflow for each gradient function mirrors that of the integral calculation functions in `libint_onebody.cpp` but is adapted for first-order derivatives:
1.  Initializes `libint2` structures for first derivatives (e.g., `Libint_overlap1_t`).
2.  Loops over primitive Gaussian pairs, setting up geometric and basis parameters.
3.  Calls the appropriate `libint2_build_...1` routine (e.g., `libint2_build_overlap1`) to compute all Cartesian components of the integral derivatives.
4.  Extracts and stores the relevant gradient components (typically d/dAx, d/dAy, d/dAz for the first center, and similarly for the second center in `elecpot_grad`) into the output arrays.
5.  Cleans up `libint2` structures.

## Important Variables/Constants

*   **Input Parameters**: Same as those in `libint_onebody.cpp` (angular momenta, contraction data, coordinates, exponents, coefficients).
*   **Output Arrays**:
    *   `overlapABx, overlapABy, overlapABz`: For x,y,z components of overlap gradients.
    *   `kineticABx, kineticABy, kineticABz`: For x,y,z components of kinetic energy gradients.
    *   `elecpotAx, ..., elecpotBz`: For x,y,z components of nuclear attraction gradients with respect to centers A and B.
*   **`Libint_...1_t` (`inteval`)**: `libint2` data structures for first-derivative integral computations.
*   **`LIBINT2_MAX_AM_...1`**: `libint2` constants for maximum supported angular momentum for first-derivative integrals.
*   `ni`: Calculated as `nint(amA) * nint(amB)`, representing the number of integral components for a given pair of shell angular momenta.

## Usage Examples

These C++ functions are intended for invocation from Fortran code within MOLGW, likely when computing atomic forces or optimizing molecular geometries.

Conceptual Fortran call for overlap gradients:
```fortran
! Assuming necessary variables (AM_A, NPRIM_A, COORDS_A, etc.) are defined
DOUBLE PRECISION, ALLOCATABLE :: OVERLAP_GRAD_X(:), OVERLAP_GRAD_Y(:), OVERLAP_GRAD_Z(:)
INTEGER :: NUM_CART_A, NUM_CART_B, NUM_INTS

! NUM_CART_A = NINT_FUNC(AM_A)
! NUM_CART_B = NINT_FUNC(AM_B)
NUM_INTS = NUM_CART_A * NUM_CART_B
ALLOCATE(OVERLAP_GRAD_X(NUM_INTS), OVERLAP_GRAD_Y(NUM_INTS), OVERLAP_GRAD_Z(NUM_INTS))

CALL LIBINT_OVERLAP_GRAD(AM_A, NPRIM_A, COORDS_A, EXP_A, COEFF_A, &
                         AM_B, NPRIM_B, COORDS_B, EXP_B, COEFF_B, &
                         OVERLAP_GRAD_X, OVERLAP_GRAD_Y, OVERLAP_GRAD_Z)
! Output arrays now contain the gradient components of the overlap integrals.
```

## Dependencies and Interactions

*   **`libint2` Library**: This is the core dependency. The code uses `libint2` data structures and API functions specifically designed for first-order integral derivatives (e.g., `Libint_overlap1_t`, `libint2_init_overlap1`, `libint2_build_overlap1`).
*   **`molgw.h`**: For general MOLGW project settings.
*   **`libint_molgw.h`**: Provides:
    *   Declaration for `boys_function_c` (used by `libint_elecpot_grad`).
    *   The `nint(am)` inline function.
*   **`boys_function_c`**: The `libint_elecpot_grad` function relies on `boys_function_c` for Boys function values, which are needed for computing the potential integrals and their derivatives.
*   **Standard C/C++ Libraries**: `stdlib.h`, `stdio.h`, `iostream`, `math.h`, `assert.h`.
*   **Conditional Compilation**:
    *   `#if !defined(NO_LIBINT)`: Guards the entire file.
    *   `#if defined(LIBINT2_SUPPORT_ONEBODY) && (LIBINT2_DERIV_ONEBODY_ORDER > 0)`: Ensures these gradient functions are only compiled if `libint2` is available, supports one-body integrals, and is configured to compute their first derivatives. This is critical as derivative computation capabilities can be an optional part of `libint2`.

This file extends the one-body integral interface to include the calculation of their gradients, which are essential for geometry optimizations and molecular dynamics simulations where forces on nuclei are required.
```
