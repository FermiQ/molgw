# `libint_twobody_grad.cpp`

## Overview

This C++ source file is part of MOLGW and provides wrapper functions to compute the **gradients** of two-body electron repulsion integrals (ERIs) with respect to atomic coordinates, utilizing the `libint2` library. It includes functions for 2-center, 3-center, and 4-center ERI gradients, which are essential for calculating atomic forces in molecular dynamics or geometry optimization procedures. These functions are designed with `extern "C"` linkage for Fortran compatibility and handle contracted Gaussian basis functions. The `rcut` parameter indicates support for screened or range-separated Coulomb interactions.

The compilation of these gradient functions is conditional on `libint2` being available and configured to support the necessary first-order ERI derivatives.

**Important Note**: The `libint_2center_grad` and `libint_3center_grad` functions are currently **not implemented** and contain `assert(false)` statements, meaning they will cause program termination if called. Only `libint_4center_grad` is functional.

## Key Components

*   **`libint_2center_grad(...)`** (Not Implemented):
    *   **Purpose (Intended)**: To compute the gradient of 2-center ERIs <&phi;<sub>A</sub>| 1/r<sub>12</sub> |&phi;<sub>C</sub>> with respect to the coordinates of centers A and C.
    *   **Status**: Contains an `assert(false)` and will abort if called. The code structure is present but appears to be a copy of the non-gradient `libint_2center` function and is not correctly set up for gradient calculations.

*   **`libint_3center_grad(...)`** (Not Implemented):
    *   **Purpose (Intended)**: To compute the gradient of 3-center ERIs <&phi;<sub>A</sub>| 1/r<sub>12</sub> |&phi;<sub>C</sub>&phi;<sub>D</sub>> with respect to the coordinates of centers A, C, and D.
    *   **Status**: Contains an `assert(false)` and will abort if called. Similar to `libint_2center_grad`, the body seems to be a placeholder.

*   **`libint_4center_grad(amA, ..., cA, amB, ..., cB, amC, ..., cC, amD, ..., cD, rcut, eriAx, eriAy, eriAz, eriBx, ..., eriBz, eriCx, ..., eriCz, eriDx, ..., eriDz)`**:
    *   **Purpose**: Computes the gradient of the 4-center, 2-electron ERIs (AB|CD) = <&phi;<sub>A</sub>&phi;<sub>B</sub>| 1/r<sub>12</sub> |&phi;<sub>C</sub>&phi;<sub>D</sub>> with respect to the Cartesian coordinates of all four atomic centers A, B, C, and D.
    *   **Parameters**:
        *   `amA`, `amB`, `amC`, `amD`: Angular momenta.
        *   `contrdepthA`, ..., `contrdepthD`: Contraction depths.
        *   `A[]`, `B[]`, `C[]`, `D[]`: Coordinates of centers.
        *   `alphaA[]`, ..., `alphaD[]`: Primitive exponents.
        *   `cA[]`, ..., `cD[]`: Contraction coefficients.
        *   `rcut` (Input, `double`): Cutoff radius for screened Coulomb interaction.
        *   `eriAx[], eriAy[], eriAz[]` (Output, `double[]`): Gradient components with respect to center A's x, y, z coordinates.
        *   `eriBx[], eriBy[], eriBz[]` (Output, `double[]`): Gradient components with respect to center B's x, y, z coordinates.
        *   `eriCx[], eriCy[], eriCz[]` (Output, `double[]`): Gradient components with respect to center C's x, y, z coordinates.
        *   `eriDx[], eriDy[], eriDz[]` (Output, `double[]`): Gradient components with respect to center D's x, y, z coordinates.

**Common Internal Workflow for `libint_4center_grad`**:
1.  Initializes `libint2` structures for first-derivative 4-center ERIs (`Libint_eri1_t`).
2.  Loops over all combinations of primitive Gaussians.
3.  For each primitive combination:
    *   Calculates product Gaussian parameters.
    *   Sets up `Libint_eri1_t` with geometric data and terms specific to derivative calculations (e.g., `alpha1_rho_over_zeta2`).
    *   Adjusts interaction parameters (`gammapq_rc2`, `gammapq_ratio`) based on `rcut`.
    *   Calls `boys_function_c` for F<sub>m</sub>(U) values (up to order `am+1` as derivatives require higher-order Boys functions).
    *   Scales Boys function results by factors including `gammapq_ratio`.
4.  Invokes `libint2_build_eri1` to compute the full tensor of ERI first derivatives.
5.  Extracts and stores the 12 sets of Cartesian gradient components (d/dAx, d/dAy, d/dAz, d/dBx, ..., d/dDz) into the respective output arrays.
6.  Cleans up `libint2` resources.

## Important Variables/Constants

*   **Input Parameters**: As described for each function. `rcut` is crucial for screened interactions.
*   **Output Arrays**: For `libint_4center_grad`, 12 arrays store the gradient components of the ERI tensor with respect to the x,y,z coordinates of each of the four centers.
*   **`Libint_...1_t` (`inteval`)**: `libint2` data structures for first-derivative ERI computations (e.g., `Libint_eri1_t`).
*   **`LIBINT2_MAX_AM_...1`**: `libint2` constants for maximum angular momentum for derivative integrals.
*   `ni`: Total number of integral components for a given set of shell angular momenta.

## Usage Examples

These functions are intended for use by MOLGW's Fortran components, particularly in modules calculating atomic forces.

Conceptual Fortran call for 4-center ERI gradients:
```fortran
! Assume basis parameters (am, contrdepth, coords, etc.) for A,B,C,D and rcut_val are defined
DOUBLE PRECISION, ALLOCATABLE :: GRAD_AX(:), GRAD_AY(:), GRAD_AZ(:), &
                                 GRAD_BX(:), GRAD_BY(:), GRAD_BZ(:), &
                                 GRAD_CX(:), GRAD_CY(:), GRAD_CZ(:), &
                                 GRAD_DX(:), GRAD_DY(:), GRAD_DZ(:)
INTEGER :: NUM_INTS_TOTAL

! NUM_INTS_TOTAL = NINT_FUNC(AMA)*NINT_FUNC(AMB)*NINT_FUNC(AMC)*NINT_FUNC(AMD)
ALLOCATE(GRAD_AX(NUM_INTS_TOTAL), ..., GRAD_DZ(NUM_INTS_TOTAL)) ! For all 12 output arrays

CALL LIBINT_4CENTER_GRAD(AMA, ..., CDA, AMB, ..., CDB, &
                         AMC, ..., CDC, AMD, ..., CDD, &
                         RCUT_VAL, &
                         GRAD_AX, GRAD_AY, GRAD_AZ, &
                         GRAD_BX, GRAD_BY, GRAD_BZ, &
                         GRAD_CX, GRAD_CY, GRAD_CZ, &
                         GRAD_DX, GRAD_DY, GRAD_DZ)
! Output arrays now contain the ERI gradient components.
```

## Dependencies and Interactions

*   **`libint2` Library**: The fundamental dependency for the ERI gradient calculations. Uses `libint2` API for first-derivative ERIs.
*   **`molgw.h`**: General MOLGW project definitions.
*   **`libint_molgw.h`**: Provides `boys_function_c` declaration, `nint(am)` function, and `pi_2p5` constant.
*   **`boys_function_c`**: External function (likely Fortran) used for Boys function evaluations, which are needed for the underlying ERI calculations and their derivatives.
*   **Standard C/C++ Libraries**: `stdlib.h`, `stdio.h`, `iostream`, `math.h`, `assert.h`.
*   **Conditional Compilation**:
    *   `#if !defined(NO_LIBINT)`: Guards the entire file.
    *   `#if (LIBINT2_DERIV_ERI_ORDER > 0) && (LIBINT2_DERIV_ERI2_ORDER > 0) && (LIBINT2_DERIV_ERI3_ORDER > 0)`: The ERI gradient functions are only compiled if `libint2` supports first derivatives for 4-center, 2-center, and 3-center ERIs respectively. The `assert(false)` in the 2-center and 3-center gradient functions makes them unusable regardless of these flags until implemented.

This file provides the necessary interface for MOLGW to compute ERI gradients, which are essential for simulations involving changes in nuclear geometry.
```
