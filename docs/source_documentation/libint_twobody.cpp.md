# `libint_twobody.cpp`

## Overview

This C++ source file provides a suite of wrapper functions for the `libint2` library, enabling the calculation of two-body electron repulsion integrals (ERIs) within MOLGW. It supports 2-center, 3-center, and 4-center ERIs involving contracted Gaussian-type orbitals (GTOs). The functions are exposed with C linkage (`extern "C"`) for compatibility with Fortran. A notable feature is the inclusion of an `rcut` parameter, suggesting support for range-separated or screened Coulomb interactions. The compilation of this file's content is conditional on `libint2` being available (`!defined(NO_LIBINT)`).

## Key Components

*   **`libint_init(int *ammax, bool *has_onebody, bool *has_gradient)`**:
    *   **Purpose**: Initializes static data for the `libint2` library. It also queries and returns capabilities of the linked `libint2` library, such as the maximum supported angular momentum (`ammax`), whether one-body integrals are supported (`has_onebody`), and whether integral gradients are available (`has_gradient`).
    *   This function should be called once before any other `libint2` related functions.

*   **`libint_2center(amA, contrdepthA, A, alphaA, cA, amC, contrdepthC, C, alphaC, cC, rcut, eriAC)`**:
    *   **Purpose**: Computes 2-center electron repulsion integrals, often denoted as (A|C) or <&phi;<sub>A</sub>| 1/r<sub>12</sub> |&phi;<sub>C</sub>>. These are less common than 4-center ERIs but can be used in specific contexts like density fitting or specialized integral evaluation schemes.
    *   **Parameters**:
        *   `amA`, `amC`: Angular momenta for basis functions on centers A and C.
        *   `contrdepthA`, `contrdepthC`: Contraction depths.
        *   `A[]`, `C[]`: Coordinates of centers A and C.
        *   `alphaA[]`, `alphaC[]`: Primitive exponents.
        *   `cA[]`, `cC[]`: Contraction coefficients.
        *   `rcut` (Input, `double`): Cutoff radius for range-separated or screened Coulomb interaction.
        *   `eriAC[]` (Output, `double[]`): Array to store computed 2-center ERIs.

*   **`libint_3center(amA, contrdepthA, A, alphaA, cA, amC, contrdepthC, C, alphaC, cC, amD, contrdepthD, D, alphaD, cD, rcut, eriACD)`**:
    *   **Purpose**: Computes 3-center ERIs, e.g., (A|CD) = <&phi;<sub>A</sub>| 1/r<sub>12</sub> |&phi;<sub>C</sub>&phi;<sub>D</sub>>, where the charge distribution on the right is a product of two GTOs centered at C and D. These are common in density fitting (RI) approximations.
    *   **Parameters**: Similar structure, adding parameters for basis function &phi;<sub>D</sub> (amD, contrdepthD, D, alphaD, cD).
        *   `eriACD[]` (Output, `double[]`): Array for 3-center ERIs.

*   **`libint_4center(amA, contrdepthA, A, alphaA, cA, amB, contrdepthB, B, alphaB, cB, amC, contrdepthC, C, alphaC, cC, amD, contrdepthD, D, alphaD, cD, rcut, eriABCD)`**:
    *   **Purpose**: Computes the standard 4-center, 2-electron ERIs (AB|CD) = <&phi;<sub>A</sub>&phi;<sub>B</sub>| 1/r<sub>12</sub> |&phi;<sub>C</sub>&phi;<sub>D</sub>>. These are the most computationally intensive integrals in many quantum chemistry methods.
    *   **Parameters**: Parameters for four basis functions (&phi;<sub>A</sub>, &phi;<sub>B</sub>, &phi;<sub>C</sub>, &phi;<sub>D</sub>).
        *   `eriABCD[]` (Output, `double[]`): Array for 4-center ERIs.

**Common Internal Workflow for ERI functions**:
1.  Initialization of `libint2` specific data structures (e.g., `Libint_eri_t`).
2.  Looping over all combinations of primitive Gaussians in the input contracted GTOs.
3.  For each primitive combination:
    *   Calculation of product Gaussian centers (P from A & B, Q from C & D) and exponents.
    *   Setup of parameters for `libint2` evaluators.
    *   Calculation of `gammapq_rc2` and `gammapq_ratio` from `rcut` to handle potentially screened interactions.
    *   Evaluation of Boys function F<sub>m</sub>(U) using `boys_function_c`, where U is adjusted by `gammapq_rc2`.
    *   The fundamental integral values (Boys function results) are scaled by factors involving `gammapq_ratio`.
4.  If total angular momentum is zero, direct summation over primitive contributions. Otherwise, invocation of `libint2_build_...eri` routines for full Cartesian tensor evaluation.
5.  Storing results in the output array.
6.  Cleanup of `libint2` resources.

## Important Variables/Constants

*   **Input Parameters (per function)**: See descriptions above. `rcut` is a key parameter for modifying the Coulomb interaction.
*   **Output Arrays**: `eriAC`, `eriACD`, `eriABCD` store the ERI tensor elements.
*   **`Libint_...eri_t` (`inteval`)**: `libint2` data structures for ERI computations.
*   `pi_2p5`: Constant &pi;<sup>2.5</sup>, typically used in ERI normalization factors.
*   `gammapq_rc2`, `gammapq_ratio`: Variables used to implement screened/range-separated Coulomb interactions based on `rcut`.

## Usage Examples

These functions are primarily intended for internal use by MOLGW, called from its Fortran core.

Conceptual Fortran call for 4-center ERIs:
```fortran
! Assuming basis function parameters (am, contrdepth, coords, exponents, coeffs) for A, B, C, D are defined
! And rcut_value is defined.
DOUBLE PRECISION, ALLOCATABLE :: ERI_VALUES(:)
INTEGER :: NUM_INTS_TOTAL

! Calculate total number of integral components based on angular momenta
! NUM_INTS_TOTAL = NINT_FUNC(AMA)*NINT_FUNC(AMB)*NINT_FUNC(AMC)*NINT_FUNC(AMD)
ALLOCATE(ERI_VALUES(NUM_INTS_TOTAL))

CALL LIBINT_4CENTER(AMA, CONTR_A, COORDS_A, EXP_A, COEFF_A, &
                    AMB, CONTR_B, COORDS_B, EXP_B, COEFF_B, &
                    AMC, CONTR_C, COORDS_C, EXP_C, COEFF_C, &
                    AMD, CONTR_D, COORDS_D, EXP_D, COEFF_D, &
                    RCUT_VALUE, ERI_VALUES)
! ERI_VALUES array now holds the computed 4-center ERIs.
```

## Dependencies and Interactions

*   **`libint2` Library**: The core computational engine for evaluating the ERIs. This file heavily relies on `libint2`'s API and data structures.
*   **`molgw.h`**: For MOLGW-specific configurations.
*   **`libint_molgw.h`**: Provides:
    *   Declaration for `boys_function_c`.
    *   The `nint(am)` inline function for Cartesian counts.
    *   The `pi_2p5` constant.
*   **`boys_function_c`**: Externally implemented (likely in Fortran as `boys_function.f90`) and used by all ERI computation routines in this file for evaluating Boys functions, which are essential for ERIs over GTOs.
*   **Standard C/C++ Libraries**: `stdlib.h`, `stdio.h`, `iostream`, `math.h`, `assert.h`.
*   **Conditional Compilation**:
    *   `#if !defined(NO_LIBINT)`: The entire file's content is conditional on this macro, allowing MOLGW to be compiled without `libint2` if desired.
    *   `#if !defined(LIBINT2_CONTRACTED_INTS)`: Asserts that contraction depth is 1 if `libint2` is not compiled to handle contracted integrals directly (implying the wrapper would only pass primitive Gaussians in that scenario, though the current code structure seems to always pass contracted info).

This file is central to MOLGW's ability to calculate two-electron integrals, a cornerstone of many quantum chemical methods.
```
