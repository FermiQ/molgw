# `libint_onebody.cpp`

## Overview

This C++ source file serves as a wrapper to the `libint2` library for computing one-body molecular integrals within the MOLGW software package. It provides C-callable functions (compatible with Fortran) to calculate overlap, kinetic energy, and nuclear attraction integrals between contracted Gaussian basis functions. The functionality of this file is conditional on the `libint2` library being available and configured to support these one-body integrals.

## Key Components

The file defines three main functions, each responsible for calculating a specific type of one-body integral:

*   **`libint_overlap(amA, contrdepthA, A, alphaA, cA, amB, contrdepthB, B, alphaB, cB, overlapAB)`**:
    *   **Purpose**: Computes the overlap integrals S<sub>&mu;&nu;</sub> = <&phi;<sub>&mu;</sub> | &phi;<sub>&nu;</sub>> between two contracted Gaussian-type orbitals (GTOs), &phi;<sub>&mu;</sub> centered at A and &phi;<sub>&nu;</sub> centered at B.
    *   **Parameters**:
        *   `amA`, `amB` (Input, `int`): Angular momenta of the basis functions on centers A and B.
        *   `contrdepthA`, `contrdepthB` (Input, `int`): Number of primitive Gaussians in the contraction for functions A and B.
        *   `A[]`, `B[]` (Input, `double[3]`): Cartesian coordinates of the centers for basis functions A and B.
        *   `alphaA[]`, `alphaB[]` (Input, `double[]`): Arrays of exponents for the primitive Gaussians.
        *   `cA[]`, `cB[]` (Input, `double[]`): Arrays of contraction coefficients for the primitive Gaussians.
        *   `overlapAB[]` (Output, `double[]`): Array to store the computed overlap integral values. The size depends on `amA` and `amB`.

*   **`libint_kinetic(amA, contrdepthA, A, alphaA, cA, amB, contrdepthB, B, alphaB, cB, kineticAB)`**:
    *   **Purpose**: Computes the kinetic energy integrals T<sub>&mu;&nu;</sub> = <&phi;<sub>&mu;</sub> | -1/2 &nabla;<sup>2</sup> | &phi;<sub>&nu;</sub>> between two contracted GTOs.
    *   **Parameters**: Similar to `libint_overlap`, with `kineticAB[]` as the output array for kinetic energy integrals.

*   **`libint_elecpot(amA, contrdepthA, A, alphaA, cA, amB, contrdepthB, B, alphaB, cB, C, elecpotAB)`**:
    *   **Purpose**: Computes the nuclear attraction (electron-potential) integrals V<sub>&mu;&nu;</sub> = <&phi;<sub>&mu;</sub> | 1/|**r** - **R**<sub>C</sub>| | &phi;<sub>&nu;</sub>> between two contracted GTOs, due to a point charge at center C.
    *   **Parameters**: Similar to `libint_overlap`, with the addition of:
        *   `C[]` (Input, `double[3]`): Cartesian coordinates of the center C (e.g., a nucleus).
        *   `elecpotAB[]` (Output, `double[]`): Array to store the computed nuclear attraction integral values.

**Common Internal Workflow for each function**:
1.  Initializes `libint2` specific structures (`Libint_overlap_t`, `Libint_kinetic_t`, `Libint_elecpot_t`).
2.  Loops over all pairs of primitive Gaussians from the two input contracted GTOs.
3.  For each primitive pair, calculates necessary geometric and exponential parameters.
4.  Sets up the `libint2` evaluation structure with these parameters. For `libint_elecpot`, this step includes calling `boys_function_c` to compute required Boys function values.
5.  If the total angular momentum is zero (S-S type interaction), the integral is computed by summing contributions from primitives.
6.  Otherwise, calls the appropriate `libint2_build_...` routine to compute all Cartesian components of the integrals for the given angular momenta.
7.  Stores the final summed/contracted integral values in the output array.
8.  Cleans up `libint2` structures.

## Important Variables/Constants

*   **Input Parameters (per function)**:
    *   `amA`, `amB`: Angular momenta of the two basis functions.
    *   `contrdepthA`, `contrdepthB`: Contraction depths (number of primitives).
    *   `A`, `B`, `C`: Cartesian coordinates of atomic centers (A, B) and potential center (C).
    *   `alphaA`, `alphaB`: Arrays of primitive Gaussian exponents.
    *   `cA`, `cB`: Arrays of contraction coefficients.
*   **Output Arrays**: `overlapAB`, `kineticAB`, `elecpotAB` store the resulting integral matrices.
*   **`Libint_..._t` (`inteval`)**: `libint2` specific data structures used to manage integral computation parameters and results for each primitive pair.
*   **`LIBINT2_MAX_AM_...`**: Constants from `libint2` defining the maximum supported angular momentum for each type of integral.
*   `M_PI`: The mathematical constant &pi;, used in Gaussian normalization and other factors.

## Usage Examples

These functions are designed to be called from Fortran, so their C `extern "C"` linkage is crucial. A typical usage pattern from the Fortran side would involve:
1.  Defining the basis set parameters (atomic coordinates, exponents, contraction coefficients, angular momenta).
2.  Allocating Fortran arrays to receive the integral values.
3.  Calling these C functions via Fortran's interoperability features.

Conceptual Fortran call:
```fortran
! Assuming necessary variables (AM_A, NPRIM_A, COORDS_A, etc.) are defined
DOUBLE PRECISION, ALLOCATABLE :: OVERLAP_VALUES(:)
INTEGER :: NUM_CART_A, NUM_CART_B

! NUM_CART_A = NINT_FUNC(AM_A) ! Using a function similar to nint() in libint_molgw.h
! NUM_CART_B = NINT_FUNC(AM_B)
ALLOCATE(OVERLAP_VALUES(NUM_CART_A * NUM_CART_B))

CALL LIBINT_OVERLAP(AM_A, NPRIM_A, COORDS_A, EXP_A, COEFF_A, &
                    AM_B, NPRIM_B, COORDS_B, EXP_B, COEFF_B, &
                    OVERLAP_VALUES)
! OVERLAP_VALUES now contains the computed overlap integrals
```

## Dependencies and Interactions

*   **`libint2` Library**: This file is fundamentally a wrapper around the `libint2` library. It makes extensive use of `libint2`'s data structures and API functions (e.g., `libint2::malloc`, `libint2_init_overlap`, `libint2_build_overlap`, etc.).
*   **`molgw.h`**: Included for general MOLGW project configurations or definitions.
*   **`libint_molgw.h`**: This header provides:
    *   The declaration for `boys_function_c` (used by `libint_elecpot`).
    *   The inline function `nint(am)` for calculating the number of Cartesian functions.
*   **`boys_function_c`**: The `libint_elecpot` function has an external dependency on `boys_function_c` (likely the Fortran implementation in `boys_function.f90`) for evaluating Boys functions.
*   **Standard C/C++ Libraries**:
    *   `stdlib.h`, `stdio.h`, `iostream`, `math.h` (for `sqrt`, `exp`, `M_PI`), `assert.h`.
*   **Conditional Compilation**:
    *   `#if !defined(NO_LIBINT)`: The entire content of the file is guarded by this preprocessor directive. If `NO_LIBINT` is defined, the file effectively becomes empty, meaning MOLGW would be compiled without `libint2` support for these one-body integrals.
    *   `#if defined(LIBINT2_SUPPORT_ONEBODY)`: The three integral functions are only compiled if `libint2` is configured to support one-body integrals.
    *   `#if !defined(LIBINT2_CONTRACTED_INTS)`: Includes an assertion to ensure contraction depth is 1 if `libint2` is not compiled for contracted integrals, implying that the wrapper expects `libint2` to handle contractions internally if this macro is not defined.

This file acts as a crucial interface layer, translating MOLGW's requests for one-body integrals into calls understandable by the `libint2` library.
```
