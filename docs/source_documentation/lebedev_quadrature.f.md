# `lebedev_quadrature.f`

## Overview

This Fortran 77 source file provides a collection of subroutines for generating Lebedev quadrature grids on the surface of a sphere. These grids are used for numerical integration of functions over spherical surfaces and are constructed to be exact for polynomials up to a certain algebraic order. The code is a translation from an original C implementation by Dr. Dmitri N. Laikov, translated by Dr. Christoph van Wuellen. The file contains routines for generating grids of various sizes, from 6 to 5810 points.

The comments within the code strongly request users to cite the original publication:
*V.I. Lebedev, and D.N. Laikov, "A quadrature formula for the sphere of the 131st algebraic order of accuracy," Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.*

## Key Components

*   **`gen_oh(code, num, x, y, z, w, a, b, v)`**:
    This is a worker subroutine that generates a set of symmetrically equivalent points on the sphere under octahedral (Oh) symmetry.
    *   `code` (Integer, Input): An integer from 1 to 6 that specifies the type of generating point and the symmetry operations to apply.
        *   `code=1`: Generates 6 points of type (a,0,0) and its permutations (e.g., (1,0,0), (-1,0,0), (0,1,0), ...). `a` is set to 1.0 internally.
        *   `code=2`: Generates 12 points of type (0,a,a) and its permutations. `a` is set to `1/sqrt(2)` internally.
        *   `code=3`: Generates 8 points of type (a,a,a) and its permutations. `a` is set to `1/sqrt(3)` internally.
        *   `code=4`: Generates 24 points of type (a,a,b) and its permutations. `b` is calculated as `sqrt(1-2a^2)`. Requires `a` as input.
        *   `code=5`: Generates 24 points of type (a,b,0) and its permutations. `b` is calculated as `sqrt(1-a^2)`. Requires `a` as input.
        *   `code=6`: Generates 48 points of type (a,b,c) and its permutations. `c` is calculated as `sqrt(1-a^2-b^2)`. Requires `a` and `b` as input.
    *   `num` (Integer, Input/Output): On input, it's the starting index in the arrays `x,y,z,w` where new points will be stored. On output, it's incremented by the number of points generated in the call.
    *   `x(*), y(*), z(*)` (Double Precision Arrays, Output): Arrays to store the Cartesian coordinates of the generated grid points.
    *   `w(*)` (Double Precision Array, Output): Array to store the (identical) weight `v` for the set of generated points.
    *   `a, b` (Double Precision, Input): Parameters defining the generating point coordinates (usage depends on `code`).
    *   `v` (Double Precision, Input): The weight associated with the set of points being generated.

*   **`LDxxxx(X,Y,Z,W,N)`**:
    This is a family of subroutines where `xxxx` represents the number of points in the specific Lebedev grid (e.g., `LD0006`, `LD0014`, ..., `LD5810`). Each such subroutine generates a complete Lebedev grid of a predefined order and number of points.
    *   `X(xxxx), Y(xxxx), Z(xxxx)` (Double Precision Arrays, Output): Store the Cartesian coordinates (x, y, z) of the `xxxx` grid points.
    *   `W(xxxx)` (Double Precision Array, Output): Stores the corresponding weights for each grid point.
    *   `N` (Integer, Output): Set to the total number of points generated for the specific grid (i.e., `xxxx`).
    Each `LDxxxx` subroutine internally calls `gen_oh` one or more times with specific, hardcoded values for `code`, `a`, `b`, and `v` that define that particular Lebedev grid.

## Important Variables/Constants

*   Within `gen_oh`:
    *   `a, b, v`: Input parameters defining the generator point and its weight.
    *   `c`: Calculated internally based on `a` and `b` for normalization (points lie on the unit sphere).
*   Within `LDxxxx` subroutines:
    *   `A, B, V`: Hardcoded double precision parameters passed to `gen_oh` to generate specific sets of points for that grid. These are the defining constants for each Lebedev rule.
    *   The number of points (e.g., 6 for `LD0006`, 5810 for `LD5810`) is implicitly a constant for each specific subroutine, defining array sizes and loop bounds.

## Usage Examples

To generate a 110-point Lebedev grid:
```fortran
PROGRAM USE_LEBEDEV_GRID
  IMPLICIT NONE
  INTEGER, PARAMETER :: MAX_POINTS = 110 ! Max points for LD0110
  DOUBLE PRECISION X_PTS(MAX_POINTS), Y_PTS(MAX_POINTS)
  DOUBLE PRECISION Z_PTS(MAX_POINTS), W_PTS(MAX_POINTS)
  INTEGER :: N_ACTUAL_POINTS
  INTEGER :: I

  CALL LD0110(X_PTS, Y_PTS, Z_PTS, W_PTS, N_ACTUAL_POINTS)

  PRINT *, 'Number of grid points generated: ', N_ACTUAL_POINTS
  PRINT *, 'Points (x, y, z) and weights (w):'
  DO I = 1, N_ACTUAL_POINTS
    PRINT '(4F12.8)', X_PTS(I), Y_PTS(I), Z_PTS(I), W_PTS(I)
  ENDDO

END PROGRAM USE_LEBEDEV_GRID
```
The user must ensure that the arrays passed to the `LDxxxx` subroutines are large enough to hold the coordinates and weights for that specific grid. The output `N` will indicate the actual number of points generated.

## Dependencies and Interactions

*   **Internal Dependency**: All `LDxxxx` subroutines depend on the `gen_oh` subroutine.
*   **Fortran Intrinsics**: Uses standard Fortran functions like `SQRT` (for `DSQRT` with double precision arguments).
*   **Self-Contained**: The file is self-contained for the generation of the grid points and weights. It does not rely on other modules or external files from the MOLGW project for its core calculations.
*   **Calling Code**: Expected to be called by other parts of a larger program (like MOLGW) that require spherical quadrature grids, for example, in Density Functional Theory (DFT) calculations for integrating angular parts of matrix elements or electron densities.
*   The `implicit logical(a-z)` statement is present but does not affect the behavior of the explicitly typed variables used in the core logic.
```
