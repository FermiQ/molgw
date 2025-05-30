# `libint_molgw.h`

## Overview

This C/C++ header file acts as a bridge or wrapper for integrating the `libint` library (a specialized library for computing molecular integrals) with the MOLGW software package. It declares specific functions and constants that are utilized by MOLGW's C/C++ components when interfacing with `libint` or performing related calculations, such as evaluating Boys functions.

## Key Components

*   **`boys_function_c(double* fnt, int nn, double tt)` (extern "C")**:
    *   This is an external function declaration, indicating that `boys_function_c` is likely implemented in another language (Fortran, in this case, as seen in `boys_function.f90`) and compiled separately. It is made available to C/C++ code through this declaration.
    *   **Purpose**: Calculates Boys function values, F_m(T), which are essential for the evaluation of Gaussian electron repulsion integrals.
    *   `fnt` (Output, `double*`): Pointer to an array where the computed Boys function values will be stored.
    *   `nn` (Input, `int`): The maximum order `m` of the Boys function to be calculated.
    *   `tt` (Input, `double`): The argument `T` of the Boys function.

*   **`nint(int am)` (inline function)**:
    *   **Purpose**: Calculates the number of Cartesian Gaussian basis functions for a given angular momentum `am`. For example, for `am=0` (s-shell), it returns 1; for `am=1` (p-shell), it returns 3 (px, py, pz); for `am=2` (d-shell), it returns 6 (dxx, dxy, dxz, dyy, dyz, dzz).
    *   The formula used is `(am+1)*(am+2)/2`.

## Important Variables/Constants

*   **`pi_2p5` (const double)**:
    *   A constant defined as &pi;<sup>2.5</sup> (i.e., `pow(M_PI, 2.5)`). `M_PI` is a standard macro representing the value of &pi;, typically found in `<cmath>` or `<math.h>`.
    *   This constant is likely used in integral evaluation formulas where this specific power of &pi; appears.
    *   The preprocessor directives `#ifndef PI_2P5 #define PI_2P5 ... #endif` ensure that this constant is defined only once, even if the header is included multiple times.

## Usage Examples

This header file is not meant to be run directly but included in C/C++ source files.

To use `nint`:
```cpp
#include "libint_molgw.h"
// ...
int angular_momentum = 2; // d-shell
int num_cartesian_functions = nint(angular_momentum);
// num_cartesian_functions will be 6
```

To use `boys_function_c` (conceptual, actual call would be from C++ code performing integral calculations):
```cpp
#include "libint_molgw.h"
// ...
// Assume tt_value and max_order are defined
double T = tt_value;
int N = max_order;
double* F_values = new double[N + 1];

boys_function_c(F_values, N, T);

// F_values array now contains F_0(T), F_1(T), ..., F_N(T)
// ...
delete[] F_values;
```

## Dependencies and Interactions

*   **`boys_function.f90`**: The `boys_function_c` declaration relies on the external Fortran implementation of this function, which must be compiled and linked with the C/C++ object files that use it.
*   **`libint` library**: This header is a component of the interface to the `libint` library. C++ source files that include `libint_molgw.h` are likely those that make direct calls to `libint` routines for computing molecular integrals.
*   **Standard C/C++ Libraries**:
    *   Requires `<cmath>` (or `<math.h>`) for the `pow()` function and the `M_PI` constant used in the definition of `pi_2p5`.
    *   The `using namespace std;` statement is present, indicating a C++ context and usage of the `std` namespace.
*   **MOLGW C/C++ Source Files**: This header is included by C/C++ source files within MOLGW that are responsible for tasks like integral computation or basis set management.

```
