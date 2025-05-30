# `boys_function.f90`

## Overview

This file provides a Fortran subroutine for the evaluation of the Boys function, F_m(T). The Boys function is crucial in quantum chemistry for calculating the values of two-electron repulsion integrals (Coulomb integrals) over Gaussian basis functions. This implementation is noted as a Fortran version of an example from the LIBINT library and is part of the MOLGW software package.

## Key Components

*   **`boys_function_c(fnt, nn, tt)`**: This is the primary subroutine in the file.
    *   `nn`: An input integer representing the maximum order `m` of the Boys function to be calculated.
    *   `tt`: An input double-precision floating-point number representing the argument `T` of the Boys function.
    *   `fnt`: An output array of double-precision numbers, where `fnt(m)` will store the value of F_m(T) for `m` from 0 to `nn`.
    The subroutine uses different algorithms based on the value of `tt`:
    *   For `tt > 20.0`, it uses an upward recursion formula.
    *   For smaller `tt` ( `tt <= 20.0`), it computes F_nn(T) using an asymptotic series expansion and then uses a downward recursion formula to compute F_m(T) for `m < nn`.
    The `BIND(C)` attribute indicates that this Fortran subroutine is intended to be callable from C code.

## Important Variables/Constants

*Within `boys_function_c`:*

*   **`maxfac` (Parameter, Integer, value: 100)**: Defines the maximum limit for certain series expansions or pre-calculated factorial-like terms.
*   **`eps` (Parameter, Double, value: 1.0e-17)**: A small tolerance value used as a convergence criterion in the asymptotic series expansion for small `tt`.
*   **`kk` (Parameter, Double, value: 0.8862269254527579)**: Represents the constant `0.5 * sqrt(pi)`.
*   **`df` (Save, Double Array, size: 2*`maxfac`)**: An array used to store pre-calculated values, likely double factorials `(2k-1)!!` or similar, to speed up computations. It's initialized once and reused in subsequent calls.
*   **`et` (Double)**: Stores the value of `EXP(-tt)`.
*   **`t2` (Double)**: Stores the value of `2.0 * tt`.

## Usage Examples

While no direct runnable example is provided within the file, the subroutine would typically be called from other parts of a quantum chemistry program where Coulomb integrals are computed. For instance, if one needed the Boys functions up to order 2 for an argument `T = 5.0`:

```fortran
! Assuming this code is within a larger Fortran program/module
! that has access to boys_function_c or its interface.

program test_boys
  implicit none
  integer, parameter :: nn_max = 2
  real(8) :: t_val = 5.0_8
  real(8) :: f_values(0:nn_max)
  integer :: i

  ! Call the Boys function subroutine
  call boys_function_c(f_values, nn_max, t_val)

  ! Print the results
  do i = 0, nn_max
    print *, "F_", i, "(", t_val, ") = ", f_values(i)
  enddo

end program test_boys
```
*(Note: The above is a hypothetical usage example. The actual integration into MOLGW would involve more complex data structures and control flow.)*

## Dependencies and Interactions

*   **`ISO_C_BINDING`**: This intrinsic Fortran module is used to ensure interoperability with C code, as indicated by the `BIND(C)` attribute and the use of C-compatible types like `C_INT` and `C_DOUBLE`.
*   **`molgw.h`**: This is a C-style header file included via `#include "molgw.h"`. It likely contains common definitions, type kinds (though `C_INT`, `C_DOUBLE` are usually from `ISO_C_BINDING`), or macros used across the MOLGW project.
*   **Intrinsic Fortran Functions**: The subroutine uses standard Fortran intrinsic functions such as `EXP` (exponential), `SQRT` (square root), `ERF` (error function), and `ABS` (absolute value).
*   **Coulomb Integral Calculations**: The primary interaction of this subroutine is with the parts of MOLGW responsible for calculating two-electron Coulomb integrals. The values computed by `boys_function_c` are essential inputs for these integral calculations.

```
