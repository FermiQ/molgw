# `m_definitions.f90`

## Overview

The `m_definitions` module is a foundational part of the MOLGW software package. Its primary role is to define and provide access to a wide range of global parameters, constants, and type kinds that are used throughout the entire application. This includes numerical precision kinds, standard I/O units, physical and mathematical constants, numerical thresholds, and global flags related to library capabilities (like Libcint). It is intended to be `USE`d by most, if not all, other Fortran modules and subroutines within MOLGW to ensure consistency.

## Key Components

*   **Module `m_definitions`**: The sole module in this file.
    *   **Type Kind Definitions**:
        *   `dp`: Double precision real kind (`KIND(0.0d0)`).
        *   `sp`: Single precision real kind (`KIND(0.0)`).
        *   `int8`: Integer kind parameter (set to 8, nominally for 8-byte integers, though standard Fortran kind selection is usually preferred e.g., `SELECTED_INT_KIND(R)` for a certain range R).
    *   **Standard I/O Units**:
        *   `stdout`: Integer variable for standard output, initialized to `OUTPUT_UNIT` from `ISO_FORTRAN_ENV`. Can be redirected using `set_standard_output`.
        *   `stderr`: Integer parameter for standard error, set to `ERROR_UNIT` from `ISO_FORTRAN_ENV`.
    *   **Libcint Capability Placeholders**:
        *   `MOLGW_LMAX` (Integer, Protected): Stores the maximum angular momentum supported by the linked Libcint library. Set by `set_molgw_lmax`.
        *   `MOLGW_has_onebody` (Logical, Protected): Flag indicating if Libcint supports one-body integrals. Set by `set_molgw_lmax`.
        *   `MOLGW_has_gradient` (Logical, Protected): Flag indicating if Libcint supports integral gradients. Set by `set_molgw_lmax`.
    *   **Physical Constants**: Defines parameters for various physical constants (e.g., `au_debye`, `Ha_eV`, `bohr_A`, `Ha_K`), sourced from NIST CODATA 2010.
    *   **Mathematical Constants**: Defines parameters for &pi;, &pi;<sup>2</sup>, common numerical values (0.0, 0.5, 1.0, etc.), and complex numbers `im`, `COMPLEX_ONE`, `COMPLEX_ZERO`.
    *   **Numerical Thresholds**: A set of predefined `real(dp)` parameters for various tolerance levels used in comparisons (e.g., `tol5 = 1.0e-5_dp`).
    *   **Quality Levels**: Integer parameters (`low`, `medium`, `high`, etc.) for defining quality of grids or integral calculations.
    *   **`debug` (Logical, Parameter)**: Set to `.TRUE.` if the `DEBUG` preprocessor macro is defined during compilation, otherwise `.FALSE.`.
    *   `output_name` (Character(len=140)): String variable, likely for default output filenames in NOFT calculations.

*   **`set_molgw_lmax(lmax, has_onebody, has_gradient)` (Subroutine)**:
    *   **Purpose**: Sets the global module variables `MOLGW_LMAX`, `MOLGW_has_onebody`, and `MOLGW_has_gradient`. This is typically called after initializing Libcint to reflect its capabilities.
    *   **Parameters**:
        *   `lmax` (Input, Integer): Maximum angular momentum.
        *   `has_onebody` (Input, Logical(C_BOOL)): True if one-body integrals are available.
        *   `has_gradient` (Input, Logical(C_BOOL)): True if integral gradients are available.

*   **`set_standard_output(unit_stdout)` (Subroutine)**:
    *   **Purpose**: Allows MOLGW's standard output to be redirected from the default `OUTPUT_UNIT` to a different Fortran unit number specified by `unit_stdout`.
    *   If `unit_stdout` is different from the current `OUTPUT_UNIT`, the default `OUTPUT_UNIT` is closed, and the module's `stdout` variable is reassigned and opened.

## Important Variables/Constants

*   **`dp` (Integer, Parameter)**: The most critical constant, defining the kind for double-precision real numbers. This ensures consistent precision across the codebase.
*   **`stdout` (Integer, Protected)**: The unit number for standard output. Its re-definable nature allows output redirection.
*   **`MOLGW_LMAX` (Integer, Protected)**: Defines the maximum angular momentum the program (via Libcint) can handle, affecting array dimensions and loop bounds in basis set and integral routines.
*   **Physical and Mathematical Constants**: Provide standardized values for common calculations.
*   **`debug` (Logical, Parameter)**: Enables or disables conditional debug code blocks.
*   **Numerical Tolerances**: Standardize thresholds for convergence checks, comparisons, etc.

## Usage Examples

This module is implicitly used in almost every other Fortran file in MOLGW via a `USE m_definitions` statement.

```fortran
MODULE my_calculation_module
  USE m_definitions ! Essential for type kinds and constants
  IMPLICIT NONE

  REAL(dp) :: my_energy_variable
  REAL(dp) :: distance_bohr, distance_angstrom
  COMPLEX(dp) :: complex_value

  my_energy_variable = -76.4_dp * Ha_eV ! Convert Hartrees to eV
  distance_bohr = 5.0_dp
  distance_angstrom = distance_bohr * bohr_A

  IF (ABS(my_energy_variable) < tol8) THEN
    WRITE(stdout, *) "Energy is close to zero."
  END IF

  complex_value = (1.0_dp, 0.0_dp) / im ! Example complex arithmetic

END MODULE my_calculation_module
```
The subroutines `set_molgw_lmax` and `set_standard_output` are typically called during the initial setup phases of the program.

## Dependencies and Interactions

*   **`ISO_FORTRAN_ENV` (Intrinsic Module)**: Used for `OUTPUT_UNIT` and `ERROR_UNIT`.
*   **`ISO_C_BINDING` (Intrinsic Module)**: Used for C interoperability type kinds. Although not all imported kinds are directly exposed as public parameters in this file's API, they are made available for modules that `USE m_definitions`.
*   **`OMP_LIB` (Intrinsic Module)**: Conditionally included if `_OPENMP` is defined, providing OpenMP runtime library routines. These are not directly exposed by `m_definitions` but become available.
*   **`molgw.h` (Header File)**: Included via the C preprocessor. This header likely defines the `MOLGW_LMAX` preprocessor macro (which might be different from the Fortran variable `MOLGW_LMAX` set at runtime) and the `DEBUG` macro.
*   **Compiler/Preprocessor**: The `debug` constant's value depends on the `DEBUG` preprocessor macro being defined at compile time.
*   **Libcint**: The values set by `set_molgw_lmax` are obtained from the Libcint initialization, making this module a conduit for Libcint's capabilities to the rest of MOLGW.
*   **All other MOLGW modules**: Virtually all other parts of MOLGW will `USE m_definitions` to access these fundamental parameters and type kinds, ensuring project-wide consistency.
```
