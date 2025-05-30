# `m_elements.f90`

## Overview

The `m_elements` Fortran module in MOLGW serves as a centralized repository for data and functions related to chemical elements. It provides a periodic table lookup for element symbols and atomic numbers, functions to retrieve element-specific properties like covalent radii and the number of core electrons (relevant for Effective Core Potential, ECP, setups), and a routine to generate parameters for simple atomic density representations.

## Key Components

*   **`element_list(nelement_max)` (Character(len=2), Parameter)**:
    *   An array storing the 2-character chemical symbols for elements up to atomic number `nelement_max` (currently 103, Lawrencium).

*   **`element_core(zval, zatom)` (Function)**:
    *   **Purpose**: Determines the number of core electrons for a given element, typically for deciding which electrons can be replaced by an ECP.
    *   If `zval` (valence charge) is significantly different from `zatom` (nuclear charge), it assumes an ECP is in use for all core electrons and returns 0.
    *   Otherwise, it returns a predefined number of core electrons based on `zval` (e.g., 0 for Z ≤ 4, 1 for Z ≤ 12, etc.).

*   **`element_covalent_radius(zatom)` (Function)**:
    *   **Purpose**: Returns the covalent radius of an element in Bohr units.
    *   Input `zatom` is the atomic number.
    *   Uses a hardcoded list of covalent radii (data sourced from Cambridge Structural Database via Wikipedia, originally in picometers and converted to Bohr). Provides a default value if the element is not in the list.

*   **`element_number(element_name)` (Function)**:
    *   **Purpose**: Converts a chemical symbol (e.g., "H", "He") into its corresponding atomic number (integer).
    *   Case-insensitive comparison is not explicitly mentioned but often handled by `ADJUSTL` before comparison.

*   **`element_name(zatom)` (Function)**:
    *   **Purpose**: Converts an atomic number (`zatom`, can be integer or real) into its 2-character chemical symbol.
    *   Handles Z=0 by returning "X ". For antinuclei (negative `zatom`), it uses the absolute value for lookup.

*   **`element_name_long(zatom)` (Function)**:
    *   **Purpose**: Similar to `element_name` but returns a longer character string (length 8).
    *   If `zatom` is an integer, it returns the symbol. If `zatom` is real and significantly different from its nearest integer (e.g., for non-integer nuclear charges in some models), it formats the real number as a string.

*   **`element_atomicdensity(zval, zatom, coeff, alpha)` (Subroutine)**:
    *   **Purpose**: Provides parameters for a simple atomic density representation using a sum of four Gaussian functions.
    *   `zval`: Valence charge.
    *   `zatom`: Nuclear charge.
    *   `coeff(4)` (Output): Coefficients for the four Gaussians.
    *   `alpha(4)` (Output): Exponents for the four Gaussians.
    *   Uses hardcoded parameters for H, C, N, O, and Au (with special handling for Au ECPs). For other elements, it uses a generic two-Gaussian fit based on `zval`. The coefficients are normalized to sum to `zval`.

## Important Variables/Constants

*   **`nelement_max` (Integer, Parameter)**: Defines the size of the `element_list` array and thus the highest atomic number directly supported by symbol lookup (103).
*   Physical constants like `bohr_A` (conversion factor from Bohr to Angstrom) are used from `m_definitions` for unit conversions in `element_covalent_radius`.

## Usage Examples

```fortran
USE m_elements
USE m_definitions ! For dp kind and bohr_A if used directly
IMPLICIT NONE

INTEGER :: z_carbon, core_e_carbon
REAL(DP) :: radius_carbon_bohr, radius_carbon_pm
CHARACTER(LEN=2) :: symbol_carbon
CHARACTER(LEN=8) :: long_symbol_carbon
REAL(DP) :: carbon_density_coeffs(4), carbon_density_exps(4)

! Get atomic number from symbol
z_carbon = element_number('C') ! Returns 6

! Get symbol from atomic number
symbol_carbon = element_name(z_carbon) ! Returns "C "
long_symbol_carbon = element_name_long(REAL(z_carbon, dp)) ! Returns "C       "

! Get covalent radius
radius_carbon_bohr = element_covalent_radius(z_carbon)
radius_carbon_pm = radius_carbon_bohr * bohr_A * 100.0_dp ! Convert Bohr to pm

! Get number of core electrons (assuming Z_valence = Z_atom for this example)
core_e_carbon = element_core(REAL(z_carbon, dp), REAL(z_carbon, dp)) ! Returns 0 for C

! Get atomic density parameters
CALL element_atomicdensity(REAL(z_carbon, dp), REAL(z_carbon, dp), &
                           carbon_density_coeffs, carbon_density_exps)

WRITE(*,*) "Carbon (Z=", z_carbon, ", Symbol='", TRIM(symbol_carbon), "')"
WRITE(*,*) "Covalent Radius (Bohr): ", radius_carbon_bohr
WRITE(*,*) "Core electrons: ", core_e_carbon
WRITE(*,*) "Atomic Density Gaussians (coeffs, exponents):"
DO I = 1, 4
  IF (ABS(carbon_density_coeffs(I)) > 1e-6_dp) THEN
    WRITE(*,*) carbon_density_coeffs(I), carbon_density_exps(I)
  END IF
END DO

END
```

## Dependencies and Interactions

*   **`m_definitions`**: For the `dp` kind parameter and physical constants like `bohr_A`.
*   **`m_warning`**: Uses the `die` subroutine for error handling if an unknown element symbol is provided or an unsupported atomic number is encountered.
*   **Other MOLGW Modules**: This module provides fundamental atomic data that is widely used across MOLGW:
    *   `m_atoms`: Uses `element_name`, `element_number`, `element_covalent_radius`, `element_core` to set up and characterize the atomic system.
    *   `m_basis_set`: Uses `element_name` to construct basis set file names.
    *   `m_ecp`: Uses `element_name` to construct ECP file names and `element_number` to identify elements.
    *   Any module performing element-specific analysis or requiring element properties will likely interact with `m_elements`.
```
