# `m_inputparam.f90`

## Overview

The `m_inputparam` Fortran module is the central hub for managing all user-defined input parameters within the MOLGW software package. It is responsible for reading the main input file (typically using a Fortran `NAMELIST`), setting default values for parameters not specified by the user, performing consistency checks on these parameters, and then making them available to other modules throughout the program. It also initializes derived type objects like `calculation_type` and `excitation_type` based on the parsed input strings, effectively determining the nature and specifics of the calculation to be performed.

## Key Components

*   **Module Variables (Input Parameters)**:
    *   A large collection of variables corresponding to keywords that can be set in the MOLGW input file. These are declared with the `protected` attribute.
    *   Their default values are assigned directly in their declaration or through an included file (`input_variables.f90`).
    *   Examples include: `scf` (SCF method), `postscf` (post-SCF method), `basis` (basis set name), `auxil_basis` (auxiliary basis name), `alpha_hybrid` (fraction of HF exchange), `nspin` (number of spin channels), `charge`, `xyz_file`, various convergence thresholds, grid quality settings, print control flags, etc.
    *   These variables are made accessible to other modules via `USE m_inputparam`.

*   **`calculation_type` (Derived Type)**:
    *   Stores a set of logical flags and integer codes that categorize the overall calculation type (e.g., `is_dft`, `is_gw`, `is_ci`, `selfenergy_approx`). This type is populated by `init_calculation_type`.

*   **`excitation_type` (Derived Type)**:
    *   Stores parameters defining an external perturbation if a real-time simulation is requested (e.g., name, form, field strength `kappa`, width, start time `time0`, direction `dir`). Populated by `init_excitation_type`.

*   **`dft_xc(:)` (Array of `dft_xc_info` type from `m_libxc_tools`)**:
    *   Stores detailed information about the DFT exchange-correlation functionals to be used, including their Libxc identifiers, mixing coefficients, and parameters (like &omega; for range-separated hybrids). Populated by `init_dft_type`.

*   **`read_inputfile_namelist()` (Subroutine)**:
    *   **Purpose**: This is the primary routine for reading and processing user inputs.
    *   **Workflow**:
        1.  Includes `input_variables.f90`, which defines the `NAMELIST /molgw/` and sets default values for all input parameters.
        2.  Determines the input file name (either from a command-line argument or defaults to standard input, though stdin is deprecated).
        3.  Reads the `&molgw ... &end` namelist from the specified input file, overriding default parameter values with any user-provided settings.
        4.  Includes `echo_input_variables.f90` to print the final values of all input parameters (after defaults and user input are processed).
        5.  Performs post-processing on input strings (e.g., capitalization, converting 'yes'/'no' to logical flags).
        6.  Calls `init_excitation_type()` to set up `excit_type`.
        7.  Calls `setup_nuclei()` to read atomic coordinates and initialize atomic/ECP information.
        8.  Calls `init_calculation_type()` to parse `scf` and `postscf` strings and set flags in `calc_type`.
        9.  Calls `init_dft_type()` if a DFT calculation is specified, to parse functional names and set up `dft_xc`.
        10. Performs extensive consistency checks on the validity and compatibility of various input parameters.
        11. Calls `summary_input()` to print a concise summary of key settings.
        12. If `print_yaml_` is true, it opens a YAML output file (`yaml_output` filename) and writes input parameters and system information to it using `echo_input_variables_yaml.f90`.

*   **`init_calculation_type(scf, postscf)` (Subroutine)**:
    *   Parses the `scf` and `postscf` character strings to determine the type of calculation and sets corresponding flags in the global `calc_type` object (e.g., `is_dft`, `is_gw`, `selfenergy_approx`, `is_bse`).

*   **`init_excitation_type()` (Subroutine)**:
    *   Parses `excit_name` and related input variables (`excit_kappa`, `excit_width`, etc.) to configure the `excit_type` object, defining external perturbations for real-time simulations.

*   **`init_dft_type(key)` (Subroutine)**:
    *   Parses the DFT functional specified by `key` (e.g., 'PBE', 'B3LYP', or a Libxc identifier string like 'LIBXC:101+130').
    *   Populates the `dft_xc` array with Libxc functional IDs and associated parameters (e.g., `alpha_hybrid`, `beta_hybrid`, `gamma_hybrid` for hybrid/range-separated functionals). It can also retrieve these parameters directly from Libxc for certain functionals.

*   **`summary_input()` (Subroutine)**:
    *   Prints a formatted summary of the most important input parameters and basic information about the molecular system (number of atoms, electrons, symmetry) to the standard output.

*   **`setup_nuclei(...)` (Subroutine)**:
    *   Reads atomic coordinates, either from the input file directly or from an external XYZ file specified by `xyz_file`.
    *   Assigns basis set names (orbital, auxiliary, ECP) to each atom.
    *   Calls `m_atoms::init_atoms` to store atomic information and `m_ecp::init_ecp` to load ECP data if specified.
    *   Calculates `zvalence` based on `zatom` and ECP core electron counts.

*   **`interpret_quality(quality)` (Function)**:
    *   Converts string keywords for quality levels (e.g., "LOW", "MEDIUM", "HIGH") into corresponding integer parameter values defined in `m_definitions`.

*   **`standardize_basis_name(basis_name_in)` (Function)**:
    *   Normalizes basis set name strings by replacing common problematic characters (e.g., '*' with 's', '+' with 'p') to ensure consistent file name lookups.

## Important Variables/Constants

*   The module itself primarily *defines* and *manages* input parameters rather than defining constants for general physical use (those are in `m_definitions`).
*   **`calc_type` (Type `calculation_type`, Protected)**: A global object holding flags that define the overall nature of the requested calculation.
*   **`dft_xc(:)` (Array of `dft_xc_info`, Protected)**: Global array storing detailed Libxc functional definitions.
*   **`excit_type` (Type `excitation_type`, Protected)**: Global object storing parameters for external field excitations.
*   **Namelist Variables**: All variables within the `NAMELIST /molgw/` (included from `input_variables.f90`) are effectively global parameters for the run, accessible via `USE m_inputparam`.
*   `unit_yaml` (Integer, Protected): Fortran unit for the YAML output file.

## Usage Examples

This module is central to any MOLGW execution. The user prepares an input file (e.g., `molgw.in`) like:
```fortran
&molgw
  xyz_file = "h2o.xyz"
  basis    = "cc-pVDZ"
  scf      = "PBE"
  postscf  = "GW"
  nspin    = 1
  charge   = 0.0
  eta      = 0.005
  ! ... other parameters ...
&end
```
MOLGW is then run, typically as `molgw molgw.in > molgw.out`. The `read_inputfile_namelist` subroutine parses `molgw.in`.

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_mpi`, `m_warning`, `m_string_tools`.
*   **System Setup**: `m_atoms` (for `init_atoms`), `m_elements` (for element data used in `setup_nuclei`), `m_ecp` (for `init_ecp`).
*   **Libxc**: `m_libxc_tools` (for `dft_xc_info` type and `init_libxc_info`). Conditionally uses Libxc via `#include <xc_funcs.h>`.
*   **Included Auto-generated Files**:
    *   `input_variable_declaration.f90`: Contains declarations and default values for all namelist variables.
    *   `input_variables.f90`: Contains the `NAMELIST /molgw/` definition.
    *   `echo_input_variables.f90`: Contains `WRITE` statements to print all input variables.
    *   `echo_input_variables_yaml.f90`: Contains `WRITE` statements for YAML output.
    *   `basis_path.f90`: Contains the default basis path.
*   This module is foundational; its output (the processed input parameters and derived settings like `calc_type`) dictates the behavior of almost all subsequent modules in a MOLGW run.
```
