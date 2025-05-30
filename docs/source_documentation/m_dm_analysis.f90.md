# `m_dm_analysis.f90`

## Overview

The `m_dm_analysis` Fortran module in MOLGW provides tools for analyzing a given one-body reduced density matrix (1-RDM). Its primary purpose is to take a density matrix, typically from an external source (like a Gaussian fchk file) or a previous MOLGW calculation that saved it, and perform various analyses. These include reconstructing natural orbitals and their occupations, calculating the electron density on a real-space grid, dumping this density to a file, and computing various electronic properties like multipole moments. This module often serves as a standalone analysis utility that terminates the program after execution.

## Key Components

*   **`dm_dump(basis)`**:
    *   **Purpose**: This is the main and only public subroutine in the module. It orchestrates the reading, processing, and analysis of a density matrix.
    *   **Workflow**:
        1.  Reads the density matrix:
            *   If `read_fchk` input parameter is specified, it attempts to read the density matrix from a Gaussian fchk file via `read_gaussian_fchk`.
            *   Otherwise, it looks for a MOLGW-specific binary file named `DENSITY_MATRIX`.
        2.  Reconstructs natural orbitals and occupation numbers: Calls `get_c_matrix_from_p_matrix` to diagonalize the input density matrix, yielding `c_matrix_test` (natural orbitals) and `occupation_test` (natural orbital occupation numbers).
        3.  Initializes a DFT integration grid using `init_dft_grid` from the `m_dft_grid` module.
        4.  Calculates the electron density on this grid using `calc_density_r_batch` from `m_density_tools`, utilizing the reconstructed orbitals and occupations.
        5.  Writes the grid coordinates (implicitly, via point index) and corresponding density values to a file named `rho_grid.dat`.
        6.  Calculates and prints the total number of electrons by integrating the density over the grid (normalization check).
        7.  Optionally, if controlled by input parameters (`print_multipole_`, `print_cube_`, `print_wfn_files_`):
            *   Calculates static electric dipole and quadrupole moments using routines from `m_multipole`.
            *   Generates a `.cube` file for visualizing the electron density using `plot_cube_wfn`.
            *   Prints wavefunction information to a file using `print_wfn_file`.
        8.  The subroutine typically calls `this_is_the_end` to terminate the program after completing the analysis.

## Important Variables/Constants

*   **`basis` (Type `basis_set`, Intent(in))**: Defines the atomic orbital basis set within which the density matrix is represented and analyses are performed.
*   **`p_matrix_test(:, :, :)` (Real(dp), Allocatable)**: Internal array holding the one-body reduced density matrix that is read from a file.
*   **`c_matrix_test(:, :, :)` (Real(dp), Allocatable)**: Array storing the natural orbital coefficients obtained from diagonalizing `p_matrix_test`.
*   **`occupation_test(:, :)` (Real(dp), Allocatable)**: Array storing the occupation numbers of the natural orbitals.
*   **Input Parameters (from `m_inputparam`)**:
    *   `read_fchk` (Character string): Specifies the path to a Gaussian fchk file. If not 'NO', the density matrix is read from this file.
    *   `print_multipole_` (Logical): If true, triggers calculation and printing of multipole moments.
    *   `print_cube_` (Logical): If true, triggers generation of a `.cube` file for the electron density.
    *   `print_wfn_files_` (Logical): If true, triggers printing of wavefunction file(s).
    *   `grid_level` (Integer): Controls the quality (density) of the DFT grid used for density evaluation.
*   **File Names**:
    *   `'DENSITY_MATRIX'`: Default filename for reading a MOLGW-formatted binary density matrix.
    *   `'rho_grid.dat'`: Output file where the electron density on the grid is written.

## Usage Examples

This module is typically executed as a specific task in MOLGW, often for post-processing or analyzing a density matrix obtained from a different source or a prior complex calculation.

To use this, a MOLGW input file would be structured to:
1.  Define the molecular geometry (`xyz_file`).
2.  Specify the basis set (`basis`).
3.  Indicate the source of the density matrix:
    *   Set `read_fchk = 'path/to/your/file.fchk'` OR
    *   Ensure a `DENSITY_MATRIX` file is present in the working directory.
4.  Set flags for desired outputs, e.g.:
    *   `print_multipole = YES`
    *   `print_cube = YES`
    *   `task = 'DM_ANALYSIS'` (or similar, depending on how MOLGW's task system invokes this).

The program would then run `dm_dump`.

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_timing`, `m_warning`, `m_memory`, `m_inputparam`.
*   **ScaLAPACK**: `m_scalapack` is included, though direct parallel execution logic for `dm_dump` itself is noted as "not coded" and a warning for single-core execution is present.
*   **SCF Routines**: `m_scf` (or a similar module) for `get_c_matrix_from_p_matrix` to convert the density matrix into orbitals and occupations.
*   **Atomic/Basis Information**: `m_atoms`, `m_basis_set`, `m_gaussian`, `m_ecp`.
*   **Grid and Density**: `m_dft_grid` (for `init_dft_grid`, `w_grid`, `get_basis_functions_r_batch`), `m_density_tools` (for `calc_density_r_batch`).
*   **Property Calculations**: `m_multipole` (for `static_dipole`, `static_quadrupole`).
*   **I/O and Plotting**: `m_io` (or equivalent modules) for `read_gaussian_fchk`, `plot_cube_wfn`, `print_wfn_file`.
*   **Program Flow**: `dm_dump` often acts as an endpoint for a calculation, calling `this_is_the_end` upon completion.
*   **File System**: Relies on the presence of either the specified fchk file or a `DENSITY_MATRIX` file. Outputs `rho_grid.dat` and potentially other files.
```
