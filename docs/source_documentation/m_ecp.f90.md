# `m_ecp.f90`

## Overview

The `m_ecp` Fortran module in MOLGW is responsible for managing Effective Core Potentials (ECPs), also known as pseudopotentials. ECPs are used in electronic structure calculations to replace the chemically inert core electrons of an atom with an effective potential. This approximation reduces the number of electrons explicitly treated in the calculation, thereby lowering computational cost, especially for systems containing heavy elements. This module handles the reading of ECP data from various standard file formats and stores this information in a structured way for use by other parts of MOLGW.

## Key Components

*   **`effective_core_potential` (Type)**:
    *   A derived data type designed to store all parameters defining an ECP for a particular element. Its members accommodate different ECP formats:
        *   `ecp_format` (Integer): An identifier for the ECP format (e.g., `ECP_NWCHEM`, `ECP_PSP6`, `ECP_PSP8`, `ECP_GTH`).
        *   `ncore` (Integer): The number of core electrons replaced by the ECP.
        *   `necp` (Integer): The number of projectors or distinct terms in the ECP expansion.
        *   `lk(:)` (Integer, Allocatable): Array storing the angular momentum (e.g., s=0, p=1, d=2) for each projector. A value of -1 typically denotes the local part of the ECP.
        *   **For Gaussian-type ECPs (e.g., NWChem format)**:
            *   `nk(:)` (Integer, Allocatable): Exponents of `r` in the prefactor `r^(nk-2)`.
            *   `dk(:)` (Real(dp), Allocatable): Coefficients for each Gaussian term.
            *   `zetak(:)` (Real(dp), Allocatable): Exponents of the Gaussian functions.
        *   **For Goedecker-Teter-Hutter (GTH) type ECPs**:
            *   `gth_rpploc` (Real(dp)): Range of the local part.
            *   `gth_nloc` (Integer): Number of coefficients for the local part.
            *   `gth_cipp(:)` (Real(dp), Allocatable): Coefficients for the local part.
            *   `gth_nl` (Integer): Number of non-local projectors by angular momentum.
            *   `gth_npl(:)` (Integer, Allocatable): Number of projectors for each angular momentum.
            *   `gth_rl(:)` (Real(dp), Allocatable): Range parameters for non-local projectors.
            *   `gth_hijl(:, :)` (Real(dp), Allocatable): Projector strength coefficients.
        *   **For numerical grid-based ECPs (e.g., ABINIT PSP6/PSP8, ONCVPSP formats)**:
            *   `mmax` (Integer): Number of radial grid points.
            *   `rad(:)` (Real(dp), Allocatable): Radial grid points.
            *   `wfll(:, :)` (Real(dp), Allocatable): Values of pseudo-wavefunctions on the radial grid.
            *   `vpspll(:, :)` (Real(dp), Allocatable): Values of pseudopotential projectors on the radial grid.
            *   `ekb(:)` (Real(dp), Allocatable): Kleinman-Bylander energies or normalizations.
            *   `rhocore(:, :)` (Real(dp), Allocatable): Core density on the radial grid (for Non-Linear Core Correction, NLCC).

*   **`init_ecp(ecp_elements, ecp_path, ecp_name, ecp_level_in)`**:
    *   **Purpose**: Initializes ECPs for the elements specified in the `ecp_elements` string.
    *   It parses `ecp_elements` to identify which atomic species will use ECPs.
    *   Constructs filenames for ECP data based on `ecp_path`, element name, and `ecp_name`.
    *   Determines the ECP format from `ecp_name` (e.g., if it contains 'PSP6', 'GTH').
    *   Calls the appropriate file reading subroutine (`read_ecp_file`, `read_psp6_file`, etc.) for each element requiring an ECP.
    *   Sets up the quality (`nradial_ecp`, `nangular_ecp`) of the numerical grid used for integrating ECP matrix elements, unless analytical GTH ECPs are used.

*   **File Reading Subroutines**:
    *   **`read_ecp_file(...)`**: Parses ECP files in the NWChem format.
    *   **`read_psp6_file(...)`**: Parses ECP files in the ABINIT PSP6 format (older ABINIT norm-conserving).
    *   **`read_psp8_file(...)`**: Parses ECP files in the ABINIT PSP8 format (e.g., ONCVPSP).
    *   **`read_gth_file(...)`**: Parses ECP files in the Goedecker-Teter-Hutter (GTH) format, commonly used by CP2K.

## Important Variables/Constants

*   **`ECP_NWCHEM`, `ECP_PSP6`, `ECP_PSP8`, `ECP_GTH` (Integer, Parameters)**: Named constants representing the different supported ECP file formats.
*   **`nelement_ecp` (Integer, Protected)**: The number of distinct elements for which ECPs are loaded.
*   **`element_ecp(:)` (Integer, Protected, Allocatable)**: An array storing the atomic numbers of the elements that utilize ECPs.
*   **`ecp(:)` (Type(`effective_core_potential`), Allocatable)**: The main array storing the detailed ECP parameters for each element in `element_ecp`.
*   **`nradial_ecp`, `nangular_ecp` (Integer, Protected)**: Define the number of radial and angular points for the numerical integration grid used when evaluating matrix elements of non-local ECP projectors (except for GTH, which are often handled analytically or with specific schemes).

## Usage Examples

The `m_ecp` module is used internally by MOLGW. Its functions are not typically called directly by the user in the input file, but its behavior is controlled by input parameters.

Example MOLGW input snippet that would trigger `init_ecp`:
```
ecp_elements = "Si Ge"  ! Specifies that Silicon and Germanium will use ECPs
ecp_name     = "BFD"     ! Common part of the ECP name (e.g., BFD, CRENBL)
basis_path   = "./basis_sets/" ! Path where ECP files (and basis files) are located
```
MOLGW would then attempt to read files like `./basis_sets/Si_BFD` and `./basis_sets/Ge_BFD`.

## Dependencies and Interactions

*   **`m_definitions`**: For `dp` (double precision kind) and other fundamental type kinds.
*   **`m_string_tools`**: For string manipulation functions like `capitalize`, `orbital_momentum_name`, `orbital_momentum_number`, and `append_to_list`.
*   **`m_warning`**: For error handling (`die`) and issuing warnings.
*   **`m_elements`**: For mapping element symbols to atomic numbers (`element_number`, `element_name`).
*   **`ISO_FORTRAN_ENV`**: For `IOSTAT_END` used in file reading operations.
*   **File System**: This module heavily interacts with the file system to read ECP data from external files. The paths and names of these files are determined by user input.
*   **Hamiltonian Construction**: The ECP data loaded and processed by this module is subsequently used by Hamiltonian construction modules (e.g., `m_hamiltonian_onebody`) to compute the matrix elements of the ECP operators, which then become part of the effective core potential in the Kohn-Sham or Hartree-Fock equations.
*   **Numerical Integration Grid (`m_dft_grid`)**: The parameters `nradial_ecp` and `nangular_ecp` set here define the quality of the grid used for integrating the non-local ECP projectors if they are not of GTH type (which have their own integration schemes).
```
