# `m_dm_mbpt.f90`

## Overview

The `m_dm_mbpt` Fortran module in MOLGW is responsible for obtaining and utilizing one-body reduced density matrices (1-RDMs) derived from Many-Body Perturbation Theory (MBPT). It can calculate the 1-RDM based on various MBPT methods (e.g., Møller-Plesset perturbation theory like MP2, or GW approximation) or read a pre-existing density matrix from external files (like Gaussian fchk or a MOLGW-specific format). Once the correlated density matrix is obtained, this module can facilitate the calculation of various electronic properties, generate outputs for visualization (like cube files), and optionally update the Fock matrix for use in self-consistent procedures that incorporate correlation effects via the density matrix.

## Key Components

*   **`get_dm_mbpt(basis, occupation, energy, c_matrix, s_matrix, hamiltonian_kinetic, hamiltonian_nucleus, hamiltonian_fock)`**:
    *   **Purpose**: This is the main subroutine that orchestrates the acquisition and subsequent processing of the correlated density matrix.
    *   **Workflow**:
        1.  **Obtain Correlated Density Matrix (`p_matrix_corr`)**:
            *   If `read_fchk` (input parameter) is specified, it attempts to read `p_matrix_corr` from a Gaussian fchk file.
            *   If `pt_density_matrix` (input parameter) specifies a method (e.g., 'ONE-RING', 'PT2', 'GW', 'EVGW', 'GW_IMAGINARY', 'GW_DYSON', 'HF', 'HF_SECOND_ORDER'), it calls the appropriate subroutine to calculate `p_matrix_corr`. This may involve:
                *   `fock_density_matrix`: First-order perturbative correction using the Fock matrix.
                *   `onering_density_matrix`, `pt2_density_matrix`: Density matrices from specific second-order theories.
                *   `gw_density_matrix`, `gw_density_matrix_imag`, `gw_density_matrix_dyson_imag`: Density matrices derived from GW self-energy calculations, utilizing polarizabilities from `m_linear_response`.
                *   `fock_density_matrix_second_order`: Second-order perturbative correction.
            *   If `p_matrix_corr` has not been obtained by the above methods, it attempts to read it from a MOLGW binary file named `DENSITY_MATRIX`.
        2.  **Natural Orbital Analysis**: Diagonalizes the obtained `p_matrix_corr` (in the MO basis, after transformation `S*P*S`) to get natural orbital occupation numbers (`natural_occupation`) and the corresponding natural orbital coefficients (`c_matrix_tmp` in AO basis).
        3.  **Optional Outputs**:
            *   If `print_cube_`, `print_wfn_`, `print_wfn_files_` flags are set, it generates cube files, density plots, or wavefunction files using the natural orbitals and occupations.
            *   If `print_hartree_` or `use_correlated_density_matrix_` flags are set, it calculates various energy components (kinetic, nuclear attraction, Hartree, exchange) using `p_matrix_corr`.
            *   If `print_multipole_` is set, it computes static electric dipole and quadrupole moments.
        4.  **Update Fock Matrix**: If `use_correlated_density_matrix_` is true, the input `hamiltonian_fock` is updated with a new Fock matrix constructed from `p_matrix_corr` (including its Hartree and exchange contributions).

*   **`fock_density_matrix(basis, occupation, energy, c_matrix, hfock, p_matrix)`**:
    *   **Purpose**: Calculates the first-order perturbative correction to the density matrix. The off-diagonal elements (occupied-virtual block) are computed as P<sup>(1)</sup><sub>ia</sub> = F<sub>ia</sub> / (&epsilon;<sub>i</sub> - &epsilon;<sub>a</sub>), where F is the Fock matrix and &epsilon; are orbital energies. The diagonal elements are the input `occupation` numbers.
    *   Stores the result in `p_matrix` (AO basis) after transformation from MO basis.

*   **`fock_density_matrix_second_order(basis, occupation, energy, c_matrix, hfock, p_matrix)`**:
    *   **Purpose**: Calculates the second-order perturbative correction to the density matrix. This involves terms derived from the difference between the Fock matrix (`hfock`) and the Kohn-Sham potential (implicit in `energy`), effectively using (Sigma - V<sub>xc/HF</sub>) as the perturbation.
    *   Adds the computed correction to the input `p_matrix`.

## Important Variables/Constants

*   **`p_matrix_corr(:, :, :)` (Real(dp), Allocatable)**: The primary working array holding the correlated one-body reduced density matrix in the atomic orbital (AO) basis.
*   **`pt_density_matrix` (Character string from `m_inputparam`)**: Input parameter controlling which MBPT method is used to calculate `p_matrix_corr`. Examples: 'PT2', 'GW', 'HF_SECOND_ORDER'.
*   **`read_fchk` (Character string from `m_inputparam`)**: Input parameter specifying the path to a Gaussian formatted checkpoint file from which to read the density matrix.
*   **`use_correlated_density_matrix_` (Logical from `m_inputparam`)**: If true, the main Fock matrix of the SCF cycle (`hamiltonian_fock`) is updated using the correlated density matrix computed in this module. This allows for a degree of self-consistency with the correlated method.
*   **Input Arrays**: `occupation`, `energy`, `c_matrix` (MO coefficients from a preceding SCF calculation), `s_matrix` (AO overlap matrix), `hamiltonian_kinetic`, `hamiltonian_nucleus`, `hamiltonian_fock` (Fock matrix from SCF).
*   **Natural Orbitals/Occupations**: `natural_occupation`, `c_matrix_tmp` are used to store and process the natural orbitals derived from `p_matrix_corr`.

## Usage Examples

This module is typically invoked internally within MOLGW's workflow, usually after an initial SCF calculation, if specific options are set in the input file.

Example scenario in MOLGW input:
```
task = ENERGY
...
scf = HF
...
pt_density_matrix = GW  ! Request GW density matrix calculation
use_correlated_density_matrix = YES ! Update Fock matrix with this DM
print_multipole = YES
```
In this case, after the HF SCF, `get_dm_mbpt` would be called. It would trigger a GW calculation to get the GW self-energy, from which the GW density matrix is computed. This density matrix would then be used to update the Fock matrix (potentially for a subsequent "GW0-like" cycle if MOLGW supports it) and to calculate multipole moments.

## Dependencies and Interactions

*   **Core/Utility Modules**: `m_definitions`, `m_timing`, `m_warning`, `m_memory`, `m_atoms`, `m_basis_set`, `m_inputparam`.
*   **SCF & Hamiltonian**: `m_scf` (for initial electronic structure data and DM diagonalization), `m_hamiltonian_tools`, `m_hamiltonian_wrapper`.
*   **MBPT/GW Components**:
    *   `m_spectral_function` (`wpol` type for polarizabilities).
    *   `m_selfenergy_tools`.
    *   `m_gw_selfenergy_grid`.
    *   `m_linear_response` (for calculating polarizability needed for GW DM).
    *   `m_pt_density_matrix` (provides specific subroutines like `onering_density_matrix`, `gw_density_matrix`).
*   **Integral Calculations**: `m_eri_ao_mo` (for transforming integrals to MO basis if needed by specific DM calculation routines).
*   **Property Calculations**: `m_multipole` (for `static_dipole`, `static_quadrupole`).
*   **I/O**: `m_io` (implicitly for `read_gaussian_fchk`, `plot_cube_wfn`, `print_wfn_file`).
*   **File System**: May read from `'gaussian.fchk'` or `'DENSITY_MATRIX'`. May write cube files or other analysis outputs.
```
