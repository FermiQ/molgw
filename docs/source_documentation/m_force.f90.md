# `m_force.f90`

## Overview

The `m_force` Fortran module in MOLGW is dedicated to calculating the forces acting on the atomic nuclei. These forces are essential for performing geometry optimizations (finding equilibrium structures) and ab initio molecular dynamics (AIMD) simulations. The calculation relies on the availability of analytical gradients of the one- and two-electron integrals with respect to nuclear coordinates, which are typically provided by an underlying integral library like Libint.

## Key Components

*   **`calculate_force(basis, occupation, energy, c_matrix)`**:
    *   **Purpose**: This is the main subroutine that orchestrates the computation of the total atomic forces. It assembles forces from various contributions according to the Hellmann-Feynman theorem and Pulay corrections (terms arising from the derivative of the basis set with respect to nuclear coordinates).
    *   **Workflow**:
        1.  **Prerequisite Check**: Verifies if the linked integral library supports gradient calculations (`MOLGW_has_gradient`). If not, it issues a warning and returns.
        2.  **Density Matrix Setup**: Constructs the one-particle density matrix (`p_matrix`) and the energy-weighted density matrix (`r_matrix`) from the molecular orbital coefficients (`c_matrix`) and occupations/energies.
        3.  **Screening**: Implements a screening mechanism (`skip_shellpair`) to identify pairs of atomic orbital shells whose contributions to the two-electron part of the force can be neglected if the corresponding density matrix elements are below `TOL_DENSITY_MATRIX`.
        4.  **Force Component Calculation**:
            *   **Nucleus-Nucleus Repulsion (`force_nuc_nuc`)**: Calls `nucleus_nucleus_force` (from `m_atoms`).
            *   **Overlap Pulay Forces (`force_ovp`)**: Calculates terms involving derivatives of overlap integrals, &sum;<sub>&mu;&nu;</sub> R<sub>&mu;&nu;</sub> &nabla;S<sub>&mu;&nu;</sub>. Uses `setup_overlap_grad`.
            *   **Kinetic Energy Forces (`force_kin`)**: Calculates terms involving derivatives of kinetic energy integrals, &sum;<sub>&mu;&nu;</sub> P<sub>&mu;&nu;</sub> &nabla;T<sub>&mu;&nu;</sub>. Uses `setup_kinetic_grad`.
            *   **Electron-Nucleus Forces (`force_nuc`)**: Calculates terms involving derivatives of electron-nucleus attraction integrals. This includes both Hellmann-Feynman contributions (where the derivative acts on the 1/r operator) and Pulay terms (derivative acts on basis functions). Uses `setup_nucleus_grad`. The Hellmann-Feynman part is stored in `force_hellfeyn`.
            *   **Two-Electron (Hartree & Exchange) Forces (`force_har`, `force_exx`)**: Calculates terms involving derivatives of 4-center two-electron repulsion integrals (ERIs). This is the most computationally intensive part. It iterates over shell quartets, calls `calculate_eri_4center_shell_grad` (from `m_eri_calculate`) to get ERI derivatives, and contracts them with appropriate products of density matrix elements. The exchange contribution is scaled by `alpha_hybrid`.
        5.  **Total Force Assembly**: Sums all the above contributions to get the total force on each atom, stored in `force(:, :)` (from `m_atoms`).
        6.  **Output**: Prints a breakdown of Hellmann-Feynman forces, Pulay forces, and total forces for each atom.

## Important Variables/Constants

*   **Input Parameters to `calculate_force`**:
    *   `basis` (Type `basis_set`): Basis set information.
    *   `occupation(:, :)` (Real(dp)): Molecular orbital occupation numbers.
    *   `energy(:, :)` (Real(dp)): Molecular orbital energies.
    *   `c_matrix(:, :, :)` (Real(dp)): Molecular orbital coefficients.
*   **Force Arrays (mostly from `m_atoms`)**:
    *   `force(:, :)`: Stores the final total Cartesian forces on each atom.
    *   `force_nuc_nuc(:, :)`: Nucleus-nucleus repulsion component.
    *   `force_ovp(:, :)`: Overlap derivative (Pulay) component.
    *   `force_kin(:, :)`: Kinetic energy derivative component.
    *   `force_nuc(:, :)`: Electron-nucleus attraction derivative component.
    *   `force_har(:, :)`: Hartree (two-electron Coulomb) derivative component.
    *   `force_exx(:, :)`: Exchange (two-electron) derivative component.
    *   `force_hellfeyn(:, :)`: Pure Hellmann-Feynman component of electron-nucleus forces.
*   **`MOLGW_has_gradient` (Logical, from `m_definitions`)**: Global flag indicating if the linked integral engine provides gradients. Essential for this module's operation.
*   **`TOL_DENSITY_MATRIX` (Real(dp), Parameter)**: Threshold (1.0e-2) used to screen out contributions from shell pairs with small density matrix elements, aiming to reduce the cost of ERI derivative calculations.
*   **`alpha_hybrid` (Real(dp), from `m_inputparam`)**: Fraction of exact exchange included, used to scale the exchange force contribution.

## Usage Examples

The `calculate_force` subroutine is called internally by MOLGW when a task requiring atomic forces is performed, such as:
*   Geometry optimization (`task = OPTIMIZE`)
*   Ab initio molecular dynamics
*   Calculation of vibrational frequencies (via finite differences of forces or analytical second derivatives if available elsewhere).

A typical call sequence would be after an SCF (or other electronic structure) calculation has converged:
```fortran
! Conceptual: Within a geometry optimization step in MOLGW
USE m_force
USE m_atoms
USE m_basis_set
! ... (SCF converged, density matrix P and energy-weighted DM R are available) ...
! ... (basis, occupation_scf, energy_scf, c_matrix_scf are populated) ...

CALL calculate_force(basis, occupation_scf, energy_scf, c_matrix_scf)

! The 'force' array in m_atoms module now contains the updated forces.
! This can be used by an optimizer (e.g., LBFGS in m_lbfgs via m_atoms%relax_atoms)
! to update atomic positions.
```

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_warning`, `m_timing`, `m_inputparam`.
*   **Atomic Structure**: `m_atoms` (provides storage for force components and atomic coordinates).
*   **Basis Set & ERI**: `m_basis_set` (basis information), `m_eri` (shell pair indexing and screening data from `negligible_shellpair`, `index_shellpair`), `m_eri_calculate` (provides `calculate_eri_4center_shell_grad` for ERI derivatives).
*   **Hamiltonian Components**: `m_hamiltonian_tools` (for `setup_density_matrix`, `setup_energy_density_matrix`), `m_hamiltonian_onebody` (provides routines like `setup_overlap_grad`, `setup_kinetic_grad`, `setup_nucleus_grad` which compute one-electron integral derivatives).
*   **Integral Engine (Libint/Libcint)**: Critically depends on the underlying integral library's ability to provide analytical gradients of one- and two-electron integrals. This is interfaced via `m_libint_tools` or `m_libcint_tools` which are called by `m_hamiltonian_onebody` and `m_eri_calculate`.
*   **Geometry Optimizer**: The computed forces are used by geometry optimization algorithms (e.g., in `m_lbfgs` as called by `m_atoms%relax_atoms`) to predict new atomic positions.
```
