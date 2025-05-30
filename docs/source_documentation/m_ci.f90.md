# `m_ci.f90`

## Overview

The `m_ci` Fortran module in MOLGW implements Configuration Interaction (CI) calculations. CI is a post-Hartree-Fock method used to determine the electronic structure of atoms and molecules by constructing the N-electron wavefunction as a linear combination of Slater determinants (configurations). This module focuses on "full CI" within a defined active space of molecular orbitals, suitable for systems with a small number of active electrons. It provides functionalities for generating CI configurations, building and diagonalizing the CI Hamiltonian matrix, and using the resulting CI wavefunctions and energies to compute properties such as self-energies via Green's function methods. Occupation number vectors (bit-encoded integers) are used to efficiently represent and manipulate Slater determinants.

## Key Components

*   **`configurations` (Type)**:
    *   A derived data type that stores essential information about a set of CI configurations for a given electronic state.
    *   Members include: `nelec` (total electrons), `nelec_valence` (active electrons), `nelec_valence_up`/`down`, `nconf` (number of configurations), `sz` (target S<sub>z</sub> value), and `keyud(:, :)` (allocatable array of 2x`nconf` integers of kind `key_int`, where each pair `keyud(1,iconf)` and `keyud(2,iconf)` stores the bit-encoded occupation numbers for spin-up and spin-down orbitals for configuration `iconf`).

*   **`prepare_ci(basis, nstate_in, nfrozen_in, c_matrix)`**:
    *   Initializes common data for CI calculations, such as the number of active orbitals (`nstate_ci - nfrozen_ci`) and the number of frozen core orbitals (`nfrozen_ci`).
    *   Computes and stores the one-electron Hamiltonian (`h_1body`) in the molecular orbital (MO) basis provided by `c_matrix`.

*   **`setup_configurations_ci(nelec, spinstate, ci_type_in, conf)`**:
    *   Generates the list of Slater determinants (configurations) to be included in the CI expansion.
    *   `nelec`, `spinstate`: Target number of electrons and S<sub>z</sub>.
    *   `ci_type_in`: String specifying the CI level (e.g., 'ALL' for FCI within active space, 'CISD', 'CISDTQ').
    *   Populates the `conf%keyud` array with the bit-encoded occupation number vectors for all selected configurations.
    *   Can read reference configurations from an external file (`manual_ci_ref_xx`).

*   **`hamiltonian_ci(keyudi, keyudj)` (Function)**:
    *   Calculates the Hamiltonian matrix element H<sub>ij</sub> = <D<sub>i</sub>|H|D<sub>j</sub>> between two Slater determinants D<sub>i</sub> and D<sub>j</sub>, represented by their occupation number vectors `keyudi` and `keyudj`.
    *   Uses Slater-Condon rules to evaluate matrix elements involving one- and two-electron integrals (`h_1body` and `eri_eigen`).
    *   Handles cases where determinants are identical (diagonal elements), differ by one spin-orbital, or differ by two spin-orbitals. Matrix elements are zero if they differ by more than two spin-orbitals.

*   **`build_ci_hamiltonian(conf, desc_hci, h_ci)`**:
    *   Constructs the full CI Hamiltonian matrix. It iterates over pairs of configurations from `conf%keyud` and calls `hamiltonian_ci` to compute each element.
    *   Supports distributed parallel construction using ScaLAPACK descriptors (`desc_hci`).

*   **`build_ci_hamiltonian_sparse(conf, desc, h)`**:
    *   Constructs the CI Hamiltonian in a compressed sparse row (CSR) or similar sparse format (`type(sparse_matrix)`). This is useful for very large CI expansions where storing the full matrix is infeasible.

*   **`full_ci_nelectrons(save_coefficients, nelectron, spinstate, nuc_nuc)`**:
    *   The main driver subroutine for performing a CI calculation.
    *   `save_coefficients`: Integer flag (0 for N-electron, -1 for N+1, +1 for N-1) to determine which global CI object (`conf_0`, `conf_m`, `conf_p`) to populate.
    *   Sets up configurations, builds the CI matrix (full or sparse), and then diagonalizes it using either full diagonalization (LAPACK/ScaLAPACK) or the Davidson iterative method (`diagonalize_davidson_ci`).
    *   Stores the resulting eigenvalues (energies) and eigenvectors (CI coefficients).
    *   Handles reading/writing of eigenvectors for restarts.

*   **`full_ci_nelectrons_selfenergy(energy_gks)`**:
    *   Calculates the self-energy using CI results for N, N-1, and N+1 electron systems (stored in `conf_0`, `conf_p`, `conf_m`).
    *   Computes Lehmann amplitudes (overlaps like <N|a<sup>+</sup>|N-1>) and energies for poles of the Green's function.
    *   Outputs the exact Green's function and self-energy on a frequency grid.

*   **Occupation Number Vector Utilities**:
    *   `gamma_sign_keyud`: Calculates the sign factor from fermionic anti-commutation rules.
    *   `get_keyud`: Converts lists of occupied spin-orbitals into a bit-encoded key.
    *   `increment_sporb`: Increments an occupation list to the next valid configuration.
    *   `get_spinz_from_keyud`: Calculates S<sub>z</sub> from a key.
    *   `get_spins_from_keyud`, `get_states_from_keyud`: Decode a key back into lists of occupied spin-orbitals.
    *   `keyud_diff_order`, `key_diff_order`, `keysud_diff_order`, `keys_diff_order`: Functions to determine the number of differing spin-orbitals between configurations.

*   **`diagonalize_davidson_ci(...)`**:
    *   Implements the Davidson iterative algorithm for finding a few lowest eigenvalues and eigenvectors of the (potentially sparse) CI matrix.

## Important Variables/Constants

*   **`key_int` (Integer, Parameter)**: Kind parameter for integers used as occupation number vectors (bit strings). Its size (e.g., 8 for 64-bit integers) limits the number of active orbitals that can be represented.
*   **`nfrozen_ci` (Integer)**: Number of core orbitals excluded from the CI active space.
*   **`nstate_ci` (Integer)**: Total number of molecular orbitals considered (frozen + active).
*   **`h_1body(:, :)` (Real(dp), Allocatable)**: Stores the one-electron Hamiltonian (kinetic + nuclear attraction + ECPs) in the MO basis.
*   **`eri_eigen` (Function from `m_eri_ao_mo`)**: Used to retrieve two-electron repulsion integrals <pq|rs> in the MO basis.
*   **`conf_0`, `conf_p`, `conf_m` (Type(configurations))**: Global objects storing CI configuration data for N, N-1 (hole/cationic), and N+1 (electron/anionic) systems, respectively.
*   **`energy_0`, `energy_p`, `energy_m` (Real(dp), Allocatable)**: Arrays storing CI eigenvalues for N, N-1, N+1 systems.
*   **`eigvec_0`, `eigvec_p`, `eigvec_m` (Real(dp), Allocatable)**: Arrays storing CI eigenvectors for N, N-1, N+1 systems.

## Usage Examples

The `m_ci` module is typically invoked when a CI calculation is requested via the `postscf` input variable in MOLGW.
```fortran
! Conceptual flow within MOLGW:

! 1. Perform SCF (HF or DFT) to get MOs (c_matrix) and orbital energies.
! CALL run_scf(...)

! 2. Prepare for CI calculations
CALL prepare_ci(orbital_basis, total_mos, num_frozen_core_mos, mo_coefficients)

! 3. Perform CI for the N-electron ground state (and possibly excited states)
CALL full_ci_nelectrons(0, num_electrons_neutral, spin_multiplicity_neutral, nuclear_repulsion_energy)

! 4. If self-energy/Green's function is needed:
!    Perform CI for (N-1) electron states
CALL full_ci_nelectrons(1, num_electrons_neutral - 1, spin_multiplicity_cation, nuclear_repulsion_energy)
!    Perform CI for (N+1) electron states
CALL full_ci_nelectrons(-1, num_electrons_neutral + 1, spin_multiplicity_anion, nuclear_repulsion_energy)

! 5. Calculate self-energy
CALL full_ci_nelectrons_selfenergy(kohn_sham_energies)

! 6. Clean up
CALL destroy_ci()
```

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_mpi`, `m_warning`, `m_linear_algebra`, `m_memory`, `m_timing`, `m_scalapack`, `m_inputparam`.
*   **Basis Set & Integrals**: `m_basis_set` (for basis info if needed, though mostly MO-based), `m_eri_ao_mo` (crucial for providing MO two-electron integrals `eri_eigen`).
*   **Hamiltonian**: `m_hamiltonian_onebody` (used in `prepare_ci` to set up the one-electron part of the Hamiltonian in the AO basis before transformation to MO).
*   **Self-Energy**: `m_selfenergy_tools` (defines `selfenergy_grid` type used in `full_ci_nelectrons_selfenergy`).
*   **Fortran Intrinsics**: Uses bit manipulation intrinsics extensively (e.g., `POPCNT`, `IEOR`, `BTEST`, `SHIFTL`, `TRAILZ`, `LEADZ`).
*   **LAPACK/ScaLAPACK**: For diagonalization of the CI matrix.
*   The module forms a high-level computational capability. Its accuracy is determined by the size of the active space and the CI expansion level.
```
