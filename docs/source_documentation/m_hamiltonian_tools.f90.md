# `m_hamiltonian_tools.f90`

## Overview

The `m_hamiltonian_tools` Fortran module in MOLGW provides a suite of utility subroutines and functions for manipulating, analyzing, and transforming matrices and orbitals commonly encountered in electronic structure calculations. This includes routines for constructing density matrices, transforming matrices between Atomic Orbital (AO) and Molecular Orbital (MO) bases, setting orbital occupations (e.g., via Fermi-Dirac statistics for finite temperature calculations), performing level shifting, calculating S<sup>2</sup> expectation values, and diagonalizing Hamiltonians (including generalized eigenvalue problems).

## Key Components

*   **`get_number_occupied_states(occupation)` (Function)**:
    *   Determines the highest index of an orbital considered occupied based on the `occupation` number array.

*   **`setup_density_matrix` (Interface for `setup_density_matrix_real` and `setup_density_matrix_cmplx`)**:
    *   **`setup_density_matrix_real(c_matrix, occupation, p_matrix)`**: Constructs the real one-particle density matrix (1-RDM) P<sub>&mu;&nu;</sub> = &sum;<sub>i</sub> C<sub>&mu;i</sub> n<sub>i</sub> C<sub>&nu;i</sub><sup>*</sup> in the AO basis.
    *   **`setup_density_matrix_cmplx(c_matrix_cmplx, occupation, p_matrix_cmplx)`**: Constructs the complex 1-RDM.

*   **Density Matrix Transformation Routines**:
    *   **`setup_density_matrix_MO_cmplx(c_matrix, s_matrix, p_matrix_cmplx, p_matrix_MO_cmplx)`**: Transforms a complex AO density matrix P<sup>AO</sup> to the MO basis: P<sup>MO</sup> = C<sup>&dagger;</sup> S P<sup>AO</sup> S C, where S is the AO overlap matrix.
    *   **`setup_density_matrix_MO_real(c_matrix, s_matrix, p_matrix_real, p_matrix_MO_real)`**: Same for a real AO density matrix.

*   **`setup_energy_density_matrix(c_matrix, occupation, energy, q_matrix)`**:
    *   Constructs the energy-weighted density matrix Q<sub>&mu;&nu;</sub> = &sum;<sub>i</sub> C<sub>&mu;i</sub> n<sub>i</sub> &epsilon;<sub>i</sub> C<sub>&nu;i</sub><sup>*</sup> in the AO basis.

*   **`test_density_matrix(p_matrix, s_matrix)`**:
    *   A debugging routine to check the idempotency condition PSP=P (or PSP=PS for fractional occupations, though the comment says it's for integer occupations).

*   **`set_occupation(temperature, electrons_in, magnetization, energy, occupation)`**:
    *   Sets the molecular orbital `occupation` numbers.
    *   At zero `temperature`, it applies Aufbau filling based on `electrons_in` and `magnetization`.
    *   At finite `temperature`, it uses Fermi-Dirac statistics, iteratively finding the chemical potential `mu` to match `electrons_in`.
    *   Can also read occupations directly from a file named `manual_occupations`.

*   **Output/Dump Routines**:
    *   `dump_out_occupation(...)`: Prints formatted occupation numbers.
    *   `dump_out_energy(...)`: Prints formatted orbital energies along with occupations.
    *   `dump_out_energy_yaml(...)`: Prints energies in YAML format.
    *   `output_homolumo(...)`: Calculates and prints HOMO, LUMO energies and the HOMO-LUMO gap.

*   **AO/MO Transformation Utilities**:
    *   `matrix_ao_to_mo_diag(c_matrix, matrix_in, diag_out)`: Transforms an AO-basis matrix (`matrix_in`) to the MO basis and returns only its diagonal elements: diag_out<sub>i</sub> = (C<sup>&dagger;</sup> H<sup>AO</sup> C)<sub>ii</sub>.
    *   `matrix_ao_to_mo(c_matrix, matrix_in, matrix_out)`: Performs a full transformation of an AO-basis matrix to the MO basis: H<sup>MO</sup> = C<sup>&dagger;</sup> H<sup>AO</sup> C.
    *   `matrix_mo_to_ao(c_matrix, matrix_in, matrix_out)`: Performs a full transformation of an MO-basis matrix to the AO basis: H<sup>AO</sup> = C H<sup>MO</sup> C<sup>&dagger;</sup> (note: this assumes C is unitary or S-orthogonal, otherwise S<sup>-1</sup> factors are needed for a proper back-transformation of operators).

*   **`evaluate_s2_operator(occupation, c_matrix, s_matrix)`**:
    *   Calculates the expectation value of the total spin squared operator, <S<sup>2</sup>>, for unrestricted Hartree-Fock or Kohn-Sham wavefunctions.

*   **Level Shifting Routines**:
    *   `level_shifting_up(...)`: Applies a positive energy shift to virtual orbitals to aid SCF convergence. Modifies the Hamiltonian matrix.
    *   `level_shifting_down(...)`: Applies a negative energy shift to virtual orbitals (or positive to occupied, depending on implementation details not fully clear but typically virtuals are shifted down if occupieds are fixed). Modifies orbital energies and the Hamiltonian.

*   **Basis Orthonormalization/Transformation**:
    *   `setup_x_matrix(TOL_OVERLAP, s_matrix, nstate, x_matrix)`: Constructs the transformation matrix X = U s<sup>-1/2</sup>, where S = U s U<sup>&dagger;</sup> is the eigendecomposition of the AO overlap matrix `s_matrix`. X is used to transform to an orthonormal AO basis. Orbitals corresponding to eigenvalues of S smaller than `TOL_OVERLAP` are removed to handle linear dependencies. `nstate` becomes the number of retained (orthonormal) basis functions.
    *   `setup_sqrt_overlap(s_matrix, s_matrix_sqrt)`: Calculates S<sup>1/2</sup>.
    *   `setup_sqrt_density_matrix(p_matrix, p_matrix_sqrt, p_matrix_occ)`: Calculates P<sup>1/2</sup> by diagonalizing P, taking the square root of its eigenvalues (occupations), and transforming back.
    *   `get_c_matrix_from_p_matrix(p_matrix, c_matrix, occupation)`: Diagonalizes a given density matrix `p_matrix` to obtain its eigenvectors (natural orbitals `c_matrix`) and eigenvalues (natural orbital `occupation` numbers).

*   **`diagonalize_hamiltonian_scalapack(hamiltonian, x_matrix, energy, c_matrix)`**:
    *   Solves the generalized eigenvalue problem H C' = S C' E by transforming it to a standard eigenvalue problem H' C'' = C'' E, where H' = X<sup>&dagger;</sup> H X and C = X C''. Here, `x_matrix` is typically S<sup>-1/2</sup> or the transformation to an orthonormal basis.
    *   Uses ScaLAPACK for parallel diagonalization if available.

## Important Variables/Constants

*   **`completely_empty` (Real(dp), from `m_definitions`)**: A small threshold (e.g., 1.0e-5) used to determine if an orbital is effectively unoccupied.
*   **`TOL_OVERLAP` (Real(dp), Intent(in))**: Threshold for eigenvalues of the overlap matrix when constructing the `x_matrix` (S<sup>-1/2</sup> transformation), used to discard functions causing linear dependencies.
*   **`scf_diago_flavor` (from `m_inputparam`)**: A string that can influence the choice of diagonalization algorithm (e.g., 'STD' for standard, 'DC' for divide-and-conquer) in LAPACK/ScaLAPACK calls.

## Usage Examples

This module is a utility library used extensively throughout MOLGW.

Setting occupations at finite temperature:
```fortran
USE m_hamiltonian_tools
USE m_definitions
IMPLICIT NONE
REAL(DP) :: orbital_energies(100, 1), occupations(100, 1)
REAL(DP) :: temperature_K, temperature_au, total_electrons, current_magnetization
! ... initialize orbital_energies, temperature_K, total_electrons, current_magnetization ...
temperature_au = temperature_K / Ha_K
CALL set_occupation(temperature_au, total_electrons, current_magnetization, orbital_energies, occupations)
```

Transforming a Fock matrix from AO to MO basis:
```fortran
USE m_hamiltonian_tools
USE m_definitions
IMPLICIT NONE
REAL(DP) :: F_ao(100,100,1), C_mo(100,50,1), F_mo(50,50,1)
! ... F_ao and C_mo are populated ...
CALL matrix_ao_to_mo(C_mo, F_ao, F_mo)
! F_mo now contains the Fock matrix in the MO basis
```

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_timing`, `m_mpi`, `m_scalapack`, `m_warning`, `m_memory`, `m_inputparam`.
*   **Basis Set**: `m_basis_set` (for basis set information, though many routines here operate on matrices assuming dimensions are already consistent).
*   **Linear Algebra**: `m_linear_algebra` (provides `diagonalize`, `matrix_lower_to_full` for serial operations).
*   **BLAS/LAPACK/ScaLAPACK**: This module relies heavily on these libraries for efficient matrix operations such as matrix multiplication (`DSYRK`, `ZHERK`, `DSYMM`, `ZHEMM`, `DGEMM`, `PDGEMM`) and diagonalization (`DSYEV`, ScaLAPACK equivalents). MKL-specific `DGEMMT` is also used.
*   **Other MOLGW Modules**:
    *   `m_scf`: Likely calls routines from this module for density matrix construction, occupation setting, and diagonalization.
    *   Property modules might use AO/MO transformation routines.
    *   Any module dealing with Hamiltonians or density matrices in different representations will likely use these tools.
```
