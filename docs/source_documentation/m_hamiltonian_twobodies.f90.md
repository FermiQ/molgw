# `m_hamiltonian_twobodies.f90`

## Overview

The `m_hamiltonian_twobodies` Fortran module in MOLGW is responsible for constructing the two-body components of the Fock or Kohn-Sham matrix. This primarily involves calculating the Hartree potential (Coulomb term) and the exchange potential (for Hartree-Fock or hybrid DFT calculations) from the two-electron repulsion integrals (ERIs) and the density matrix. It also includes routines for computing the exchange-correlation (XC) potential and energy by numerical integration on a DFT grid, interfacing with the Libxc library for the XC functional definitions. Additionally, this module provides methods to initialize the molecular orbital (MO) coefficient matrix using various guess Hamiltonians or by reading from external sources like Gaussian fchk files.

## Key Components

*   **Hartree Potential Calculation**:
    *   **`setup_hartree(p_matrix, hartree_ao, ehartree)`**: Computes the Hartree potential matrix V<sup>H</sup><sub>&mu;&nu;</sub> = &sum;<sub>&lambda;&sigma;</sub> P<sub>&lambda;&sigma;</sub> (&mu;&nu;|&lambda;&sigma;) and the Hartree energy E<sub>H</sub> using pre-calculated 4-center ERIs. Employs 8-fold permutation symmetry for efficiency.
    *   **`setup_hartree_oneshell(basis, p_matrix, hartree_ao, ehartree)`**: An "out-of-core"-like approach for the Hartree potential, where 4-center ERIs are computed for one shell quartet at a time, contracted with the density matrix, and then discarded. This is more memory-efficient for very large systems if ERIs are not stored in memory.
    *   **`setup_hartree_ri(p_matrix, hartree_ao, ehartree)`**: Computes the Hartree potential and energy using the Resolution of Identity (RI) approximation: V<sup>H</sup><sub>&mu;&nu;</sub> = &sum;<sub>Q,P</sub> (&mu;&nu;|Q) (Q|P)<sup>-1</sup> (P|&lambda;&sigma;) P<sub>&lambda;&sigma;</sub>. Requires 3-center and 2-center ERIs involving an auxiliary basis.
    *   **`calculate_density_auxilbasis(p_matrix, rho_coeff)`**: Computes the coefficients R<sub>I</sub> of the electron density when expanded in the auxiliary basis set, used in conjunction with "genuine" 3-center RI.
    *   **`setup_hartree_genuine_ri(p_matrix, rho_coeff, hartree_ao, ehartree)`**: Computes the Hartree potential using R<sub>I</sub> and genuine 3-center integrals (Q|&mu;&nu;): V<sup>H</sup><sub>&mu;&nu;</sub> = &sum;<sub>Q</sub> (&mu;&nu;|Q) R<sub>Q</sub>.

*   **Exchange Potential Calculation**:
    *   **`setup_exchange(p_matrix, exchange_ao, eexchange)`**: Computes the exact exchange potential matrix K<sub>&mu;&nu;</sub> = -&sum;<sub>&lambda;&sigma;</sub> P<sub>&lambda;&sigma;</sub> (&mu;&lambda;|&nu;&sigma;) and the exchange energy using 4-center ERIs.
    *   **`setup_exchange_longrange(p_matrix, exchange_ao, eexchange)`**: Similar to `setup_exchange`, but uses the long-range part of 4-center ERIs (for range-separated hybrids).
    *   **`setup_exchange_ri(...)`**: Computes the exchange potential and energy using the RI approximation. Handles real and complex (`_cmplx`) MO coefficients, as well as specialized versions for x2c calculations (`_x2c_1`, `_x2c_2`).
    *   **`setup_exchange_longrange_ri(...)`**: Similar for long-range RI exchange.

*   **DFT Exchange-Correlation (XC) Calculation**:
    *   **`dft_exc_vxc_batch(batch_size, basis, occupation, c_matrix, vxc_ao, exc_xc, ...)`**: Calculates the XC potential matrix V<sup>XC</sup><sub>&mu;&nu;</sub> and the XC energy E<sub>XC</sub> by numerical integration on a DFT grid.
        *   It evaluates the electron density (&rho;) and its gradient (&nabla;&rho;, for GGAs) on grid batches.
        *   Calls routines from `m_libxc_tools` (which interfaces Libxc) to get XC energy densities and potential derivatives (v<sub>&rho;</sub>, v<sub>&sigma;</sub>) on the grid.
        *   Constructs V<sup>XC</sup><sub>&mu;&nu;</sub> by integrating basis functions with these XC potential derivatives.
        *   Supports both real and complex `c_matrix` (via `CLASS(*)`).

*   **Initial Guess for Molecular Orbitals**:
    *   **`dft_approximate_vhxc(basis, vhxc_ao)`**: Computes an approximate sum of Hartree and XC potentials. The density is approximated as a sum of atomic densities, and a simple LDA functional is used on a coarse grid. Used for initial guesses.
    *   **`init_c_matrix(basis, occupation, x_matrix, hkin, hnuc, c_matrix)`**: Initializes the MO coefficient matrix `c_matrix`.
        *   `init_hamiltonian` (input parameter) controls the method:
            *   'GUESS' or 'MIX': Diagonalizes H<sub>core</sub> + V<sub>approx_HXC</sub>. 'MIX' additionally mixes HOMO/LUMO for spin-compensated systems.
            *   'CORE' or 'NOFT': Diagonalizes H<sub>core</sub> = T + V<sub>nuc</sub>.
            *   'GAUSSIAN': Attempts to read MO coefficients from a Gaussian fchk file. Falls back to 'CORE' if issues arise.
    *   **`init_c_matrix_cmplx(...)`**: Initializes complex MOs by applying random phases to real MOs.
    *   **`init_c_matrix_x2c(...)`**: Initializes MOs for x2c calculations, potentially using an approximate Hamiltonian.

## Important Variables/Constants

*   **`p_matrix` (Real(dp) or Complex(dp))**: The one-body reduced density matrix in the AO basis.
*   **`hartree_ao(:, :)` (Real(dp))**: AO-basis Hartree potential matrix.
*   **`exchange_ao(:, :, :)` (Real(dp) or Complex(dp))**: AO-basis exchange potential matrix (per spin).
*   **`vxc_ao(:, :, :)` (Real(dp))**: AO-basis XC potential matrix (per spin).
*   **`ehartree`, `eexchange`, `exc_xc` (Real(dp))**: Hartree, exchange, and XC energy contributions.
*   **`eri_4center(:)` (from `m_eri`)**: Array of 4-center AO ERIs.
*   **`eri_3center(:, :)` (from `m_eri`)**: Array of 3-center AO ERIs (Aux|AO AO) or RI-transformed integrals.
*   **`eri_2center_inv(:, :)` (from `m_eri_calculate`)**: Inverse of the 2-center auxiliary basis ERI matrix, (Q|P)<sup>-1</sup>, used in RI.
*   **`init_hamiltonian` (Character string from `m_inputparam`)**: Controls the method for generating the initial guess for MO coefficients.

## Usage Examples

This module is primarily used internally by the SCF (Self-Consistent Field) procedure in MOLGW.
1.  **Initialization**: `init_c_matrix` is called at the beginning of an SCF calculation to get an initial set of MO coefficients.
2.  **SCF Iteration**: Inside an SCF loop:
    *   `setup_density_matrix` (from `m_hamiltonian_tools`) is called to compute `p_matrix` from current MOs.
    *   `setup_hartree` or `setup_hartree_ri` is called to compute the Hartree potential.
    *   If Hartree-Fock or hybrid DFT: `setup_exchange` or `setup_exchange_ri` is called for the exchange potential.
    *   If DFT: `dft_exc_vxc_batch` is called for the XC potential.
    *   These potentials are combined with one-body terms (kinetic, nuclear, ECPs from `m_hamiltonian_onebody`) to form the new Fock/Kohn-Sham matrix.
    *   This new matrix is diagonalized (e.g., by `m_hamiltonian_tools::diagonalize_hamiltonian_scalapack`) to get updated MOs and energies.

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_timing`, `m_mpi`, `m_scalapack`, `m_warning`, `m_memory`, `m_cart_to_pure`, `m_inputparam`.
*   **Basis Set & ERI**: `m_basis_set`, `m_eri` (for AO integrals and indexing), `m_eri_calculate` (for `eri_2center_inv`).
*   **Density & Grid**: `m_density_tools` (for `calc_density_r_batch`, etc.), `m_dft_grid` (provides the grid and basis function values/gradients on it for `dft_exc_vxc_batch`).
*   **XC Library**: `m_libxc_tools` (provides the interface to Libxc for evaluating XC functionals in `dft_exc_vxc_batch`).
*   **Hamiltonian Tools**: `m_hamiltonian_tools` (for `diagonalize_hamiltonian_scalapack`, `setup_density_matrix`).
*   **I/O**: `m_io` (for reading Gaussian fchk files in `init_c_matrix`).
*   This module is central to constructing the Fock or Kohn-Sham matrix, which is the heart of SCF calculations. It relies heavily on efficient ERI handling (from `m_eri`) and numerical integration (from `m_dft_grid`).
```
