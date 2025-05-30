# `m_linear_response.f90`

## Overview

The `m_linear_response` Fortran module in MOLGW is responsible for calculating the frequency-dependent polarizability &chi;(&omega;) and related quantities. This is achieved by solving the Casida equations (for Time-Dependent Density Functional Theory - TDDFT, or Time-Dependent Hartree-Fock - TDHF) or the Bethe-Salpeter Equation (BSE). The computed polarizability or screened Coulomb interaction W(&omega;) is fundamental for obtaining optical excitation spectra, and W(&omega;) is a key ingredient for GW self-energy calculations. The module handles different approximations (RPA, TDDFT, BSE), spin multiplicities (singlet, triplet), and can incorporate a coupling constant for Adiabatic Connection Fluctuation-Dissipation Theorem (ACFD) calculations.

## Key Components

*   **`polarizability(...)`**:
    *   **Purpose**: This is the main driver subroutine for calculating the dynamic polarizability and/or the screened Coulomb interaction W(&omega;).
    *   **Inputs**: Takes SCF results (`occupation`, `energy`, `c_matrix`), basis set information, and control flags like `enforce_rpa`, `calculate_w`, `enforce_spin_multiplicity`, and `lambda`.
    *   **Workflow**:
        1.  Determines the type of calculation (RPA, TDDFT, BSE, TDHF) based on input flags and global parameters.
        2.  Constructs the (A-B) and (A+B) matrices (or A and B separately) by calling routines in `m_build_bse`. `A` contains orbital energy differences and two-electron integrals (Coulomb and exchange-type terms), while `B` contains only two-electron integrals.
        3.  Solves the Casida/BSE eigenvalue problem `[(A-B)(A+B)]^(1/2) Z = Z &Omega;` or `(A-B)X = &Omega;X`, `(A+B)Y = &Omega;Y` by calling diagonalization routines from `m_block_diago`. This yields excitation energies (`eigenvalue`) and corresponding eigenvectors (`xpy_matrix` for X+Y, `xmy_matrix` for X-Y).
        4.  If requested (`calculate_w` or `print_w_`), transforms these eigenvalues and eigenvectors into a spectral representation of the polarizability &chi; or the screened Coulomb interaction W, storing the results (poles and residues) in the `wpol_out` object of type `spectral_function`. This transformation is done by `chi_to_sqrtvchisqrtv_auxil` (for RI) or `chi_to_vchiv` (for 4-center integrals).
        5.  Optionally calculates and prints optical spectra and stopping power via routines in `m_spectra`.
    *   **Outputs**: `en_rpa` (RPA correlation energy part), `en_gw` (GW correlation energy part if RI is used), and `wpol_out` (spectral function for &chi; or W). Optional outputs include the A, B, X, and Y matrices themselves.

*   **`coupled_perturbed(...)`**:
    *   **Purpose**: Calculates the static inverse of the (A+B) matrix, i.e., (A+B)<sup>-1</sup>. This is relevant for static Coupled-Perturbed Hartree-Fock/Kohn-Sham (CPHF/CPKS) theory, used for static polarizabilities or other linear response properties at zero frequency.
    *   It calls `polarizability` to construct the A and B matrices and then inverts their sum.

*   **`polarizability_onering(...)`**:
    *   **Purpose**: Calculates the "one-ring" approximation to the polarizability. This often refers to the non-interacting response function &chi;<sub>0</sub>, but dressed with v<sup>1/2</sup> factors if using an auxiliary basis for RI: (v<sup>1/2</sup>&chi;<sub>0</sub>v<sup>1/2</sup>). The poles are orbital energy differences, and residues are products of 3-center integrals.
    *   Stores results in the `vchi0v` spectral function object.

*   **`get_energy_qp(energy, occupation, energy_qp)`**:
    *   **Purpose**: A utility to obtain quasiparticle (QP) energies for use in BSE calculations.
    *   If `scissor` (input parameter) is non-zero, it applies a scissor correction to the input `energy` (typically Kohn-Sham energies).
    *   Otherwise, it attempts to read QP energies from an `ENERGY_QP` file.
    *   If neither is available, it defaults to using the input `energy` as QP energies.

*   **`chi_to_vchiv(c_matrix, xpy_matrix, eigenvalue, wpol)`**:
    *   **Purpose**: Converts the Casida eigenvalues (&Omega;<sub>I</sub>) and eigenvectors (X<sub>I</sub>+Y<sub>I</sub>) into the spectral representation of the screened Coulomb interaction W = v + v&chi;v, when using the 4-center ERI formalism.
    *   The residues stored in `wpol%residue_left` are of the form (&mu;&nu;|XY)<sub>I</sub> = &sum;<sub>ia,jb</sub> (&mu;&nu;|ij)(X<sub>ia</sub>+Y<sub>ia</sub>)<sub>I</sub>, where (ij) are MO pairs. (Note: This is a simplification; the actual form is more complex, involving contractions with four MOs).

*   **`chi_to_sqrtvchisqrtv_auxil(desc_x, xpy_matrix, eigenvalue, wpol, energy_gm)`**:
    *   **Purpose**: Converts Casida eigenvalues and eigenvectors into the spectral representation v<sup>1/2</sup>&chi;v<sup>1/2</sup>, when using the Resolution of Identity (RI) with an auxiliary basis.
    *   The residues stored in `wpol%residue_left` are of the form (Q|XY)<sub>I</sub> = &sum;<sub>ia</sub> (Q|ia)(X<sub>ia</sub>+Y<sub>ia</sub>)<sub>I</sub>, where Q is an auxiliary basis function.
    *   Also calculates the Galitskii-Migdal correlation energy (`energy_gm`) if RI is used.

*   **`static_polarizability(occupation, energy, wpol_out)`**:
    *   **Purpose**: Calculates the static RPA polarizability &chi;(&omega;=0) using the RI formalism. It constructs and inverts (1 - v<sup>1/2</sup>&chi;<sub>0</sub>v<sup>1/2</sup>).

## Important Variables/Constants

*   **`wpol_out` (Type `spectral_function`)**: The primary output for routines calculating dynamic response, storing poles and residues of &chi;(&omega;) or W(&omega;).
*   **`enforce_rpa` (Logical, Intent(in))**: If `.TRUE.`, forces the calculation to use the RPA kernel (bare Coulomb, no exchange or XC kernel), overriding other settings.
*   **`calculate_w` (Logical, Intent(in))**: If `.TRUE.`, the routine explicitly computes and stores the screened Coulomb interaction W in `wpol_out`. Otherwise, `wpol_out` might contain &chi; or be partially filled.
*   **`enforce_spin_multiplicity` (Integer, Intent(in), Optional)**: Allows forcing singlet (1) or triplet (3) calculations, overriding the global `is_triplet` flag.
*   **`lambda` (Real(dp), Intent(in), Optional)**: Adiabatic connection parameter, scales interaction terms. Defaults to 1.0.
*   **`xpy_matrix`, `xmy_matrix` (Real(dp))**: Arrays storing the (X+Y) and (X-Y) solution vectors of the Casida/BSE equations.
*   **`eigenvalue` (Real(dp))**: Array storing the excitation energies (&Omega;).
*   **`is_bse`, `is_tdhf`, `is_tddft`, `is_rpa` (Logicals)**: Internal flags set based on input parameters to control the specific type of linear response calculation.

## Usage Examples

This module is a core part of MOLGW's post-SCF capabilities. It is invoked when calculations of excitation spectra (TDDFT/BSE) or GW self-energies are requested.

```fortran
! Conceptual call from a GW driver routine:
USE m_linear_response
USE m_spectral_function
! ... (SCF results: occupation_scf, energy_scf, c_matrix_scf, orbital_basis are available) ...
! ... (en_rpa_corr, en_gw_corr are scalar energy outputs) ...

TYPE(spectral_function) :: w_screener
CALL w_screener%init(n_mo_total, occupation_scf, num_freq_points_for_W, grid_type=POLE_SUM)

! Calculate W(omega) within RPA and store it in w_screener
CALL polarizability(enforce_rpa=.TRUE., calculate_w=.TRUE., &
                    basis=orbital_basis, occupation=occupation_scf, &
                    energy=energy_scf, c_matrix=c_matrix_scf, &
                    en_rpa=en_rpa_corr, en_gw=en_gw_corr, wpol_out=w_screener)

! w_screener now contains the poles and residues of the RPA W(omega),
! which can be used in m_gw_selfenergy_analytic or m_gw_selfenergy_grid.
! ...
CALL w_screener%destroy()
```

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_timing`, `m_warning`, `m_memory`, `m_inputparam`, `m_mpi`, `m_scalapack`.
*   **Basis Set & ERI**: `m_basis_set`, `m_eri_ao_mo` (for `calculate_eri_3center_eigen` if RI is used).
*   **Hamiltonian Construction**: `m_build_bse` (critical dependency for constructing the A and B matrices that define the linear response eigenvalue problem).
*   **Diagonalization**: `m_block_diago` (for solving the Casida/BSE eigenvalue problem via `diago_4blocks_chol`, `diago_4blocks_rpa_sca`, or `diago_4blocks_davidson`).
*   **Spectral Objects**: `m_spectral_function` (defines the `wpol_out` type), `m_spectra` (for `optical_spectrum`, `stopping_power` calculations using the results).
*   The module orchestrates the complex process of solving linear response equations, from matrix construction to diagonalization and the formation of the spectral function for &chi; or W.
```
