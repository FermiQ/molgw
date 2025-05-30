# `m_gw_selfenergy_analytic.f90`

## Overview

The `m_gw_selfenergy_analytic` Fortran module in MOLGW is dedicated to calculating the electron self-energy within various analytic formulations of the GW approximation and related methods. This includes standard G0W0, eigenvalue self-consistent GnW0, fully self-consistent GnWn, COHSEX (Coulomb Hole and Screened Exchange approximation), and Quasiparticle Self-Consistent GW (QSGW). The calculations primarily rely on the spectral representation (poles and residues) of the screened Coulomb interaction W(&omega;) and Green's function G.

## Key Components

*   **`gw_selfenergy(selfenergy_approx, occupation, energy, c_matrix, wpol, se)`**:
    *   **Purpose**: This is the primary subroutine for computing the GW self-energy &Sigma;(&omega;).
    *   `selfenergy_approx`: An integer flag specifying the GW flavor (GW for G0W0, GnW0, GnWn, COHSEX, ONE_RING).
    *   It takes molecular orbital `occupation`, `energy`, coefficients `c_matrix`, and the `spectral_function` object `wpol` (containing W(&omega;)) as input.
    *   The computed self-energy is stored in the `se%sigma` array of the `selfenergy_grid` object `se`.
    *   The calculation involves summations over occupied and virtual states and the poles of W(&omega;), using the formula &Sigma; = iGW.
    *   For COHSEX, it computes the frequency-independent self-energy.

*   **`gw_selfenergy_upfolding(selfenergy_approx, occupation, energy, c_matrix, wpol, exchange_m_vxc)`**:
    *   **Purpose**: Implements a method to obtain quasiparticle energies by solving the Dyson equation via "upfolding" the frequency-dependent self-energy into a larger, frequency-independent eigenvalue problem.
    *   This creates an effective Hamiltonian matrix whose dimensions are (number of MOs + number of MOs * number of poles in W).
    *   The diagonalization of this large matrix yields the quasiparticle energies. This approach is typically used for G0W0.
    *   `exchange_m_vxc` is the <&phi;<sub>p</sub>| &Sigma;<sub>x</sub> - V<sub>xc</sub> |&phi;<sub>q</sub>> matrix.

*   **`gw_selfenergy_scalapack(selfenergy_approx, occupation, energy, c_matrix, wpol, se)`**:
    *   **Purpose**: A ScaLAPACK-enabled version of `gw_selfenergy` for parallel computation of the GW self-energy.
    *   It distributes the calculation of self-energy contributions over MPI processes, particularly the summation over poles of W and intermediate states.
    *   Requires `has_auxil_basis` to be true.

*   **`gw_selfenergy_qs(occupation, energy, c_matrix, s_matrix, wpol, selfenergy)`**:
    *   **Purpose**: Calculates the quasiparticle self-energy matrix elements <&phi;<sub>p</sub>| &Sigma;<sub>c</sub>(&omega;=(&epsilon;<sub>p</sub>+&epsilon;<sub>q</sub>)/2) |&phi;<sub>q</sub>> for Quasiparticle Self-Consistent GW (QSGW).
    *   It computes the static (&omega; averaged or set to zero, depending on implementation details not fully clear from signature alone) correlation part of the self-energy.
    *   Applies Kotani's hermitianization trick to ensure the self-energy matrix is Hermitian: &Sigma;<sub>pq</sub><sup>QSGW</sup> = 1/2 (&Sigma;<sub>pq</sub> + &Sigma;<sub>qp</sub><sup>*</sup>).
    *   The output `selfenergy` is the static, hermitianized self-energy matrix in the MO basis.

*   **`dump_gw_ingredients(energy, c_matrix, wpol)`**:
    *   **Purpose**: A utility subroutine to write various components used in or produced by the GW calculation to files for debugging or analysis.
    *   Outputs include Kohn-Sham energies (`energies.dat`), RPA neutral excitation energies/poles from `wpol` (`omegas.dat`), 3-center MO integrals (v coefficients, `v.bin`), and W coefficients (`w.bin`).

## Important Variables/Constants

*   **`selfenergy_approx` (Integer)**: Input flag controlling the GW method (e.g., `GW`, `GnW0`, `COHSEX`).
*   **`wpol` (Type `spectral_function`)**: Contains the spectral representation of the screened Coulomb interaction W(&omega;), including its poles (`wpol%pole`) and residues (`wpol%residue_left`). This is a crucial input.
*   **`se` (Type `selfenergy_grid`)**: An object that stores the Kohn-Sham energies (`se%energy0`) on which the self-energy is evaluated, the frequency grid (`se%omega`), and holds the computed self-energy `se%sigma(&omega;, state, spin)`.
*   **`eri_3center_eigen` (from `m_eri_ao_mo`)**: Array of 3-center ERIs transformed to the MO basis, of the form (AuxBasisFunc | MO<sub>i</sub> MO<sub>j</sub>). Used extensively if `has_auxil_basis` is true.
*   **`index_prodstate` (Function from `m_selfenergy_tools`)**: Maps a pair of MO indices (p,q) to a linear index, used for accessing `wpol%residue_left` when it's stored in a product basis.
*   **`ieta` (Real(dp) from `m_definitions`)**: A small positive infinitesimal used for Green's function denominators.

## Usage Examples

These routines are typically called from a higher-level GW driver within MOLGW after an SCF calculation and the calculation of W(&omega;) (which populates `wpol`).

```fortran
! Conceptual: Inside a GW workflow
USE m_gw_selfenergy_analytic
USE m_spectral_function
USE m_selfenergy_tools
! ... (SCF results: occupation_scf, energy_scf, c_matrix_scf) ...
! ... (w_response is a populated spectral_function object for W(omega)) ...
! ... (se_gw is an initialized selfenergy_grid object) ...

INTEGER :: gw_approximation_type ! Set based on user input, e.g., to GW for G0W0

CALL gw_selfenergy(gw_approximation_type, occupation_scf, energy_scf, c_matrix_scf, &
                   w_response, se_gw)

! se_gw%sigma now contains the G0W0 self-energy.
! Quasiparticle energies can then be found by solving Dyson's equation using se_gw.
```
For QSGW:
```fortran
REAL(DP), ALLOCATABLE :: qs_sigma_matrix(:,:,:)
! ... (Populate qs_sigma_matrix usually with Sigma_x - Vxc as starting point) ...
CALL gw_selfenergy_qs(occupation_scf, energy_scf, c_matrix_scf, overlap_matrix_ao, &
                      w_response, qs_sigma_matrix)
! qs_sigma_matrix now contains the static hermitianized QSGW self-energy.
```

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_mpi`, `m_timing`, `m_inputparam`, `m_warning`.
*   **Basis Set & ERI**: `m_basis_set`, `m_eri_ao_mo` (for `calculate_eri_3center_eigen` and `eri_3center_eigen` which provide MO integrals if `has_auxil_basis` is true).
*   **Spectral Function & Self-Energy**: `m_spectral_function` (defines `wpol`), `m_selfenergy_tools` (defines `se` type and `index_prodstate`).
*   **Linear Algebra**: `m_linear_algebra` (for `diagonalize` if used in `gw_selfenergy_upfolding`), ScaLAPACK routines are used in `gw_selfenergy_scalapack`.
*   **Input Parameters**: The choice of GW method (`selfenergy_approx`), number of states (`nsemin`, `nsemax`), etc., are controlled by parameters from `m_inputparam`.
*   This module is central to performing GW calculations. It takes the pre-calculated screened Coulomb interaction W and combines it with Green's functions G (constructed from MO energies and occupations) to compute the self-energy &Sigma;.
```
