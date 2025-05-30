# `m_gw_selfenergy_grid.f90`

## Overview

The `m_gw_selfenergy_grid` Fortran module in MOLGW is dedicated to the calculation of the GW self-energy (&Sigma;) and the Random Phase Approximation (RPA) polarizability (&chi;) using numerical integration over frequency grids. This typically involves calculations on the imaginary frequency axis, which can offer numerical stability and facilitate the use of contour deformation techniques for obtaining real-axis quantities. The module provides routines for constructing &chi;(&omega;) and subsequently computing &Sigma;(&omega;) = iG(&omega;)W(&omega;), where W(&omega;) = v + v&chi;(&omega;)v.

## Key Components

*   **`polarizability_grid_scalapack(occupation, energy, c_matrix, erpa, egw, wpol)`**:
    *   **Purpose**: Calculates the RPA polarizability &chi;(&omega;) on a predefined imaginary frequency grid stored in `wpol`.
    *   It first computes the non-interacting polarizability &chi;<sub>0</sub>(&omega;) using molecular orbital energies and 3-center ERI_MOs (from `m_eri_ao_mo`).
    *   Then, it solves the Dyson-like equation &chi; = &chi;<sub>0</sub> + &chi;<sub>0</sub> v &chi;  (or equivalently &chi; = (1 - &chi;<sub>0</sub> v)<sup>-1</sup>&chi;<sub>0</sub>, though the code implements (v<sup>1/2</sup>&chi;<sub>0</sub>v<sup>1/2</sup>) type matrices). The result &chi;(&omega;) is stored in `wpol%chi`.
    *   It also computes the RPA correlation energy (`erpa`) and a GW correlation energy estimate (`egw`) by integrating expressions involving &chi;(&omega;) over the imaginary frequency grid.
    *   Uses ScaLAPACK for distributed matrix operations.

*   **`gw_selfenergy_imag_scalapack(energy, c_matrix, wpol, se)`**:
    *   **Purpose**: Computes the GW self-energy &Sigma;(&omega;) on a predefined imaginary frequency grid specified in the `se` object (`se%omega_calc`).
    *   It uses the formula &Sigma;(i&omega;<sub>n</sub>) = - &int; d(i&omega;')/(2&pi;) G(i&omega;<sub>n</sub>+i&omega;')W(i&omega;').
    *   `wpol%chi` (pre-calculated by `polarizability_grid_scalapack` or similar) provides W(i&omega;').
    *   The Green's function G is constructed from MO energies.
    *   The integration over imaginary frequencies i&omega;' is performed numerically using the quadrature weights in `wpol`.
    *   The resulting self-energy is stored in `se%sigma_calc`.

*   **`gw_selfenergy_contour(energy, occupation, c_matrix, se)`**:
    *   **Purpose**: Calculates the GW self-energy on a real frequency grid using contour deformation techniques.
    *   This method combines integration along the imaginary axis (using `wpol_imag` which it computes via RPA) with contributions from residues of poles near the real axis (using `wpol_real`, also computed via RPA on a real frequency grid).
    *   This approach is more complex but aims to accurately capture the structure of the self-energy on the real frequency axis.

*   **`gw_selfenergy_grid(basis, occupation, energy, c_matrix, se)`**:
    *   **Purpose**: Calculates the GW self-energy on an imaginary frequency grid. Similar to `gw_selfenergy_imag_scalapack`.
    *   It can either use an analytically constructed &chi;<sub>0</sub> (if `analytic_chi_` is true, via `m_linear_response`) or compute &chi;<sub>0</sub> from sums over states to build W(&omega;) stored in `wpol_imag%chi`.
    *   The self-energy is then computed by numerical integration over the imaginary frequency grid.

*   **`fsos_selfenergy_grid(basis, energy, occupation, c_matrix, se)`**:
    *   **Purpose**: Calculates a "Full Second Order Screened" (FSOS) self-energy on an imaginary frequency grid. This is likely a more advanced self-energy diagram beyond standard GW, possibly related to SOSEX or G3W2-type terms but evaluated fully on the grid.
    *   It involves two screened interactions W, implying a higher-order expansion.
    *   Boolean parameters `impose_sox` and `impose_sosex` suggest it can be reduced to simpler forms.

## Important Variables/Constants

*   **`wpol` (Type `spectral_function`)**: An object that stores the frequency-dependent polarizability &chi;(i&omega;) or screened Coulomb interaction W(i&omega;) on a grid of imaginary frequencies. It includes `wpol%omega` (grid points), `wpol%weight_quad` (quadrature weights), and `wpol%chi` (values of &chi; or related quantities).
*   **`se` (Type `selfenergy_grid`)**: An object that defines the energy/frequency grid for the self-energy (`se%omega_calc`) and stores the computed self-energy `se%sigma_calc`. It also contains the initial Kohn-Sham energies `se%energy0`.
*   **`has_auxil_basis` (Logical, from `m_inputparam`)**: If true, calculations use 3-center ERIs with an auxiliary basis (RI approximation). This is a requirement for most routines in this module.
*   **`analytic_chi_` (Logical, from `m_inputparam`)**: If true, the non-interacting polarizability &chi;<sub>0</sub> is constructed using analytical expressions (typically from `m_linear_response`). Otherwise, it's computed by direct summation over occupied-virtual pairs.
*   **`nomega_chi_imag` (Integer, from `m_inputparam`)**: Number of imaginary frequency points used for the quadrature of W and &Sigma;.
*   **`eri_3center_eigen` (from `m_eri_ao_mo`)**: Array of 3-center ERIs transformed to the MO basis, (AuxBasisFunc | MO<sub>i</sub> MO<sub>j</sub>). Essential for constructing &chi;<sub>0</sub> and W in the RI approximation.

## Usage Examples

These routines are internal to MOLGW's GW calculation workflow. They are typically called after an SCF calculation and after the necessary 3-center MO integrals have been prepared.

```fortran
! Conceptual flow within a GW calculation:
USE m_gw_selfenergy_grid
USE m_spectral_function
USE m_selfenergy_tools
! ... (SCF results: occupation_scf, energy_scf, c_matrix_scf) ...
! ... (eri_3center_eigen are computed and available via m_eri_ao_mo) ...

TYPE(spectral_function) :: w_chi_on_grid
TYPE(selfenergy_grid) :: sigma_on_grid

! Initialize w_chi_on_grid with desired imaginary frequency grid
CALL w_chi_on_grid%init(n_mo_total, occupation_scf, num_imag_freq_points, grid_type=IMAGINARY_QUAD)

! Calculate Chi(iw) and store in w_chi_on_grid%chi
CALL polarizability_grid_scalapack(occupation_scf, energy_scf, c_matrix_scf, &
                                   rpa_corr_energy, gw_corr_energy, w_chi_on_grid)

! Initialize sigma_on_grid with the target real or imaginary frequency grid for Sigma
CALL sigma_on_grid%init(use_contour_deformation, kohn_sham_energies_for_sigma_eval)

! Calculate Sigma(iw) and store in sigma_on_grid%sigma_calc
CALL gw_selfenergy_imag_scalapack(energy_scf, c_matrix_scf, w_chi_on_grid, sigma_on_grid)
! OR
! CALL gw_selfenergy_contour(energy_scf, occupation_scf, c_matrix_scf, sigma_on_grid)


! ... (Then use sigma_on_grid to find quasiparticle energies) ...

CALL w_chi_on_grid%destroy()
CALL sigma_on_grid%destroy()
```

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_mpi`, `m_timing`, `m_inputparam`, `m_warning`, `m_memory`, `m_scalapack`, `m_linear_algebra`.
*   **Basis Set & ERI**: `m_basis_set`, `m_eri_ao_mo` (for `calculate_eri_3center_eigen` and accessing `eri_3center_eigen`).
*   **Spectral Function & Self-Energy**: `m_spectral_function` (defines `wpol`), `m_selfenergy_tools` (defines `se`).
*   **Response Theory**: `m_linear_response` (potentially for `wpol_analytic%evaluate` if `analytic_chi_` is true).
*   **ScaLAPACK/BLAS**: Heavily used for distributed matrix algebra (e.g., `PDSYRK`, `PDGEMM`, `PDGEMR2D`) in the `_scalapack` versions of routines, and `DSYRK`, `DGEMM` in serial or OpenMP sections.
*   This module provides an alternative route to the GW self-energy compared to `m_gw_selfenergy_analytic`, focusing on numerical integration on frequency grids. This can be more robust for systems where pole-based summations become difficult.
```
