# `m_g3w2_selfenergy.f90`

## Overview

The `m_g3w2_selfenergy` Fortran module in MOLGW is dedicated to calculating higher-order corrections to the electron self-energy, going beyond the standard GW approximation. It implements several advanced diagrammatic approaches, including Second-Order Screened Exchange (SOSEX), its variant 2SOSEX, and the G3W2 self-energy (a third-order diagram involving two screened Coulomb interactions W and three Green's functions G). These methods aim to provide more accurate quasiparticle energies by incorporating more sophisticated electron correlation effects.

## Key Components

*   **`sosex_selfenergy(basis, occupation, energy, c_matrix, wpol, se)`**:
    *   **Purpose**: Calculates the GW self-energy augmented with the SOSEX correction. The SOSEX term can be schematically represented as G*W*G*v*G, where 'v' is the bare Coulomb interaction and 'W' is the dynamically screened Coulomb interaction.
    *   If `gwgamma_tddft_` is enabled, it can also include a TDDFT exchange-correlation kernel contribution to the vertex function within the SOSEX diagrams.
    *   The `factor_sosex` input parameter can scale the SOSEX contribution (e.g., for 2SOSEX).

*   **`sosex_selfenergy_analyzed(basis, occupation, energy, c_matrix, wpol, se)`**:
    *   **Purpose**: A more detailed version of `sosex_selfenergy` designed for analysis. It calculates and explicitly stores various components of the SOSEX self-energy, such as ring GPP diagram contributions and SOX-like terms, categorizing them by the nature of intermediate states (e.g., hole-plasmon, particle-plasmon).
    *   Writes detailed weight analyses of different diagrammatic contributions to output files.

*   **`gwgw0g_selfenergy(occupation, energy, c_matrix, wpol, se)`**:
    *   **Purpose**: Computes self-energy corrections involving terms with both dynamic (W(&omega;)) and static (W(&omega;=0) or W<sup>RPA</sup>(&omega;=0)) screened Coulomb interactions.
    *   Specifically, it calculates terms like G*W*G*W<sub>0</sub>*G. Depending on `calc_type%selfenergy_approx`, W<sub>0</sub> can be the static limit of the full W used in GW, or the static RPA screened interaction.
    *   Implemented methods: `GW0GW0G` (G*W<sub>0</sub>*G*W<sub>0</sub>*G), `GWGW0G` (2*G*v*G*W<sub>0</sub>*G - G*W<sub>0</sub>*G*W<sub>0</sub>*G + 2*G*(W-v)*G*W<sub>0</sub>*G), and `GWGW0RPAG` (G*v*G*W<sub>0</sub><sup>RPA</sup>*G + G*(W-v)*G*W<sub>0</sub><sup>RPA</sup>*G).

*   **`g3w2_selfenergy(occupation, energy, c_matrix, wpol, se)`**:
    *   **Purpose**: Implements the G3W2 self-energy, which includes third-order diagrams with two screened interactions W and three Green's functions G (schematically G*W*G*W*G).
    *   It sums over various time-orderings of these interactions, which correspond to different energy denominators in the perturbation theory expressions.
    *   Input flags `g3w2_skip_vvv_` and `g3w2_skip_vv_` allow neglecting computationally expensive diagram subsets involving three virtual lines or two virtual lines, respectively.

*   **`g3w2_selfenergy_real_grid(basis, occupation, energy, c_matrix, se)`**:
    *   **Purpose**: Calculates the G3W2 self-energy by performing a direct numerical double integration over real frequency axes for the two W lines. This is an alternative to the pole-based summation used in `g3w2_selfenergy`.

*   **`g3w2_selfenergy_imag_grid(basis, occupation, energy, c_matrix, se)`**:
    *   **Purpose**: Calculates the G3W2 self-energy by performing a numerical double integration over imaginary frequency axes for the two W lines. This method can offer better numerical stability.

*   **`sox_selfenergy_imag_grid(occupation, energy, c_matrix, se)`**:
    *   **Purpose**: Calculates the Second-Order Exchange (SOX) self-energy (G*v*G*v*G) on an imaginary frequency grid and adds it to an existing self-energy (presumably a GW self-energy) stored in the `se` object.

## Important Variables/Constants

*   **`wpol` (Type `spectral_function`)**: Input object containing the spectral representation (poles and residues) or frequency-dependent values of the screened Coulomb interaction W(&omega;).
*   **`se` (Type `selfenergy_grid`)**: Input/Output object that stores the self-energy &Sigma;(E) on a defined energy/frequency grid. The corrections computed by this module are added to `se%sigma`.
*   **`calc_type%selfenergy_approx` (from `m_inputparam`)**: A string that determines which specific self-energy approximation (e.g., 'GWSOX', 'GWSOSEX', 'G3W2', 'GWGW0G') is to be computed.
*   **`factor_sosex` (from `m_inputparam`)**: A scaling factor for the SOSEX contribution (e.g., 1.0 for SOSEX, 2.0 for 2SOSEX).
*   **`gwgamma_tddft_` (Logical, from `m_inputparam`)**: If `.TRUE.`, a TDDFT-derived exchange-correlation kernel is included in the vertex function for SOSEX calculations.
*   **`alpha_hybrid` (Real(dp), from `m_inputparam`)**: Fraction of exact exchange to be used if `gwgamma_tddft_` is active.
*   **`eri_eigen(...)` (Function from `m_eri_ao_mo`)**: Used to retrieve two-electron repulsion integrals in the Molecular Orbital (MO) basis.

## Usage Examples

These subroutines are advanced post-GW methods and are called internally within MOLGW's workflow when specified by the user. The typical sequence involves:
1.  An initial SCF (DFT or HF) calculation.
2.  A GW calculation to obtain an initial `wpol` and `se%sigma` (the GW self-energy).
3.  Then, one of the subroutines from `m_g3w2_selfenergy` is called to compute higher-order corrections, which are added to `se%sigma`.

Example MOLGW input snippet:
```
task = GW
...
selfenergy_approx = G3W2 ! Request G3W2 calculation after initial G0W0
```
MOLGW's GW driver would first compute G0W0, populate `wpol` and `se`, and then call `g3w2_selfenergy` (or `sosex_selfenergy` if `selfenergy_approx = GWSOSEX`).

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_mpi`, `m_timing`, `m_inputparam`, `m_warning`, `m_basis_set`.
*   **ERI and Transformations**: `m_eri_ao_mo` (for MO-basis ERIs, e.g., `eri_eigen` and `calculate_eri_3center_eigen`).
*   **Self-Energy & Spectral Functions**: `m_selfenergy_tools` (defines `selfenergy_grid` type), `m_spectral_function` (defines `spectral_function` type `wpol`).
*   **Response Theory**: `m_linear_response` (may be used to compute `wpol` if `analytic_chi_` is true).
*   **TDDFT Components**: `m_tddft_fxc` (if `gwgamma_tddft_` is true, for evaluating fxc kernels).
*   **Input Parameters**: The behavior is heavily controlled by `calc_type%selfenergy_approx`, `factor_sosex`, `g3w2_skip_vvv_`, `g3w2_skip_vv_`, etc., read via `m_inputparam`.
*   The module implements complex many-body perturbation theory formulas involving sums over occupied and virtual states and frequency/pole summations. The computational cost generally increases with the sophistication of the method (SOX < SOSEX < G3W2).
```
