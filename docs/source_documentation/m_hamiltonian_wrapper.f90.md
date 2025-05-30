# `m_hamiltonian_wrapper.f90`

## Overview

The `m_hamiltonian_wrapper` Fortran module in MOLGW serves as a high-level interface for constructing the two-body components of the Kohn-Sham or Fock Hamiltonian. It wraps the more detailed routines found in `m_hamiltonian_twobodies`, providing a simplified API for calculating the Hartree potential, exchange potential (including long-range versions), and the DFT exchange-correlation (XC) potential. The wrapper functions intelligently select the appropriate underlying computational strategy based on global input parameters, such as whether an auxiliary basis is used (for Resolution of Identity, RI), whether calculations are in-core, or if "genuine" 3-center integrals are employed in RI.

## Key Components

*   **`calculate_hartree(basis, p_matrix, hhartree, eh)`**:
    *   **Purpose**: Calculates the Hartree potential matrix and, optionally, the Hartree energy.
    *   **Logic**:
        *   If not using an auxiliary basis (`has_auxil_basis` is false):
            *   If `incore_` is true (integrals fit in memory), it calls `m_hamiltonian_twobodies::setup_hartree` (full 4-center ERI-based).
            *   Else (out-of-core), it calls `m_hamiltonian_twobodies::setup_hartree_oneshell` (shell-wise 4-center ERI processing).
        *   If using an auxiliary basis:
            *   If `eri3_genuine_` is false (standard RI), it calls `m_hamiltonian_twobodies::setup_hartree_ri`.
            *   Else (genuine 3-center RI), it first calls `m_hamiltonian_twobodies::calculate_density_auxilbasis` to get density coefficients in the auxiliary basis, then calls `m_hamiltonian_twobodies::setup_hartree_genuine_ri`.
    *   `p_matrix` can be real or complex (via `CLASS(*)`).

*   **`calculate_exchange` (Interface for `calculate_exchange_real`)**:
    *   **`calculate_exchange_real(basis, p_matrix, hexx, ex, occupation, c_matrix)`**:
        *   **Purpose**: Calculates the real-valued exchange potential matrix and, optionally, the exchange energy.
        *   **Logic**: Similar to `calculate_hartree`, it chooses between:
            *   `m_hamiltonian_twobodies::setup_exchange` (4-center, in-core).
            *   Issues a warning if out-of-core 4-center exchange is requested (not implemented).
            *   `m_hamiltonian_twobodies::setup_exchange_ri` or `m_hamiltonian_twobodies::setup_exchange_genuine_ri` if using an auxiliary basis, depending on `eri3_genuine_`. If MO coefficients (`c_matrix`, `occupation`) are not provided, they are reconstructed from `p_matrix`.

*   **`calculate_exchange_lr(basis, p_matrix, hexx, ex, occupation, c_matrix)`**:
    *   **Purpose**: Calculates the long-range part of the exchange potential and energy, typically for range-separated hybrid functionals.
    *   **Logic**: Similar to `calculate_exchange_real`, but calls the `_longrange` versions of the underlying routines in `m_hamiltonian_twobodies` (e.g., `setup_exchange_longrange`, `setup_exchange_longrange_ri`).

*   **`calculate_hamiltonian_hxc(basis, nstate, occupation, c_matrix, p_matrix, hamiltonian_hxc, en_inout)`**:
    *   **Purpose**: Constructs the total Hartree + Exchange + XC (HXC) part of the Kohn-Sham/Fock matrix and updates energy components.
    *   **Logic**:
        1.  Calls `calculate_hartree` to get the Hartree potential and energy.
        2.  If `calc_type%is_dft` is true, calls `m_hamiltonian_twobodies::dft_exc_vxc_batch` to get the DFT XC potential and energy.
        3.  If `calc_type%need_exchange_lr` is true, calls `calculate_exchange_lr` and adds the scaled (by `beta_hybrid`) long-range exchange potential and energy.
        4.  If `calc_type%need_exchange` is true, calls `calculate_exchange` (for the short-range or full exchange part) and adds the scaled (by `alpha_hybrid`) exchange potential and energy.
        *   Accumulates energy contributions into `en_inout`.

*   **`calculate_hamiltonian_hxc_ri_cmplx(...)`**:
    *   **Purpose**: A version of `calculate_hamiltonian_hxc` specifically for complex arithmetic and RI methods.
    *   Calls complex-valued RI exchange setup routines from `m_hamiltonian_twobodies` (e.g., `setup_exchange_ri_cmplx`).
    *   Hartree and DFT XC parts are still fundamentally real but are added to the complex HXC matrix.

*   **`calculate_hamiltonian_hartree_x2c(...)`**:
    *   **Purpose**: Calculates only the Hartree potential component for x2c (exact two-component relativistic) calculations. The Hartree potential is real even if the density matrix is complex.

*   **`calculate_hamiltonian_xc_x2c(...)`**:
    *   **Purpose**: Calculates only the DFT XC potential component for x2c calculations. The XC potential is real.

## Important Variables/Constants

*   **Input Parameters (from `m_inputparam`)**:
    *   `has_auxil_basis` (Logical): Determines if RI methods are used.
    *   `incore_` (Logical): If `has_auxil_basis` is false, determines if 4-center ERIs are stored in memory or processed shell-wise.
    *   `eri3_genuine_` (Logical): If `has_auxil_basis` is true, selects between standard RI (V<sup>-1/2</sup> approach) or using genuine 3-center integrals with a density fit.
*   **SCF Control Flags (from `m_scf::calc_type`)**:
    *   `calc_type%is_dft` (Logical): True if a DFT calculation is being performed (triggers XC potential calculation).
    *   `calc_type%need_exchange` (Logical): True if full/short-range exact exchange is needed.
    *   `calc_type%need_exchange_lr` (Logical): True if long-range exact exchange is needed.
*   **Hybrid Functional Parameters (from `m_inputparam`)**:
    *   `alpha_hybrid` (Real(dp)): Scaling factor for short-range/full exact exchange.
    *   `beta_hybrid` (Real(dp)): Scaling factor for long-range exact exchange.

## Usage Examples

This module is primarily called by the main SCF (Self-Consistent Field) routines within MOLGW.

```fortran
! Conceptual call within an SCF iteration:
USE m_hamiltonian_wrapper
USE m_hamiltonian_tools ! For setup_density_matrix
USE m_scf             ! For en_scf (energy_contributions type) and calc_type
IMPLICIT NONE

TYPE(basis_set) :: orbital_basis
INTEGER :: num_orbitals, num_spin_channels
REAL(DP), ALLOCATABLE :: mo_coeffs(:,:,:), mo_energies(:,:), mo_occupations(:,:)
REAL(DP), ALLOCATABLE :: p_density_matrix(:,:,:), hxc_matrix(:,:,:)
TYPE(energy_contributions) :: current_energies

! ... (orbital_basis, mo_coeffs, mo_occupations, mo_energies are initialized/updated) ...
! ... (current_energies is initialized) ...

CALL setup_density_matrix(mo_coeffs, mo_occupations, p_density_matrix)

CALL calculate_hamiltonian_hxc(orbital_basis, num_orbitals, mo_occupations, mo_coeffs, &
                               p_density_matrix, hxc_matrix, current_energies)

! hxc_matrix now contains V_Hartree + V_XC + alpha*V_Exchange_SR + beta*V_Exchange_LR
! current_energies is updated with E_Hartree, E_XC, E_Exchange, E_Exchange_Hybrid
! This hxc_matrix is then added to the one-body Hamiltonian (T + V_nuc) to form the new Fock/KS matrix.
```

## Dependencies and Interactions

*   **`m_hamiltonian_twobodies`**: This is the core dependency. `m_hamiltonian_wrapper` delegates all actual ERI contractions and XC potential calculations to the specialized subroutines within `m_hamiltonian_twobodies`.
*   **`m_inputparam`**: For global configuration flags (`has_auxil_basis`, `incore_`, `eri3_genuine_`, `alpha_hybrid`, `beta_hybrid`).
*   **`m_scf`**: For the `calc_type` derived type, which contains flags indicating the type of calculation (DFT, HF, hybrid) and thus which potential terms are needed.
*   **`m_hamiltonian_tools`**: For `get_c_matrix_from_p_matrix` when MO coefficients need to be reconstructed from a density matrix for some exchange calculation pathways.
*   **Core Modules**: `m_definitions`, `m_timing`, `m_mpi`, `m_scalapack`, `m_warning`, `m_memory`, `m_basis_set`.
*   This module acts as a high-level dispatcher, simplifying the logic within the main SCF procedure by abstracting the details of how different two-body potential components are computed based on the chosen methodology and settings.
```
