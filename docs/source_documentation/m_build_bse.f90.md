# `m_build_bse.f90`

## Overview

The `m_build_bse` Fortran module in MOLGW is responsible for constructing the key matrices, commonly referred to as (A-B) and (A+B), which form the core of the eigenvalue problem in linear response theories such as Time-Dependent Density Functional Theory (TDDFT), Random Phase Approximation (RPA), and the Bethe-Salpeter Equation (BSE). These matrices encapsulate the orbital energy differences and the two-electron interaction terms (Coulomb, exchange, and exchange-correlation kernels, or screened Coulomb interaction for BSE). Once constructed, these matrices are typically passed to diagonalization routines (e.g., in `m_block_diago`) to obtain excitation energies and response properties.

## Key Components

*   **`build_amb_apb_common(...)`**:
    *   **Purpose**: Constructs the common parts of the (A-B) and (A+B) matrices. This includes the diagonal orbital energy differences (&epsilon;<sub>a</sub> - &epsilon;<sub>i</sub> - &epsilon;<sub>b</sub> + &epsilon;<sub>j</sub>) and the two-electron integral contributions (Coulomb and exchange terms).
    *   It handles different scenarios:
        *   Whether an auxiliary basis set is used (for RI approximation of ERIs) or direct 4-center ERIs are computed.
        *   Spin considerations for singlet or triplet excitations.
        *   Scaling of interaction terms by the adiabatic connection parameter `lambda`.
        *   Inclusion of exact exchange scaled by `alpha_local` (from hybrid DFT functionals).
    *   **Output**: Populates `amb_matrix` (A-B) and `apb_matrix` (A+B), and also returns the diagonal of (A-B) in the RPA approximation (`amb_diag_rpa`) and the RPA correlation energy part from diagonal elements.

*   **`build_amb_apb_diag_auxil(...)`**:
    *   **Purpose**: Specifically adds the diagonal orbital energy differences to pre-allocated (A-B) and (A+B) matrices when an auxiliary basis set is used. This assumes other contributions (like ERIs) are added separately.

*   **`remove_a_energy_diag(...)`**:
    *   **Purpose**: Modifies a given `a_matrix` by subtracting the diagonal orbital energy differences (&epsilon;<sub>a</sub> - &epsilon;<sub>i</sub>). This creates A' = A - (&Delta;&epsilon;), which is used in some specific ACFD energy calculations.

*   **`get_rpa_correlation(...)`**:
    *   **Purpose**: Calculates a part of the RPA correlation energy derived from the trace of the (A-B) and (A+B) matrices: E<sub>c,RPA_trace_part</sub> = -1/4 Tr(A+B) - 1/4 Tr(A-B). The full RPA correlation energy also includes +1/2 Tr(&Omega;<sub>RPA</sub>), where &Omega;<sub>RPA</sub> are the RPA excitation energies.

*   **`build_apb_hartree_auxil(...)`**:
    *   **Purpose**: Adds the Hartree kernel contribution (2 * (ia|jb) for singlets, 0 for triplets) to the (A+B) matrix using an auxiliary basis (RI approximation). It utilizes precomputed 3-center ERIs.
    *   This is for non-ScaLAPACK execution paths.

*   **`build_apb_hartree_auxil_scalapack(...)`**:
    *   **Purpose**: ScaLAPACK-enabled version of `build_apb_hartree_auxil`. It computes the Hartree kernel contribution to (A+B) in parallel using `PDSYRK`.

*   **`build_apb_tddft(...)`**:
    *   **Purpose**: Adds the exchange-correlation (xc) kernel contribution (K<sub>ia,jb</sub><sup>xc</sup> = &int;&int; &phi;<sub>i</sub>&phi;<sub>a</sub> f<sub>xc</sub>(r,r') &phi;<sub>j</sub>&phi;<sub>b</sub> dr dr') to the (A+B) matrix for TDDFT calculations.
    *   It calls routines from the `m_tddft_fxc` module to evaluate these xc kernel matrix elements.

*   **`build_amb_apb_bse(...)`**:
    *   **Purpose**: Constructs the screened Coulomb interaction (-W) terms for the Bethe-Salpeter Equation. These terms are added to both (A-B) and (A+B) matrices.
    *   This version typically uses 4-center ERIs and information from a statically screened Coulomb interaction `wpol_static` (containing residues and poles of W(ω=0)).

*   **`build_amb_apb_screened_exchange_auxil(...)`**:
    *   **Purpose**: Adds exchange (-&alpha;K) or screened exchange (-&alpha;W, where W is the screened Coulomb interaction) contributions to the (A-B) and (A+B) matrices when using an auxiliary basis set.
    *   It handles both standard exchange (scaled by `alpha_local`) and the full screened interaction for BSE (if `wpol_static` for W(ω=0) is provided).
    *   Supports long-range exchange corrections if `beta_hybrid` is non-zero and specific conditions are met.

## Important Variables/Constants

*   **`amb_matrix` (Real(dp), Intent(inout))**: The (A-B) block matrix.
*   **`apb_matrix` (Real(dp), Intent(inout))**: The (A+B) block matrix.
*   **`wpol` (Type `spectral_function`, Intent(in))**: Contains the definition of the transition space (pairs of occupied `i` and virtual `a` orbitals, and their spins `iaspin`) over which the (A&plusmn;B) matrices are built.
*   **`wpol_static` (Type `spectral_function`, Intent(in))**: Used in BSE calculations, it provides information about the statically screened Coulomb interaction W(&omega;=0), often in terms of its spectral representation (poles and residues).
*   **`alpha_local` (Real(dp), Intent(in))**: Scaling factor for the exchange kernel (K<sub>ia,jb</sub> = (ij|ab)), typically derived from the amount of exact exchange in a hybrid DFT functional.
*   **`lambda` (Real(dp), Intent(in))**: The adiabatic connection coupling constant, used to scale interaction terms in ACFD methods. For standard BSE/TDDFT, &lambda;=1.
*   **`eri_3center_eigen` (Real(dp), Allocatable, 4D array)**: Stores 3-center ERIs (Q|ia) when an auxiliary basis is used. Dimensions are typically (n_aux_basis, n_occ, n_virt, n_spin).
*   **`has_auxil_basis` (Logical)**: A global or module variable indicating if an auxiliary basis is being used for ERIs.

## Usage Examples

These subroutines are internal components of MOLGW, typically called by `m_linear_response` when a TDDFT or BSE calculation is requested.

```fortran
! Conceptual sequence within a routine in m_linear_response:
! ... (Initialize wpol, energy, c_matrix, etc.) ...
! ALLOCATE(amb_matrix(local_rows, local_cols), apb_matrix(local_rows, local_cols))
! ALLOCATE(amb_diag_rpa(nmat_global))
! amb_matrix = 0.0_dp
! apb_matrix = 0.0_dp

! 1. Build common parts (orbital energy differences, Hartree, and exact exchange)
CALL build_amb_apb_common(is_triplet, lambda_val, nmat_global, nbf_global, n_states_total, &
                          orbital_energies, mo_coeffs, wpol, hf_exchange_fraction, &
                          local_rows, local_cols, amb_matrix, apb_matrix, &
                          amb_diag_rpa, rpa_corr_energy_part)

! 2. If TDDFT, add fxc kernel
IF (is_tddft_calculation) THEN
  CALL build_apb_tddft(is_triplet, nmat_global, n_states_total, orbital_basis, mo_coeffs, &
                       occupation_numbers, wpol, local_rows, local_cols, apb_matrix)
END IF

! 3. If BSE (using auxiliary basis for W), add screened exchange
IF (is_bse_calculation .AND. has_auxil_basis) THEN
  CALL build_amb_apb_screened_exchange_auxil(hf_exchange_fraction, lambda_val, &
                                             desc_matrix_AB, wpol, static_polarizability_W, &
                                             local_rows, local_cols, amb_matrix, apb_matrix)
END IF

! ... (Matrices are now ready for diagonalization via m_block_diago) ...
```

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_timing`, `m_warning`, `m_memory`, `m_mpi`, `m_scalapack`, `m_inputparam`.
*   **Spectral Function**: `m_spectral_function` (defines `wpol` and transition tables).
*   **Basis Set**: `m_basis_set` (defines `basis` type, though often only `nbf` is directly used here, assuming orbitals `c_matrix` are already in the MO basis).
*   **ERIs**: `m_eri_ao_mo` (provides `eri_3center_eigen` and routines like `eri_eigen_ri_paral` or `calculate_eri_4center_eigen` for computing necessary two-electron integrals).
*   **TDDFT Kernel**: `m_dft_grid`, `m_tddft_fxc` (for evaluating the exchange-correlation kernel f<sub>xc</sub>).
*   **ScaLAPACK/BLAS**: These routines make extensive use of matrix operations, either via direct BLAS calls (`DSYMM`, `DGER`) or their ScaLAPACK counterparts (`PDSYRK`, etc.) for distributed parallel execution.
*   The constructed (A-B) and (A+B) matrices are subsequently passed to routines in `m_block_diago` for solving the eigenvalue problem.
```
