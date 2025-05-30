# `m_noft.f90`

## Overview

The `m_noft` Fortran module in MOLGW serves as an interface and driver for Natural Orbital Functional Theory (NOFT) calculations. NOFT is a quantum chemical method that aims to determine the ground state energy and properties of a system by minimizing a Natural Orbital Functional (NOF) with respect to the natural orbitals (NOs) and their natural occupation numbers (NONs). This module prepares the necessary inputs (one- and two-electron integrals in the NO basis, control parameters) and calls an external NOFT solver routine (`run_noft`) to perform the energy minimization. It supports various NOFs, range-separated NOFT+DFT, and both real and complex arithmetic for relativistic calculations.

## Key Components

*   **Module-Level Private Variables**:
    *   `nstate_noft`: Number of natural orbitals in the active space.
    *   `nstate_frozen`: Number of frozen core orbitals.
    *   `irs_noft`: Integer flag indicating the type of range-separation for NOFT+DFT (0: none, 1: srNOFT+lrDFT, 2: srDFT+lrNOFT).
    *   `ExcDFT`: Stores the DFT exchange-correlation energy in NOFT+DFT.
    *   `AhCORE(:, :)`, `AhCORE_cmplx(:, :)`: Stores the core Hamiltonian (T + V<sub>ext</sub>) in the AO basis (real or complex).
    *   `T_Vext(:)`: Stores the diagonal elements of the core Hamiltonian in the NO basis.
    *   `basis_pointer`: A pointer to the `basis_set` type.

*   **`noft_energy(basis, occupation, Enoft, Vnn, Aoverlap, c_matrix, c_matrix_rel, hkin, hnuc, hkin_nuc_rel)` (Subroutine)**:
    *   **Purpose**: This is the main public subroutine that drives the NOFT calculation.
    *   **Inputs**:
        *   `basis`: The atomic orbital basis set.
        *   `occupation`: Initial guess for MO occupation numbers (often used just to determine electron count).
        *   `Vnn`: Nuclear repulsion energy.
        *   `Aoverlap` (Optional): AO overlap matrix.
        *   `hkin`, `hnuc` (Optional): AO kinetic and nuclear attraction matrices.
        *   `c_matrix` (Optional, InOut): Initial guess for MO/NO coefficients (real), updated with optimized NOs.
        *   `hkin_nuc_rel` (Optional): Complex core Hamiltonian for x2c.
        *   `c_matrix_rel` (Optional, InOut): Initial guess for complex MO/NO coefficients (x2c), updated.
    *   **Outputs**:
        *   `Enoft`: The final NOFT total energy.
        *   `occupation`: Updated with optimized NONs.
        *   `c_matrix` / `c_matrix_rel`: Updated with optimized NO coefficients.
    *   **Workflow**:
        1.  Initializes parameters based on input keywords (e.g., `noft_functional`, `noft_nscf`, `noft_Lpower`, `noft_complex`, `noft_dft`).
        2.  Sets up the core Hamiltonian (`AhCORE` or `AhCORE_cmplx`).
        3.  Sets up initial guess for NOs (`NO_COEF` or `NO_COEF_cmplx`) from input `c_matrix` or H<sub>core</sub> diagonalization.
        4.  Determines active space parameters (number of frozen, active, coupled pairs).
        5.  If range-separated NOFT+DFT is requested (`irs_noft /= 0`), initializes the DFT grid.
        6.  Calls the external `run_noft` subroutine (not defined in this module), passing `mo_ints` (or `mo_ints_x2c`) as a callback function to provide integrals in the current NO basis.
        7.  If range-separated NOFT+DFT, performs a second call to `run_noft` to compute the DFT energy component using the optimized NOs/NONs, and adds `ExcDFT` to `Enoft`.
        8.  Updates output `occupation` and `c_matrix` (or `c_matrix_rel`) with optimized NONs and NOs.
        9.  Optionally prints NOs/density to files if requested by input flags.

*   **`mo_ints(...)` (Subroutine, passed as argument to `run_noft`)**:
    *   **Purpose**: Callback routine for the external NOFT solver. It provides the one-body core Hamiltonian (`hCORE`) and two-body electron repulsion integrals (`ERImol`, and `ERImolJsr`/`ERImolLsr` for range-separated parts) transformed into the current natural orbital basis.
    *   It takes the current `NO_COEF` (or `NO_COEF_cmplx`) from the NOFT solver.
    *   Transforms AO `AhCORE` to the NO basis to get `hCORE`.
    *   Transforms AO ERIs (from `m_eri_ao_mo`) to the NO basis to get `ERImol` (and range-separated components if needed).
    *   If `noft_dft='yes'` and `irs_noft /= 0` (range-separated NOFT+DFT), it calculates the DFT XC contribution using `dft_exc_vxc_batch` (from `m_hamiltonian_twobodies`) with the current NOs and NONs, and adds it to `hCORE`. The XC energy is stored in `ExcDFT`.

*   **`mo_ints_x2c(...)` (Subroutine, passed as argument to `run_noft`)**:
    *   A version of `mo_ints` specifically for x2c relativistic calculations, handling complex arithmetic for orbitals and integrals.

## Important Variables/Constants

*   **Input Parameters (from `m_inputparam`)**:
    *   `noft_functional` (Character string): Specifies the Natural Orbital Functional to use (e.g., 'PNOF7', 'GNOF', 'PCCD').
    *   `noft_nscf` (Integer): Maximum number of NOFT SCF iterations.
    *   `noft_tolE` (Real(dp)): Convergence threshold for NOFT energy.
    *   `noft_Lpower` (Real(dp)): Parameter for the 'POWER' functional or related custom functionals.
    *   `noft_complex` (Character string 'yes'/'no'): Controls whether complex arithmetic is used (for x2c or specific NOFs).
    *   `noft_dft` (Character string 'yes'/'no'): Enables NOFT+DFT calculations.
    *   `noft_rsinter` (Character string 'yes'/'no'): Enables range-separated NOFT+DFT.
    *   `alpha_hybrid`, `beta_hybrid`, `gamma_hybrid`: Parameters for range-separated hybrids if NOFT+DFT is used.
*   **`nstate_noft` (Integer)**: Number of active natural orbitals in the NOFT calculation.
*   **`irs_noft` (Integer)**: Flag indicating the type of range separation in NOFT+DFT.
    *   0: No range separation (standard NOFT or NOF part of NOFT+srDFT).
    *   1: srNOFT + lrDFT (short-range NOF, long-range DFT).
    *   2: srDFT + lrNOFT (short-range DFT, long-range NOF).

## Usage Examples

The `m_noft` module is invoked when the user specifies `scf = 'NOFT'` or a NOF for `postscf` in the MOLGW input file.
```
&molgw
  xyz_file = "h2.xyz"
  basis    = "cc-pVDZ"
  scf      = "NOFT"
  noft_functional = "PNOF7"
  noft_nscf = 200
  noft_tolE = 1.0e-8
&end
```
MOLGW's main routine would call `noft_energy` to perform the calculation.

## Dependencies and Interactions

*   **External NOFT Solver (`run_noft`)**: This is a critical dependency. The `m_noft` module prepares inputs and calls `run_noft` (which is not part of this module and is assumed to be linked externally) to perform the actual NOF minimization. The `mo_ints` or `mo_ints_x2c` subroutines are passed as arguments to `run_noft` to provide the necessary integrals in the current NO basis.
*   **Core Modules**: `m_definitions`, `m_mpi`, `m_cart_to_pure`, `m_inputparam`, `m_warning`, `m_memory`, `m_timing`.
*   **System Description**: `m_basis_set`, `m_atoms`.
*   **Hamiltonian and ERI Modules**:
    *   `m_hamiltonian_tools`: For matrix transformations and potentially initial NO guess from density matrix.
    *   `m_hamiltonian_onebody`: For `hkin` and `hnuc` (core Hamiltonian components).
    *   `m_hamiltonian_twobodies`: For `dft_exc_vxc_batch` if NOFT+DFT is used.
    *   `m_eri_ao_mo`: For transforming ERIs from AO to MO (NO) basis within `mo_ints`.
*   **DFT Components (for NOFT+DFT)**:
    *   `m_dft_grid`: For initializing the DFT grid if `noft_dft='yes'`.
    *   `m_density_tools`: (Indirectly via `dft_exc_vxc_batch`).
*   The module acts as an interface between MOLGW's general infrastructure (basis sets, AO integrals, SCF framework) and a specialized external NOFT optimization engine.
```
