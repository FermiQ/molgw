# `m_acfd.f90`

## Overview

The `m_acfd` Fortran module in MOLGW is dedicated to calculating the total energy using methods based on the Adiabatic Connection Fluctuation-Dissipation Theorem (ACFD). This primarily involves computing the correlation energy via different flavors of the Random Phase Approximation (RPA) and related formalisms (e.g., RPAX, BSE-I). It integrates these correlation energies with other energy components (Hartree-Fock or DFT energy) to yield a final total energy.

## Key Components

*   **`acfd_total_energy(basis, nstate, occupation, energy, c_matrix, en_mbpt)`**:
    *   This is the main subroutine of the module. It orchestrates the calculation of the ACFD correlation energy based on the method specified by the `postscf` input parameter.
    *   It handles various RPA-like methods including standard RPA, RPA+ (RPA plus a short-range correlation correction from Libxc), RPAX (RPA with exchange contributions), and methods involving integration over the coupling constant &lambda; (RPA-I, RPAX-I, BSE-I).
    *   The results are stored in the `en_mbpt` (energy_contributions type) variable, which aggregates different parts of the total energy.

*   **`calculate_ec_acft(desc_x, a_matrix, b_matrix, x_matrix, y_matrix, erpa)`**:
    *   A specialized subroutine used by ACFD methods that require numerical integration over the coupling constant &lambda; (e.g., RPA-I, BSE-I).
    *   It computes a specific formula for the correlation energy at a given &lambda; point: `Ec(&lambda;) = -0.5 * Tr(A(&lambda;)) + 0.5 * Tr(X(&lambda;)^T * (A(&lambda;)*X(&lambda;) + B(&lambda;)*Y(&lambda;))) + 0.5 * Tr(Y(&lambda;)^T * (A(&lambda;)*Y(&lambda;) + B(&lambda;)*X(&lambda;)))`.
    *   The matrices `A` and `B` are components of the linear response kernel (related to orbital energy differences and ERIs), while `X` and `Y` are solution vectors from the Casida-like equations.
    *   This calculation is performed using ScaLAPACK routines for parallel execution if available, otherwise, it falls back to standard BLAS/LAPACK operations.

## Important Variables/Constants

*   **`postscf` (character string)**: An input parameter (from `m_inputparam`) that dictates the specific ACFD/RPA methodology to be employed (e.g., 'RPA', 'RPA+', 'RPAX', 'RPA-I', 'BSE-I').
*   **`acfd_nlambda` (integer)**: An input parameter (from `m_inputparam`) specifying the number of quadrature points for the numerical integration over the coupling constant &lambda; in methods like RPA-I.
*   **`kappa_hybrid` (real(dp))**: An input parameter (from `m_inputparam`) used to scale the calculated correlation energy, relevant in the context of double-hybrid DFT methods.
*   **`en_mbpt` (type `energy_contributions`)**: A derived-type variable that accumulates various energy terms (nuclear repulsion, kinetic, electron-nucleus, Hartree, exchange, and the correlation energy computed by this module).
*   **`wpol` (type `spectral_function`)**: A derived-type variable holding information about the frequency-dependent polarizability, which is central to RPA calculations.
*   **`a_matrix`, `b_matrix` (real(dp), allocatable, 2D arrays)**: Represent the A and B matrices in the Casida equation formulation of linear response theory. `A` typically involves orbital energy differences and Coulomb/exchange integrals; `B` involves only Coulomb/exchange integrals.
*   **`x_matrix`, `y_matrix` (real(dp), allocatable, 2D arrays)**: Represent the (X+Y) and (X-Y) solution vectors (or a transformation thereof) of the Casida equations, used in constructing the frequency-dependent polarizability and in the `calculate_ec_acft` subroutine.
*   **`dft_xc_tmp` (type `dft_xc_info`)**: An array of derived-type variables used to configure calls to the `Libxc` library for calculating the short-range correlation correction in RPA+ methods.

## Usage Examples

The `m_acfd` module is used internally by MOLGW. When a user specifies a `postscf` method like `'RPA'` in the MOLGW input file, the `acfd_total_energy` subroutine is called after the initial SCF (HF or DFT) calculation.

```fortran
! Conceptual call within MOLGW's post-SCF driver routine
! ... (SCF calculation completed, en_scf contains SCF energy components) ...

IF (TRIM(postscf) == 'RPA' .OR. TRIM(postscf) == 'RPA+' .OR. &
    TRIM(postscf) == 'RPAX' .OR. TRIM(postscf) == 'RPA_IM' .OR. &
    TRIM(postscf) == 'RPAP_IM' .OR. TRIM(postscf) == 'RPA+_IM' .OR. &
    TRIM(postscf) == 'RPA-I' .OR. TRIM(postscf) == 'RPAX-I' .OR. &
    TRIM(postscf) == 'BSE-I' ) THEN
  
  ! Initialize en_mbpt with relevant components from en_scf
  en_mbpt%nuc_nuc = en_scf%nuc_nuc
  en_mbpt%kinetic = en_scf%kinetic
  en_mbpt%nucleus = en_scf%nucleus
  en_mbpt%hartree = en_scf%hartree
  en_mbpt%exx     = en_scf%exx      ! Exact exchange from HF or hybrid DFT
  en_mbpt%total   = en_scf%total    ! Initial total energy (could be HF or hybrid DFT)

  CALL acfd_total_energy(basis, nstate_occ, occupation_numbers, orbital_energies, &
                         mo_coefficients, en_mbpt)
  
  PRINT *, "Final ACFD Total Energy (Ha):", en_mbpt%total
ENDIF
```

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_timing`, `m_warning`, `m_memory`, `m_mpi`, `m_scalapack`, `m_inputparam`.
*   **Basis Set & Grid**: `m_basis_set`, `m_dft_grid` (for RPA+ correction).
*   **Hamiltonian Components**: `m_hamiltonian_twobodies` (for integrals feeding into A and B matrices).
*   **Response Theory**: `m_spectral_function`, `m_spectra`, `m_linear_response` (these are heavily used to calculate polarizabilities and solve Casida equations).
*   **SCF Data**: `m_scf` (provides initial orbitals, energies, and energy components).
*   **Numerical Tools**: `m_numerical_tools` (for Gauss-Legendre quadrature coefficients used in &lambda;-integration).
*   **External Libraries**:
    *   **`Libxc`**: (Optional, if compiled) Used to compute short-range correlation corrections for 'RPA+' methods.
    *   **BLAS/LAPACK**: For dense matrix algebra (e.g., `DSYMM`, `DGEMM`, matrix trace).
    *   **ScaLAPACK**: (Optional, if compiled and MPI is enabled) For parallel distributed matrix algebra (e.g., `PDSYMM`, `PDGEMM`, `PDLATRA`).
*   **Interaction Flow**:
    1.  `acfd_total_energy` is called with SCF results.
    2.  For 'RPA+' methods, it computes a short-range correction using `Libxc` via `dft_exc_vxc_batch`.
    3.  Based on `postscf` type:
        *   Direct RPA/RPAX: Calls `polarizability` to get the correlation energy directly.
        *   Imaginary frequency RPA: Calls `polarizability_grid_scalapack`.
        *   &lambda;-integrated methods (RPA-I, RPAX-I, BSE-I):
            *   Obtains A and B matrices from `polarizability` (with &lambda;=1).
            *   Loops over &lambda; points:
                *   Calls `polarizability` to get X(&lambda;) and Y(&lambda;) matrices.
                *   Calls `calculate_ec_acft` to compute energy contribution at current &lambda;.
            *   Sums contributions using Gauss-Legendre weights.
    4.  Updates `en_mbpt` with the calculated correlation energy and prints results.
```
