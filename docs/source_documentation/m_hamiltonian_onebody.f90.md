# `m_hamiltonian_onebody.f90`

## Overview

The `m_hamiltonian_onebody` Fortran module in MOLGW is responsible for constructing the one-body components of the Hamiltonian matrix in the Atomic Orbital (AO) basis. This includes the overlap matrix (S), the kinetic energy matrix (T), and the nuclear attraction potential matrix (V<sub>nuc</sub>). It also handles contributions from Effective Core Potentials (ECPs), external static electric fields, and parabolic confinement potentials. Furthermore, the module provides routines to calculate the gradients of these matrices with respect to nuclear coordinates, which are essential for force calculations. It interfaces with external libraries (Libcint or a native Libint via wrappers) for the actual computation of the underlying one-electron integrals.

## Key Components

*   **Overlap Matrix Routines**:
    *   `setup_overlap(basis, s_matrix)`: Computes the overlap matrix S<sub>&mu;&nu;</sub> = <&phi;<sub>&mu;</sub>|&phi;<sub>&nu;</sub>>.
    *   `setup_overlap_mixedbasis(basis1, basis2, s_matrix)`: Computes overlap integrals between basis functions from two different basis sets.
    *   `recalc_overlap(basis_t, basis_p, s_matrix)`: Recalculates necessary parts of the overlap matrix when a "projectile" part of the system (`basis_p`) moves relative to a "target" part (`basis_t`).
    *   `setup_overlap_grad(basis, s_matrix_grad)`: Computes the gradient of the overlap matrix, &nabla;S<sub>&mu;&nu;</sub>.
    *   `setup_overlap_hessian(basis, s_matrix_hess)`: Computes the Hessian (second derivatives) of the overlap matrix, &nabla;&nabla;S<sub>&mu;&nu;</sub>.
    *   `recalc_overlap_grad(basis_t, basis_p, s_matrix_grad)`: Recalculates overlap matrix gradients for moving projectile systems.

*   **Kinetic Energy Matrix Routines**:
    *   `setup_kinetic(basis, hamiltonian_kinetic, timing)`: Computes the kinetic energy matrix T<sub>&mu;&nu;</sub> = <&phi;<sub>&mu;</sub>|-1/2 &nabla;<sup>2</sup>|&phi;<sub>&nu;</sub>>.
    *   `recalc_kinetic(basis_t, basis_p, hamiltonian_kinetic)`: Recalculates kinetic energy matrix elements for systems with moving projectiles.
    *   `setup_kinetic_grad(basis, hamiltonian_kinetic_grad)`: Computes the gradient of the kinetic energy matrix.

*   **Nuclear Attraction and ECP Routines**:
    *   `setup_nucleus(basis, hamiltonian_nucleus, atom_list)`: Computes the nuclear attraction matrix V<sub>&mu;&nu;</sub> = <&phi;<sub>&mu;</sub>|&sum;<sub>A</sub> -Z<sub>A</sub>/|**r**-**R**<sub>A</sub>||&phi;<sub>&nu;</sub>>. `atom_list` can restrict the sum over nuclei.
    *   `setup_para_conf(basis, hamiltonian_nucleus)`: Adds contributions from a parabolic confinement potential (1/2 &omega;<sup>2</sup>r<sup>2</sup>) or a harmonium potential to `hamiltonian_nucleus`.
    *   `recalc_nucleus(basis_t, basis_p, hamiltonian_nucleus)`: Recalculates nuclear attraction matrix elements for systems with moving projectiles.
    *   `setup_nucleus_grad(basis, hamiltonian_nucleus_grad, atom_list, verbose)`: Computes the gradient of the nuclear attraction matrix. This includes Hellmann-Feynman terms and Pulay (basis set derivative) terms.
    *   `setup_nucleus_ecp(basis, hamiltonian_nucleus)`: Adds Effective Core Potential (ECP) contributions to the `hamiltonian_nucleus` matrix. It dispatches to either analytical (for GTH ECPs) or numerical quadrature routines.
    *   `setup_nucleus_ecp_quadrature(basis, hamiltonian_nucleus)`: Computes ECP matrix elements using numerical quadrature for ECPs like NWChem, PSP6, PSP8.
    *   `setup_nucleus_ecp_analytic(basis, hamiltonian_nucleus)`: Computes ECP matrix elements analytically for GTH pseudopotentials, typically using Libcint.

*   **Other One-Electron Property Matrix Routines**:
    *   `setup_rxp_ao(basis, rxp_ao)`: Computes AO matrix elements of the angular momentum operator **L** = **r** &times; **p**.
    *   `setup_giao_rxp_ao(basis, giao_rxp_ao)`: Computes AO matrix elements for Gauge-Independent Atomic Orbitals (GIAO) related to the **r** &times; **p** operator, used in magnetic property calculations.
    *   `setup_electric_field(basis, hext, eext)`: Adds the interaction with an external static electric field (-**E**&middot;**r**) to the one-body Hamiltonian matrix `hext` and calculates the nuclear contribution `eext`.
    *   `setup_dipole_ao(basis, dipole_ao)`: Computes AO matrix elements of the electric dipole operator **r**.
    *   `setup_nabla_ao(basis, nabla_ao)`: Computes AO matrix elements of the gradient operator &nabla;.

## Important Variables/Constants

*   **Input `basis` (Type `basis_set`)**: Contains all AO basis set information (shells, primitives, centers, etc.).
*   **Output Matrices**: `s_matrix` (Overlap), `hamiltonian_kinetic` (Kinetic), `hamiltonian_nucleus` (Nuclear + ECP + External Potentials), and their gradient counterparts.
*   **`MOLGW_has_gradient` (Logical, from `m_definitions`)**: A global flag indicating if the linked integral library supports gradient calculations. If `.FALSE.`, gradient routines will typically issue a warning or error.
*   **Internal C Interface Variables**: The module uses several variables (e.g., `amA`, `contrdepthA`, `A`, `alphaA`, `cA`) to pass data to C-wrapped integral routines (from `m_libint_tools` or `m_libcint_tools`).

## Usage Examples

The subroutines in this module are primarily called during the setup phase of an electronic structure calculation in MOLGW.

```fortran
! Conceptual sequence in MOLGW's SCF setup:
USE m_hamiltonian_onebody
USE m_basis_set
USE m_atoms
! ... (basis_set 'bs', atomic structure 'atoms' are initialized) ...

REAL(DP), ALLOCATABLE :: S_mat(:,:), T_mat(:,:), V_nuc_mat(:,:)

ALLOCATE(S_mat(bs%nbf, bs%nbf))
ALLOCATE(T_mat(bs%nbf, bs%nbf))
ALLOCATE(V_nuc_mat(bs%nbf, bs%nbf))

CALL setup_overlap(bs, S_mat)
CALL setup_kinetic(bs, T_mat)
CALL setup_nucleus(bs, V_nuc_mat) ! Calculates V_nuc (point charges)

IF (nelement_ecp > 0) THEN ! nelement_ecp from m_ecp
  CALL setup_nucleus_ecp(bs, V_nuc_mat) ! Adds ECP contribution to V_nuc_mat
END IF

! H_core = T_mat + V_nuc_mat
! ... (proceed with SCF using H_core and S_mat) ...
```
For force calculations, after obtaining a converged density matrix:
```fortran
! Conceptual: Within a force calculation routine
! ... (SCF converged, density matrices P and R are available) ...
REAL(DP), ALLOCATABLE :: S_grad(:,:,:), T_grad(:,:,:), Vnuc_grad(:,:,:,:)
ALLOCATE(S_grad(bs%nbf, bs%nbf, 3))
ALLOCATE(T_grad(bs%nbf, bs%nbf, 3))
ALLOCATE(Vnuc_grad(bs%nbf, bs%nbf, ncenter_nuclei + 1, 3)) ! Last index for basis deriv terms

CALL setup_overlap_grad(bs, S_grad)
CALL setup_kinetic_grad(bs, T_grad)
CALL setup_nucleus_grad(bs, Vnuc_grad)
! ... (contract these with P and R matrices to get force contributions) ...
```

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_timing`, `m_mpi`, `m_scalapack`, `m_warning`, `m_memory`, `m_inputparam`.
*   **Atomic and Basis Set Information**: `m_basis_set` (provides the `basis_set` type and details of basis functions), `m_atoms` (provides nuclear coordinates and charges).
*   **ECP Data**: `m_ecp` (provides ECP parameters if ECPs are used).
*   **Transformations**: `m_cart_to_pure` (for `transform_libint_to_molgw` and `transform_molgw_to_molgw` which convert integral blocks between Libint/Libcint Cartesian order and MOLGW's possibly pure spherical harmonic basis).
*   **Integral Libraries (Critical Dependency)**:
    *   `m_libint_tools`: Contains Fortran interfaces to a native Libint (version 1 or a MOLGW-internal version) or similar library for computing one-electron integrals (e.g., `libint_overlap`, `libint_kinetic`, `libint_elecpot`, and their gradient counterparts like `libint_overlap_grad`).
    *   `m_libcint_tools`: (If `HAVE_LIBCINT` is defined) Contains Fortran interfaces to the Libcint library (e.g., `cint1e_ovlp_cart`, `cint1e_kin_cart`, `cint1e_rinv_cart`, `cint1e_ipovlp_cart`).
    The choice between these is often handled by preprocessor directives based on how MOLGW was compiled.
*   **DFT Grid**: `m_dft_grid` (indirectly, as `setup_nucleus_ecp_quadrature` uses the grid defined there for numerical integration of ECPs).
*   **Output**: The calculated matrices (S, T, V<sub>nuc</sub>) form the one-electron part of the Fock or Kohn-Sham matrix, which is then used in SCF procedures or other electronic structure methods. The gradient matrices are used in force calculations.
```
