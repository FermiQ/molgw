# `m_density_tools.f90`

## Overview

The `m_density_tools` Fortran module in MOLGW provides a collection of subroutines for calculating various quantities related to the electron density &rho;(r). These include the density itself, its gradient &nabla;&rho;(r), the Laplacian &nabla;<sup>2</sup>&rho;(r), the kinetic energy density &tau;(r), the on-top pair density &Pi;(r), and the electron current density **j**(r). Additionally, it contains specific implementations for calculating Local Density Approximation (LDA) exchange-correlation energies and potentials (e.g., Teter's LDA, custom parameterizations) and the Heyd-Scuseria-Ernzerhof (HSE) screened exchange functional enhancement factor. These tools are essential for Density Functional Theory (DFT) calculations and for analyzing electron distributions.

## Key Components

*   **`setup_atomic_density(basis, rr, rhor)`**:
    *   **Purpose**: Calculates an initial guess for the electron density at a point `rr` by summing spherically averaged atomic densities derived from pre-defined Gaussian fits for each atom.
    *   `basis`: The basis set (used to identify atom types and positions).
    *   `rr(3)`: The point in space.
    *   `rhor`: The calculated atomic density sum at `rr`.

*   **`calc_density_r_batch(occupation, c_matrix, basis_function_r, rhor)`**:
    *   **Purpose**: Computes the electron density &rho;<sub>&sigma;</sub>(r) for each spin channel &sigma; on a batch of real-space grid points. It supports both real and complex molecular orbital coefficients (`c_matrix`).
    *   `occupation`: Occupation numbers of molecular orbitals.
    *   `c_matrix`: Molecular orbital coefficients.
    *   `basis_function_r`: Values of basis functions evaluated at the grid points.
    *   `rhor`: Output array storing &rho;<sub>&sigma;</sub>(r) at each grid point.

*   **`calc_PI_dens_grad_r_batch(...)`**:
    *   **Purpose**: Calculates the on-top pair density &Pi;(r), electron density &rho;<sub>&sigma;</sub>(r), and its gradient &nabla;&rho;<sub>&sigma;</sub>(r) on a batch of grid points. The on-top pair density is needed for some meta-GGA functionals.
    *   `dm2_JK`: Two-body density matrix components related to J and K.
    *   Other parameters are similar to `calc_density_r_batch` with additional inputs for basis function gradients (`bf_gradx`, `bf_grady`, `bf_gradz`).
    *   `PIr`: Output array for &Pi;(r).
    *   `grad_rhor`: Output array for &nabla;&rho;<sub>&sigma;</sub>(r).

*   **`calc_density_gradr_batch(...)`**:
    *   **Purpose**: Computes the electron density &rho;<sub>&sigma;</sub>(r) and its gradient &nabla;&rho;<sub>&sigma;</sub>(r) on a batch of grid points. Supports real and complex `c_matrix`.

*   **`calc_density_current_rr_cmplx(...)`**:
    *   **Purpose**: Calculates the electron current density **j**<sub>&sigma;</sub>(r) on a batch of grid points for complex molecular orbital coefficients.
    *   `jcurdens`: Output array for the three Cartesian components of **j**<sub>&sigma;</sub>(r).

*   **`calc_density_in_disc_cmplx_dft_grid(...)`**:
    *   **Purpose**: Integrates the electron density within a series of cylindrical discs oriented along the z-axis. This is useful for analyzing charge distribution and transfer, particularly in systems with 1D transport characteristics. Outputs results to files.
    *   `num`: An integer typically used for numbering output files.
    *   `time_cur`: Current simulation time, written to output files.

*   **`calc_density_gradr_laplr(...)`**:
    *   **Purpose**: Calculates the gradient of the density &nabla;&rho;(r), the kinetic energy density &tau;(r), and the Laplacian of the density &nabla;<sup>2</sup>&rho;(r) at a single point in space.
    *   Requires the density matrix `p_matrix` and values of basis functions and their first and second derivatives at the point.

*   **`teter_lda_vxc_exc(nr, rhor, vxc, exc)`**:
    *   **Purpose**: Calculates the Local Density Approximation (LDA) exchange-correlation energy density `exc` and potential `vxc` using Teter's parameterization for a given electron density `rhor`.

*   **`my_lda_exc_vxc(nspin, ixc, rhor, exc, vxc)`**:
    *   **Purpose**: Calculates LDA exchange-correlation energy `exc` and potential `vxc` based on a custom Pade-like parameterization.
    *   `ixc`: An integer ID specifying the particular LDA functional variant (e.g., Slater exchange, Teter's full LDA, RPA-based LDA).
    *   Handles both spin-unpolarized (`nspin=1`) and spin-polarized (`nspin=2`) cases.

*   **`my_lda_exc_vxc_mu(mu, rhor, exc, vxc)`**:
    *   **Purpose**: Calculates LDA exchange-correlation energy and potential for a range-separated functional characterized by the parameter `mu` (omega).

*   **`HSE08Fx(omega, ipol, rho, s, Fxhse, d10Fxhse, d01Fxhse)`**:
    *   **Purpose**: Evaluates the Heyd-Scuseria-Ernzerhof (HSE) screened Coulomb exchange functional enhancement factor `Fxhse` and its first derivatives with respect to density `rho` (`d10Fxhse`) and reduced gradient `s` (`d01Fxhse`).
    *   `omega`: Screening parameter for the HSE functional.
    *   `ipol`: Polarization flag (1 for unpolarized, presumably 2 for polarized).

## Important Variables/Constants

*   **Input data for density calculation**:
    *   `occupation(:, :)`: MO occupation numbers.
    *   `c_matrix(:, :, :)`: MO coefficients (can be real or complex).
    *   `basis_function_r(:, :)`: Basis functions evaluated on a grid.
    *   `bf_gradx/y/z(:, :)`: Gradients of basis functions on a grid.
*   **Output density-derived quantities**:
    *   `rhor(:, :)`: Electron density per spin channel.
    *   `grad_rhor(:, :, :)`: Gradient of electron density.
    *   `PIr(:)`: On-top pair density.
    *   `jcurdens(:, :, :)`: Current density.
    *   `tau(:)`: Kinetic energy density.
    *   `lapl_rhor(:)`: Laplacian of electron density.
    *   `exc`, `vxc`: Exchange-correlation energy and potential.
*   **LDA Parameters**: Hardcoded coefficients within `teter_lda_vxc_exc` and `my_lda_exc_vxc` for different LDA functional forms.
*   **HSE Constants**: Hardcoded constants within `HSE08Fx` for the HSE screened exchange functional.

## Usage Examples

These routines are primarily called internally by MOLGW's DFT and analysis modules.

Calculating electron density on a grid:
```fortran
USE m_density_tools
! ... (Assume occupation, c_matrix, basis_function_r_batch are prepared) ...
REAL(DP), ALLOCATABLE :: rho_on_grid(:,:)
ALLOCATE(rho_on_grid(nspin, num_grid_points_in_batch))

CALL calc_density_r_batch(occupation, c_matrix, basis_function_r_batch, rho_on_grid)
```

Calculating LDA exchange-correlation potential for a given density value:
```fortran
USE m_density_tools
REAL(DP) :: density_at_point(1), exc_lda, vxc_lda(1)
INTEGER :: lda_functional_id = 1100 ! Teter's LDA

density_at_point(1) = 0.1_dp ! Example density value
CALL my_lda_exc_vxc(1, lda_functional_id, density_at_point, exc_lda, vxc_lda)
! vxc_lda(1) now contains the LDA potential at that density
```

## Dependencies and Interactions

*   **`m_definitions`**: For `dp`, `pi`, `im` and other constants.
*   **`m_mpi`**: Some routines (e.g., `calc_density_in_disc_cmplx_dft_grid` via `grid%sum`) may involve MPI communication if grid data is distributed.
*   **`m_atoms`**: For atomic coordinates (`xatom`) and charges (`zvalence`, `zatom`) in `setup_atomic_density`.
*   **`m_gaussian`**: (Indirectly) The basis functions evaluated on the grid originate from Gaussian primitives defined in `m_gaussian`.
*   **`m_inputparam`**: (Indirectly) The choice of XC functional (which determines if `my_lda_exc_vxc` or `HSE08Fx` are called) comes from input parameters.
*   **`m_basis_set`**: For the `basis_set` type definition.
*   **`m_hamiltonian_tools`**: For `get_number_occupied_states`.
*   **`m_dft_grid`**: Crucial for routines operating on batches of grid points. It provides the grid coordinates (`rr_grid`), weights (`w_grid`), and the evaluated basis functions on these points (`get_basis_functions_r_batch`).
*   **DFT Modules**: This module is heavily used by the main DFT modules responsible for computing the exchange-correlation energy and potential contributions to the total energy and Kohn-Sham matrix.
*   **Analysis/Plotting Utilities**: Routines like `calc_density_in_disc_cmplx_dft_grid` produce data files for later analysis or visualization.
```
