# `m_dft_grid.f90`

## Overview

The `m_dft_grid` Fortran module in MOLGW is responsible for generating, managing, and providing access to the numerical integration grid used in Density Functional Theory (DFT) calculations. This grid is essential for evaluating exchange-correlation (XC) functionals, which depend on the electron density and potentially its derivatives at points in real space. The module implements atom-centered grids using a combination of radial (Log3 type) and angular (Lebedev-Laikov) quadratures. It also includes atomic partitioning schemes (Becke or Stratmann-Scuseria-Frisch) to assign weights to grid points around each atom and provides routines to evaluate basis functions and their gradients on these grid points, with options for precalculation to optimize performance.

## Key Components

*   **`init_dft_grid(basis, grid_level_in, needs_gradient, precalculate_wfn, batch_size)`**:
    *   **Purpose**: The main initialization routine for the DFT grid. It generates grid points (`rr_grid`) and their corresponding weights (`w_grid`).
    *   The grid density is determined by `grid_level_in` (e.g., `low`, `medium`, `high`), which maps to predefined numbers of radial and angular points.
    *   Uses Log3 radial quadrature for each atom and Lebedev-Laikov angular quadratures.
    *   Applies an atomic partitioning scheme (`partition_scheme` from input, e.g., 'BECKE' or 'SSF') to distribute the total weight among atomic centers.
    *   `needs_gradient` (Logical): If true, indicates that basis function gradients will be required.
    *   `precalculate_wfn` (Logical): If true, basis function values (and gradients if `needs_gradient` is true) are precalculated and stored for a subset of grid points (`ngrid_stored`) to speed up subsequent evaluations. The amount of precalculation is limited by `grid_memory` input parameter.
    *   `batch_size`: The size of batches for processing grid points.

*   **`setup_rhocore_grid()`**:
    *   **Purpose**: Calculates and stores the core electron density (`rhocore`) and its gradient (`rhocore_grad`) on the DFT grid. This is used if Effective Core Potentials (ECPs) with Non-Linear Core Corrections (NLCC) are employed (e.g., PSP6 or PSP8 ECP formats).

*   **`destroy_dft_grid()`**:
    *   **Purpose**: Deallocates all module-level allocatable arrays associated with the DFT grid and precalculated basis function values.

*   **`smooth_step(mu)` (Function)**:
    *   **Purpose**: A smooth step function `0.5 * mu * (3 - mu^2)`, used in Becke's partitioning scheme to create a smooth transition for weights between atoms.

*   **Basis Function Evaluation on Grid**:
    *   **`prepare_basis_functions_r(basis, batch_size)`**: Precalculates and stores values of basis functions on the `ngrid_stored` portion of the grid.
    *   **`prepare_basis_functions_gradr(basis, batch_size)`**: Precalculates and stores gradients of basis functions on the `ngrid_stored` portion of the grid.
    *   **`get_basis_functions_r_batch(basis, igrid, basis_function_r)`**: Retrieves basis function values for a batch of grid points starting at index `igrid`. Uses precalculated values if available, otherwise computes them on-the-fly.
    *   **`get_basis_functions_gradr_batch(basis, igrid, bf_gradx, bf_grady, bf_gradz)`**: Retrieves basis function gradients for a batch. Uses precalculated values if available.
    *   **`calculate_basis_functions_r(basis, rr, basis_function_r)`**: Calculates values of all basis functions at a single point `rr`.
    *   **`calculate_basis_functions_r_batch(basis, rr, basis_function_r)`**: Calculates values of all basis functions for a batch of points `rr(:, :)`. Handles Cartesian to Pure transformations.
    *   **`calculate_basis_functions_gradr(basis, rr, basis_function_gradr)`**: Calculates gradients of all basis functions at a single point `rr`.
    *   **`calculate_basis_functions_gradr_batch(basis, rr, bf_gradx, bf_grady, bf_gradz)`**: Calculates gradients of all basis functions for a batch of points `rr(:, :)`.
    *   **`calculate_basis_functions_laplr(basis, rr, basis_function_gradr, basis_function_laplr)`**: Calculates gradients and Laplacians of basis functions at a single point `rr`.

## Important Variables/Constants

*   **`BATCH_SIZE` (Integer, Parameter)**: Default size for processing grid points in batches (set to 128).
*   **`ngrid` (Integer, Protected)**: The total number of DFT grid points assigned to the current MPI process.
*   **`rr_grid(:, :)` (Real(dp), Protected, Allocatable)**: 2D array (3, `ngrid`) storing the Cartesian coordinates of the grid points.
*   **`w_grid(:)` (Real(dp), Protected, Allocatable)**: 1D array (`ngrid`) storing the quadrature weight associated with each grid point.
*   **`grid_level_in` (Integer)**: Input to `init_dft_grid`, controls the density of the grid (maps to `low`, `medium`, `high`, etc.).
*   **`partition_scheme` (Character string from `m_inputparam`)**: Determines the method for partitioning atomic contributions to the grid weights (e.g., 'BECKE', 'SSF').
*   **`bfr`, `bfgx`, `bfgy`, `bfgz` (Real(dp), Allocatable)**: Arrays used to store precalculated basis function values and their Cartesian gradients on `ngrid_stored` grid points.
*   **`rhocore(:)`, `rhocore_grad(:, :)` (Real(dp), Allocatable)**: Arrays for storing the core electron density and its gradient on the grid, used for NLCC with ECPs.
*   **Pruning and Scheme Parameters**: `pruning_limit`, `aa` (for SSF), `TOL_WEIGHT` (for discarding points with negligible weights).

## Usage Examples

The `m_dft_grid` module is primarily used internally by other modules in MOLGW, especially those involved in DFT calculations.

Initialization of the grid:
```fortran
USE m_dft_grid
USE m_basis_set
USE m_inputparam ! For grid_quality, partition_scheme etc.

TYPE(basis_set) :: orbital_basis
LOGICAL :: needs_gradients_for_xc, precalculate_wavefunctions
INTEGER :: effective_batch_size

! ... orbital_basis is initialized ...
! ... set grid_quality, needs_gradients_for_xc, precalculate_wavefunctions, effective_batch_size from inputs ...

CALL init_dft_grid(orbital_basis, grid_quality, needs_gradients_for_xc, &
                   precalculate_wavefunctions, effective_batch_size)
```

Looping over grid batches to compute density (conceptual):
```fortran
USE m_dft_grid
USE m_density_tools
! ... (grid initialized, MO coefficients c_matrix and occupations available) ...
INTEGER :: igrid_start, igrid_end, n_batch_points
REAL(DP), ALLOCATABLE :: bfr_batch(:,:), rho_batch(:,:)

DO igrid_start = 1, ngrid, BATCH_SIZE  ! ngrid is from m_dft_grid
  igrid_end = MIN(ngrid, igrid_start + BATCH_SIZE - 1)
  n_batch_points = igrid_end - igrid_start + 1

  ALLOCATE(bfr_batch(orbital_basis%nbf, n_batch_points))
  ALLOCATE(rho_batch(nspin, n_batch_points))

  CALL get_basis_functions_r_batch(orbital_basis, igrid_start, bfr_batch)
  CALL calc_density_r_batch(occupation_numbers, c_matrix, bfr_batch, rho_batch)
  
  ! ... use rho_batch and w_grid(igrid_start:igrid_end) for XC calculations ...

  DEALLOCATE(bfr_batch, rho_batch)
END DO
```

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_warning`, `m_mpi`, `m_memory`, `m_timing`, `m_inputparam`.
*   **Basis Set & Atomic Info**: `m_basis_set` (for `basis_set` type and evaluating basis functions), `m_atoms` (for atomic coordinates `xbasis` and types `zbasis`), `m_elements` (for covalent radii).
*   **Transformations**: `m_cart_to_pure` (for transforming basis functions/derivatives if 'PURE' type is used).
*   **ECPs**: `m_ecp` (for non-linear core correction data if `setup_rhocore_grid` is used).
*   **Numerical Tools**: `m_numerical_tools` (potentially for `coeffs_gausslegint`, though Log3 grid is custom).
*   **Lebedev Quadrature**: Directly calls `ldxxxx` subroutines (from `lebedev_quadrature.f`) for angular grid points and weights.
*   **DFT Calculation Engine**: This module is fundamental for DFT calculations, providing the grid points and weights for numerical integration of exchange-correlation functionals. It works in tandem with modules like `m_density_tools` (which calculates density on the grid) and XC functional evaluation routines.
*   **Memory Management**: The precalculation of basis function values is managed based on the `grid_memory` input parameter to balance speed and memory usage.
```
