# `m_atoms.f90`

## Overview

The `m_atoms` Fortran module in MOLGW is responsible for managing all information related to the atomic composition and geometry of the molecular system under study. This includes storing atomic numbers, Cartesian coordinates, valence charges, defining centers for basis functions (which can include ghost atoms), and handling projectiles. It also provides routines for calculating nucleus-nucleus repulsion energy and forces, determining geometric properties (like symmetry), and updating atomic positions during geometry optimizations.

## Key Components

*   **Module Variables**:
    *   `ncenter_nuclei` (Integer): Total number of atomic nuclei.
    *   `ncenter_basis` (Integer): Total number of centers where basis functions are placed (can differ from `ncenter_nuclei` if ghost atoms are present).
    *   `zatom(:)` (Real(dp), Allocatable): Array storing the atomic numbers (nuclear charges) of the nuclei.
    *   `xatom(:, :)` (Real(dp), Allocatable): 2D array storing the Cartesian coordinates (3, `ncenter_nuclei`) of the nuclei.
    *   `zvalence(:)` (Real(dp), Allocatable): Array storing the effective valence charge for each nucleus (used for ECPs or specific analyses).
    *   `xbasis(:, :)` (Real(dp), Allocatable): 2D array storing the Cartesian coordinates (3, `ncenter_basis`) of the centers for basis functions.
    *   `zbasis(:)` (Integer, Allocatable): Array storing the atomic numbers associated with each basis center (identifies the element type for basis set selection).
    *   `vel_nuclei(:, :)`, `vel_basis(:, :)` (Real(dp), Allocatable): Velocities associated with nuclei and basis centers, relevant for dynamics.
    *   `force(:, :)` (Real(dp), Allocatable): 2D array storing the Cartesian components of the force (3, `ncenter_nuclei`) acting on each nucleus.
    *   `nprojectile` (Integer): Number of projectile particles (typically 0 or 1).
    *   `nghost_` (Integer): Number of ghost atom centers.
    *   `nbond` (Integer): Number of identified covalent bonds.
    *   `inversion`, `linear`, `planar` (Logical): Flags indicating geometric properties of the molecule.
    *   `xcenter(3)`, `xnormal(3)` (Real(dp)): Coordinates of the center of mass/inversion and normal vector for planar molecules.

*   **`init_atoms(...)`**:
    *   Initializes the atomic system from user input (read externally). Sets up `xatom`, `zatom`, `xbasis`, `zbasis`, and related arrays.
    *   Differentiates between regular atoms, ghost atoms (basis functions without a nucleus), and projectiles.
    *   Adjusts projectile charge based on `excit_name` and `projectile_charge_scaling`.
    *   Checks for interatomic distances that are too small.
    *   Identifies covalent bonds based on covalent radii.
    *   Calls `find_inversion()` and determines linearity/planarity.

*   **`atoms_core_states()` (Function)**:
    *   Calculates and returns the total number of core electrons in the system based on `zatom` and `zvalence` using data from `m_elements`.

*   **`get_bondcenter(ibond, xbond)`**:
    *   Calculates the midpoint coordinates `xbond(3)` for a given bond index `ibond`.

*   **`change_position_one_atom(iatom, xposition)`**:
    *   Updates the coordinates in `xatom` for a specific atom `iatom`.

*   **`change_basis_center_one_atom(iatom, xposition)`**:
    *   Updates the coordinates in `xbasis` for a specific basis center `iatom`.

*   **`destroy_atoms()`**:
    *   Deallocates all allocatable arrays defined in the module to free memory.

*   **`relax_atoms(lbfgs_plan, etotal)`**:
    *   Updates atomic positions (`xatom` and `xbasis`) using forces stored in the `force` array, employing the L-BFGS algorithm provided by `m_lbfgs`.
    *   Limits the maximum displacement of any atom in a single optimization step.

*   **`output_positions()`**:
    *   Prints a formatted list of all atomic centers (nuclei, ghosts, projectile) with their labels, element names, and coordinates in both Bohr and Angstrom units.

*   **`output_projectile_position()`**:
    *   Specifically prints the position of the projectile particle if one is present.

*   **`nucleus_nucleus_energy(energy)`**:
    *   Calculates the classical electrostatic repulsion energy between all pairs of nuclei using their `zvalence` charges.

*   **`nucleus_nucleus_force()`**:
    *   Calculates the contribution to atomic forces arising from nucleus-nucleus electrostatic repulsion. Stores results in `force_nuc_nuc`.

*   **`find_inversion()`**:
    *   Determines if the molecule (excluding any projectile) possesses an inversion center. Sets the `inversion` logical flag and `xcenter`.

*   **`same_element(icenter, jcenter)` (Function)**:
    *   Returns `.TRUE.` if the nuclei at `icenter` and `jcenter` are of the same element, `.FALSE.` otherwise.

## Important Variables/Constants

*   **`tol_geom` (Real(dp), Parameter, Private)**: A small tolerance (1.0e-5 Bohr) used for geometric comparisons, such.as checking for inversion symmetry or co-linearity/co-planarity.
*   **`bohr_A` (Real(dp))**: Conversion factor from Bohr to Angstrom (implicitly used from `m_definitions` or `m_elements`).
*   The module stores various force components in separate arrays (e.g., `force_nuc_nuc`, `force_kin`, `force_exc`) although their direct usage or summation into the main `force` array is handled by other modules responsible for force calculations.

## Usage Examples

The `m_atoms` module is primarily used internally. `init_atoms` is called during program initialization based on user input.
```fortran
! Conceptual usage during MOLGW initialization
CALL init_atoms(natom_from_input, nghost_from_input, nucleus_without_basis_flags, &
                atomic_numbers_from_input, coordinates_from_input, projectile_velocity, &
                calculate_forces_flag, excitation_type_name, projectile_scaling_factor)

! Later, to get nucleus-nucleus repulsion energy:
REAL(DP) :: enuc
CALL nucleus_nucleus_energy(enuc)

! To output atomic positions to log:
CALL output_positions()
```

## Dependencies and Interactions

*   **`m_definitions`**: Provides the `dp` kind parameter for double-precision real numbers.
*   **`m_warning`**: Used for error handling (`die`) and issuing warnings (`issue_warning`).
*   **`m_elements`**: Provides element-specific data such as covalent radii (`element_covalent_radius`), names (`element_name_long`, `element_name`), and core electron counts (`element_core`).
*   **`m_linear_algebra`**: Provides vector operations like `cross_product`, `NORM2`, `DOT_PRODUCT`.
*   **`m_lbfgs`**: The `relax_atoms` subroutine depends on this module for L-BFGS geometry optimization steps.
*   **Input System**: `init_atoms` receives data that is typically read from MOLGW's main input file (e.g., atom coordinates, types, projectile information).
*   **Basis Set Modules**: This module provides the definition of atomic centers (`xbasis`, `zbasis`) which are then used by modules like `m_basis_set` to construct the actual basis functions.
*   **Hamiltonian & Force Calculation Modules**: The atomic coordinates and charges stored here are fundamental inputs for calculating various terms in the Hamiltonian (e.g., electron-nucleus attraction) and for computing forces on atoms.
```
