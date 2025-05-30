# `m_basis_set.f90`

## Overview

The `m_basis_set` Fortran module is a cornerstone of MOLGW, responsible for the definition, initialization, and management of Gaussian basis sets. It handles both orbital basis sets and auxiliary basis sets. The module can read basis set information from external files or generate even-tempered basis sets. It supports Cartesian and pure spherical harmonic Gaussian functions and provides functionalities for evaluating basis functions, their derivatives, and various one-electron integrals. It also includes routines for preparing basis set information for use with the Libcint library.

## Key Components

*   **`basis_function` (Type)**:
    *   A derived data type representing a single contracted Gaussian basis function. It stores:
        *   Angular momentum (Cartesian `nx,ny,nz` or pure `am,mm`).
        *   Atomic center coordinates (`x0`) and velocity (`v0`).
        *   Contraction depth (`ngaussian`), primitive Gaussian exponents (`g(:)%alpha`), and contraction coefficients (`coeff(:)`).
        *   An array of `gaussian` type (from `m_gaussian`) for primitive Gaussians.
        *   Shell information (`shell_index`, `index_in_shell`).

*   **`shell_type` (Type)**:
    *   Represents a shell of basis functions that share the same center, angular momentum (`am`), and set of primitive Gaussians (exponents `alpha` and original coefficients `coeff`).
    *   Stores start/end indices for its functions within the larger basis set arrays.

*   **`basis_set` (Type)**:
    *   The main derived data type representing an entire basis set. It contains:
        *   `nbf`: Total number of final basis functions (pure or Cartesian).
        *   `nbf_cart`: Total number of underlying Cartesian basis functions.
        *   `nshell`: Total number of shells.
        *   `ammax`: Maximum angular momentum in the basis.
        *   `gaussian_type`: 'CART' or 'PURE'.
        *   `bff(:)`: Array of `basis_function` type for the final (possibly pure) basis functions.
        *   `bfc(:)`: Array of `basis_function` type for the Cartesian basis functions.
        *   `shell(:)`: Array of `shell_type`.
        *   `LIBCINT_*` members: Data structured for interfacing with the Libcint library.

*   **`init_basis_set(...)`**:
    *   Initializes a `basis_set` object. It reads basis set definitions from files (path constructed from `basis_path`, element name, and `basis_name` or `ecp_basis_name`).
    *   Supports generation of "EVEN_TEMPERED" basis sets based on `even_tempered_alpha`, `even_tempered_beta`, and `even_tempered_n_list` parameters.
    *   Handles the distinction and potential transformation between Cartesian and pure spherical harmonic Gaussians based on `gaussian_type`.

*   **`init_auxil_basis_set_auto(...)`**:
    *   Generates an auxiliary basis set automatically using the "Auto" or "PAuto" method described by Yang, Rendell, and Frisch. This is typically for density fitting purposes.
    *   Parameters `auto_auxil_fsam` and `auto_auxil_lmaxinc` control the generation.

*   **`split_basis_set(basis, basis_t, basis_p)`**:
    *   Divides an existing `basis_set` into two new sets: `basis_t` (target/fixed centers) and `basis_p` (projectile/moving centers) based on the velocity `v0` of the basis function centers.

*   **`moving_basis_set(new_basis)`**:
    *   Updates the coordinates (`x0`) of basis functions within `new_basis` that are marked as moving (i.e., have non-zero `v0`). This is used for systems with moving projectiles.

*   **`init_basis_function(normalized, ng, nx, ny, nz, ...)`**:
    *   Initializes a single Cartesian `basis_function` object from primitive Gaussians, applying normalization if `normalized` is true.

*   **`init_basis_function_pure(normalized, ng, am, mm, ...)`**:
    *   Initializes a single pure spherical harmonic `basis_function` object.

*   **Integral & Evaluation Routines**:
    *   `eval_basis_function(bf, x)`: Evaluates the amplitude of basis function `bf` at point `x`.
    *   `eval_basis_function_grad(bf, x)`: Evaluates the gradient of `bf` at `x`.
    *   `eval_basis_function_lapl(bf, x)`: Evaluates the Laplacian of `bf` at `x`.
    *   `overlap_basis_function(bf1, bf2, overlap)`: Computes <bf1|bf2>.
    *   `kinetic_basis_function(bf1, bf2, kinetic)`: Computes <bf1|-1/2 &nabla;<sup>2</sup>|bf2>.
    *   `nucleus_basis_function(bf1, bf2, zatom, x, nucleus_pot)`: Computes <bf1|Z/|r-R_Nuc||bf2>.
    *   `basis_function_prod(bf1, bf2, bfprod)`: Computes the product bf1*bf2 as a sum of new (unnormalized) basis functions.
    *   `basis_function_dipole(bf1, bf2, dipole)`: Computes <bf1|r|bf2>.
    *   `basis_function_quadrupole(bf1, bf2, quad)`: Computes <bf1|r*r'|bf2>.
    *   `basis_function_gos(bf1, bf2, qvec, gos_bf1bf2)`: Computes generalized oscillator strength matrix elements.
    *   `basis_function_fourier(bf, qvec)`: Computes the Fourier transform of basis function `bf`.

*   **Serialization Routines**:
    *   `write_basis_set`, `read_basis_set`, `write_basis_function`, `read_basis_function`, `write_basis_shell`, `read_basis_shell`: Subroutines for writing and reading basis set data to/from a Fortran unit, likely for checkpointing.

## Important Variables/Constants

*   **`basis%gaussian_type` (Character(len=4))**: Determines if the final basis functions (`bff`) are 'CART' (Cartesian) or 'PURE' (spherical harmonics).
*   **`MOLGW_LMAX` (Preprocessor Constant)**: Defines the maximum angular momentum supported by the linked Libcint library. `init_basis_set` checks against this.
*   **Normalization**: The `init_basis_function` routine normalizes the contracted basis functions unless `normalized` is explicitly false. Product Gaussians from `basis_function_prod` are explicitly unnormalized.

## Usage Examples

This module is mainly used by other parts of MOLGW.

Loading a standard basis set:
```fortran
USE m_basis_set
USE m_atoms ! Assuming ncenter_basis, zbasis, xbasis, vel_basis are initialized from m_atoms
TYPE(basis_set) :: orbital_basis
CHARACTER(LEN=100), ALLOCATABLE :: basis_name_per_atom(:)
CHARACTER(LEN=100), ALLOCATABLE :: ecp_name_per_atom(:)
CHARACTER(LEN=200) :: path_to_basis_files
CHARACTER(LEN=4)   :: type_of_gaussians  ! 'CART' or 'PURE'
! ... initialize basis_name_per_atom, ecp_name_per_atom, path_to_basis_files, type_of_gaussians ...
! ... (typically from user input or m_inputparam) ...

CALL init_basis_set(path_to_basis_files, basis_name_per_atom, ecp_name_per_atom, &
                    type_of_gaussians, &
                    even_tempered_alpha_val, even_tempered_beta_val, & ! if using even-tempered
                    even_tempered_n_list_str, orbital_basis) 
```

Generating an auto-auxiliary basis:
```fortran
USE m_basis_set
TYPE(basis_set) :: main_orbital_basis, aux_basis
CHARACTER(LEN=100) :: aux_basis_name_specifier ! e.g., "Auto" or "PAuto"
REAL(DP) :: fsam_param
INTEGER :: lmaxinc_param
! ... initialize main_orbital_basis, aux_basis_name_specifier, fsam_param, lmaxinc_param ...

CALL init_auxil_basis_set_auto(aux_basis_name_specifier, main_orbital_basis, &
                               'CART', fsam_param, lmaxinc_param, aux_basis)
```

## Dependencies and Interactions

*   **`m_definitions`**: For `dp` (double precision) and `C_INT`, `C_DOUBLE` type kinds.
*   **`m_warning`**: For error handling (`die`) and warnings.
*   **`m_timing`**: For performance profiling of basis set initialization.
*   **`m_mpi`**: (Potentially for distributed loading or data, though not explicit in the provided routines).
*   **`m_elements`**: For `element_name`.
*   **`m_string_tools`**: For `orbital_momentum_name` and `append_to_list`.
*   **`m_atoms`**: Crucial dependency, providing the number of basis centers (`ncenter_basis`), their atomic numbers (`zbasis`), coordinates (`xbasis`), and velocities (`vel_basis`).
*   **`m_ecp`**: Interacts to determine if an ECP-specific basis set should be loaded for an atom.
*   **`m_gaussian`**: Fundamental dependency. The `basis_function` type contains `gaussian` objects, and many routines (overlap, kinetic, etc.) call corresponding routines in `m_gaussian` for primitive Gaussians.
*   **`m_cart_to_pure`**: Used for transformations if `gaussian_type` is 'PURE'.
*   **Libcint Interface**: The `basis_set` type includes members (`LIBCINT_*`) to store data formatted for the Libcint library, suggesting an interaction for more complex (e.g., two-electron) integral calculations handled by Libcint.
```
