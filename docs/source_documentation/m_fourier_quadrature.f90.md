# `m_fourier_quadrature.f90`

## Overview

The `m_fourier_quadrature` Fortran module in MOLGW appears to be an experimental or alternative implementation for calculating one-electron integrals—specifically overlap, kinetic energy, and electron-nucleus attraction integrals. Instead of the standard approach of evaluating these integrals in real space, this module explores their computation via numerical quadrature in Fourier (momentum) space. The Fourier transform of Gaussian basis functions and the Coulomb potential have simpler analytical forms, which can sometimes be exploited. The routines seem designed for testing and comparison against real-space analytical results, as they often take `reference` integrals as input.

## Key Components

*   **`setup_overlap_fourier(basis_p, basis_t, reference)`**:
    *   **Purpose**: Calculates generalized overlap integrals of the form <&phi;<sub>a</sub>|e<sup>i**q**&middot;**r**</sup>|&phi;<sub>b</sub>>, where **q** is related to a `velocity` vector. This is achieved by calling `setup_gos_fourier`.
    *   The `reference` argument suggests it's used for comparing results with analytically computed overlap integrals (likely S<sub>ab</sub> when velocity is zero).

*   **`setup_nucleus_fourier(basis_p, basis_t, reference)`**:
    *   **Purpose**: Calculates electron-nucleus attraction integrals V<sub>ab</sub> = <&phi;<sub>a</sub>|&sum;<sub>A</sub> -Z<sub>A</sub>/|**r**-**R**<sub>A</sub>|&phi;<sub>b</sub>> using Fourier space quadrature.
    *   **Method**:
        1.  Sets up a 3D numerical grid in q-space (momentum space) using Log3 radial quadrature and Lebedev-Laikov angular quadrature.
        2.  For each q-point on the grid:
            *   Calculates the structure factor S(q) = -&sum;<sub>A</sub> Z<sub>A</sub> e<sup>i**q**&middot;**R**<sub>A</sub></sup>.
            *   Computes generalized overlap matrix elements <&phi;<sub>a</sub>|e<sup>-i(**q**+**v**)&middot;**r**</sup>|&phi;<sub>b</sub>> using `setup_gos_fourier` (where **v** is a velocity vector, potentially zero).
            *   Accumulates the electron-nucleus integral as &sum;<sub>q</sub> w<sub>q</sub> (4&pi;/q<sup>2</sup>) * S(q) * <&phi;<sub>a</sub>|e<sup>-i(**q**+**v**)&middot;**r**</sup>|&phi;<sub>b</sub>> / (2&pi;)<sup>3</sup>.

*   **`setup_kinetic_fourier(basis_p, basis_t, reference)`**:
    *   **Purpose**: Calculates kinetic energy integrals T<sub>ab</sub> = <&phi;<sub>a</sub>|-1/2 &nabla;<sup>2</sup>|&phi;<sub>b</sub>> using Fourier space quadrature.
    *   **Method**:
        1.  Uses the same q-space grid as `setup_nucleus_fourier`.
        2.  For each q-point:
            *   Computes the Fourier transform of basis functions &phi;<sub>b</sub>(q) and &phi;<sub>a</sub>(q-v).
            *   Accumulates the kinetic energy integral as &sum;<sub>q</sub> w<sub>q</sub> ( |**q**-**v**|<sup>2</sup> / 2 ) * &phi;<sub>a</sub><sup>*</sup>(q-v) &phi;<sub>b</sub>(q) * (2&pi;)<sup>3</sup> (the (2&pi;)<sup>3</sup> factor in the provided code seems to be a typo and should likely be 1/(2&pi;)<sup>3</sup> or handled differently based on FT conventions, but the code has it as a multiplier before the sum). The formula used in the code for accumulation is `0.5 * CONJG(ekin) * (2*pi)^3` after summing `weight * NORM2(qvec)**2 * basis_function_qmv * CONJG(basis_function_q)`. The `NORM2(qvec)**2` should likely be `NORM2(qmvvec)**2` if the kinetic energy operator acts on &phi;<sub>a</sub>. The factor `0.5` is correct.

*   **`setup_gos_fourier(basis_p, basis_t, qvec, gos_ao)`**:
    *   **Purpose**: Computes a matrix of generalized overlap elements (or Generalized Oscillator Strength form factors) G<sub>ab</sub>(q) = <&phi;<sub>a</sub>|e<sup>i**q**&middot;**r**</sup>|&phi;<sub>b</sub>> for a given momentum transfer vector `qvec`.
    *   It iterates over shells of `basis_p` and `basis_t`, calls `basis_function_gos` (from `m_basis_set`) for each pair of basis functions, which in turn likely evaluates the Fourier transform of the product of two Gaussians.

*   **`evaluate_gos_cart(ga, gb, qvec, gos_cart)`**:
    *   **Purpose**: A low-level helper function (marked as "buggy" in comments) intended to calculate <g<sub>a</sub>|e<sup>i**q**&middot;**r**</sup>|g<sub>b</sub>> for two *primitive* Cartesian Gaussian functions `ga` and `gb`.
    *   Contains an inner `auxiliary` function to compute 1D components of the integral. This seems to be an analytical or semi-analytical approach for the primitive integrals.

## Important Variables/Constants

*   **`basis_p`, `basis_t` (Type `basis_set`)**: Input basis sets, potentially for a "projectile" and "target" in scattering contexts, or simply left and right functions of an integral.
*   **`velocity(3)` (Real(dp))**: A velocity vector used to introduce a Galilean shift q &rarr; q &plusmn; v in the momentum space representation of one of the basis functions. If zero, standard integrals are computed.
*   **q-space grid**: `qlist(3, nq)` (coordinates), `wq(nq)` (weights), `n1` (angular points), `nqradial` (radial points). Generated using Log3 radial scheme and Lebedev angular scheme.
*   **`alpha` (Real(dp), Parameter)**: A parameter (e.g., 5.0 or 6.0) controlling the Log3 radial grid generation for q-space.

## Usage Examples

This module is likely experimental and for testing the feasibility or accuracy of Fourier-space quadrature for one-electron integrals. The `reference` argument in the main subroutines suggests they are intended to be compared against known correct values from traditional real-space methods.

```fortran
! Conceptual usage (assuming basis_projectile and basis_target are initialized)
USE m_fourier_quadrature
USE m_hamiltonian_onebody ! For reference real-space integrals
IMPLICIT NONE

TYPE(basis_set) :: basis_projectile, basis_target
REAL(DP), ALLOCATABLE :: overlap_ref(:,:), kinetic_ref(:,:), nucleus_ref(:,:)

! ... Initialize basis_projectile, basis_target ...

! Calculate reference integrals using standard real-space methods
ALLOCATE(overlap_ref(basis_projectile%nbf, basis_target%nbf))
CALL setup_overlap(basis_projectile, basis_target, overlap_ref) ! From m_hamiltonian_onebody

! Test Fourier space overlap calculation
CALL setup_overlap_fourier(basis_projectile, basis_target, overlap_ref)

! ... Similarly for kinetic and nucleus attraction integrals ...

DEALLOCATE(overlap_ref)
! ...
```

## Dependencies and Interactions

*   **`m_definitions`**: For `dp`, `pi`, `im`, and other constants.
*   **`m_warning`**: For `die`.
*   **`m_timing`**: For performance measurement.
*   **`m_atoms`**: For atomic coordinates (`xatom`) and nuclear charges (`zvalence`) used in `setup_nucleus_fourier`.
*   **`m_basis_set`**: For the `basis_set` type and for `basis_function_fourier` and `basis_function_gos` which compute the Fourier transform of basis functions or their products.
*   **`m_inputparam`**: (Potentially, if grid parameters or velocity were read from input, though not directly shown).
*   **`m_hamiltonian_onebody`**: (Implicitly) This module seems to be an alternative to the real-space methods in `m_hamiltonian_onebody`.
*   **Lebedev Quadrature (`lebedev_quadrature.f`)**: Uses `ld0230` for generating angular q-points.
*   **MPI (`m_mpi`)**: The q-space integration loop in `setup_nucleus_fourier` and `setup_kinetic_fourier` is parallelized using `world%nproc` and `world%rank`.
*   The module's correctness depends heavily on the accurate Fourier transformation of Gaussian basis functions (handled in `m_basis_set` via `m_gaussian`) and the robustness of the q-space quadrature grid. The comment "this routine is buggy" for `evaluate_gos_cart` suggests some parts might be incomplete or incorrect.
```
