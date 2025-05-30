# `m_io.f90`

## Overview

The `m_io` Fortran module in MOLGW serves as a comprehensive suite for handling various input and output (I/O) operations. This includes printing program headers and final summaries, reporting on memory usage and timings, managing warning messages, and providing specialized routines for outputting data in different formats. Key functionalities involve dumping matrices for debugging, generating Gaussian cube files for wavefunction and density visualization, creating WFN files for AIM analysis, reading density matrices or molecular orbitals from Gaussian fchk files, and logging trajectory-specific data for real-time TDDFT simulations.

## Key Components

*   **`header()` (Subroutine)**:
    *   Prints the initial MOLGW banner, including version information, compilation options, system date/time.
    *   Initializes and reports capabilities of linked libraries like Libxc, Libint/Libcint, OpenMP, ScaLAPACK, MKL, and HDF5.

*   **`this_is_the_end()` (Subroutine)**:
    *   The final routine called before program termination.
    *   Prints total memory usage statistics and detailed timing information for different computational sections.
    *   Outputs any accumulated warning messages.
    *   Writes run summary information (version, MPI/OpenMP setup, total time, peak memory) to the YAML output file if enabled.
    *   Prints a citation request.
    *   Finalizes MPI.

*   **`dump_out_matrix` (Generic Interface)**:
    *   Provides `dump_out_matrix_nospin_dp` (for 2D real(dp) arrays), `dump_out_matrix_dp` (for 3D real(dp) arrays, typically spin-resolved), and `dump_out_matrix_cdp` (for 3D complex(dp) arrays).
    *   Prints a segment (up to MAXSIZE x MAXSIZE) of a matrix to standard output, useful for debugging. Output is conditional on `print_matrix` flag or global `debug` flag.

*   **Population Analysis**:
    *   `mulliken_pdos(...)`, `mulliken_pdos_cmplx(...)`: Calculates and prints Mulliken population analysis, projecting MOs onto atomic orbital contributions (s, p, d, etc.) for selected atoms.
    *   `lowdin_pdos(...)`, `lowdin_pdos_cmplx(...)`: Calculates and prints Lowdin population analysis.

*   **Wavefunction and Density Plotting**:
    *   `plot_wfn(basis, c_matrix)`: Plots selected MOs along a line (defined in `manual_plotwfn` or default). Outputs to `fort.101` (phi) and `fort.102` (phi^2).
    *   `plot_wfn_fourier(basis, c_matrix)`: Computes and writes the Fourier transform of selected MOs to files (`wfn_fourier_xxxx.dat`).
    *   `plot_rho(rootname, basis, occupation, c_matrix)`: Plots the total electron density along a line (defined in `manual_plotrho` or default). Outputs to `[rootname]_density_cut.dat`.
    *   `plot_rho_xy(...)`: Generates a 2D plot of electron density integrated along z, in the xy-plane. Outputs to `density_plane_xy.dat`.
    *   `plot_cube_wfn(rootname, basis, occupation, c_matrix)`: Generates Gaussian cube files for specified MOs and the total electron density for each spin channel. Cube file parameters (grid dimensions, range) are controlled by input parameters.
    *   `plot_cube_wfn_cmplx(...)`: Similar to `plot_cube_wfn` but for complex-valued MO coefficients.
    *   `plot_cube_diff_cmplx(...)`: Plots the difference between the current density and an initial density (from `cube_density_start`) to a cube file. Used in RT-TDDFT.
    *   `calc_rho_initial_cmplx(...)`, `initialize_rho_diff_cmplx(...)`: Helper routines for `plot_cube_diff_cmplx` to store/manage the initial density.

*   **Trajectory-Specific Density Output (for RT-TDDFT)**:
    *   `plot_rho_traj_bunch(...)`: Plots density integrated along lines for different impact parameters.
    *   `plot_rho_traj_bunch_contrib(...)`: Similar to above but with contributions from selected states.
    *   `plot_rho_traj_points_set_contrib(...)`: Plots density contributions along lines defined by sets of start/end points from `manual_dens_points_set`.
    *   `calc_density_in_disc_cmplx_regular(...)`: Calculates integrated electron density within discs along z-axis, outputting to `disc_dens_xxxx_sY_...dat`.

*   **File I/O for External Formats**:
    *   `print_wfn_file(rootname, basis, occupation, c_matrix, etotal, energy, print_all)`: Generates a file in the WFN format (used by AIMPAC/AIMALL). Includes basis set details, MO coefficients, occupation numbers, and energies.
    *   `read_gaussian_fchk(read_fchk_in, file_name, basis, p_matrix_out)`: Reads a density matrix from a Gaussian formatted checkpoint (fchk) file. Handles reordering of Cartesian basis functions from Gaussian's convention to MOLGW's internal (Libint-like) convention.
    *   `read_guess_fchk(...)`: A specialized version of `read_gaussian_fchk` to read MO coefficients and optionally orbital energies as an initial guess for an SCF calculation.

*   **Quasiparticle Energy I/O**:
    *   `write_energy_qp(energy_qp)`: Writes quasiparticle energies to a file named `ENERGY_QP`.
    *   `read_energy_qp(nstate, energy_qp, reading_status)`: Reads quasiparticle energies from `ENERGY_QP` or `energy_qp`.

*   **Wavefunction Utilities**:
    *   `evaluate_wfn_r(...)`: Evaluates a range of MOs at a specific point `rr` in real space.
    *   `wfn_parity(...)`, `wfn_reflection(...)`: Determines the parity or reflection symmetry of a given MO with respect to the molecular center or a plane.

*   **HDF5 Output**:
    *   `dump_matrix_cmplx_hdf5(...)`: Writes a complex 3D matrix to an HDF5 file as two datasets (real and imaginary parts).
    *   `print_restart_hdf5(...)`: Saves key information (basis set labels, atomic coordinates, overlap matrix, MO coefficients, energies) to an HDF5 restart file (`molgw_restart.h5`).

*   **YAML Output**:
    *   `calculation_parameters_yaml(...)`: Writes key calculation parameters (number of basis functions, states, etc.) to the YAML output file.

## Important Variables/Constants

*   **`stdout`, `unit_yaml` (from `m_definitions`)**: Standard output and YAML file unit numbers.
*   **Input Parameters (from `m_inputparam`)**: Many routines are controlled by boolean flags from `m_inputparam` (e.g., `print_cube_`, `print_wfn_files_`, `print_yaml_`). Cube file generation is also controlled by `cube_state_min`, `cube_state_max`, `cube_nx`, `cube_ny`, `cube_nz`.
*   **`MOLGW_VERSION` (Preprocessor Constant)**: Used in `header()`.

## Usage Examples

Most subroutines in `m_io` are called internally by MOLGW based on user input or the stage of calculation.
*   `header()` is called at the very beginning of a MOLGW run.
*   `this_is_the_end()` is called at the very end.
*   `plot_cube_wfn` is called if `print_cube = YES` is in the input file.
*   `read_gaussian_fchk` is called if `read_fchk = 'path/to/file.fchk'` is specified.

## Dependencies and Interactions

*   **Core Modules**: `m_definitions`, `m_mpi`, `m_timing`, `m_memory`, `m_warning`, `m_string_tools`, `m_inputparam`.
*   **System Description**: `m_atoms`, `m_basis_set`, `m_elements`, `m_cart_to_pure`.
*   **Hamiltonian/SCF**: `m_hamiltonian_tools`, `m_scf`.
*   **DFT Grid**: `m_dft_grid` (for `calculate_basis_functions_r` used in plotting routines).
*   **Integral Libraries**: `m_libint_tools`, `m_libcint_tools` (for `libint_init` and `check_capability_libcint` called in `header`).
*   **Libxc**: `m_libxc_tools` (for `xc_version` called in `header`).
*   **HDF5**: `m_hdf5_tools` (if `HAVE_HDF5` is defined, for HDF5 output routines).
*   **File System**: Interacts extensively with the file system for reading various input files (e.g., `manual_plotwfn`, `manual_pdos`, fchk files, `ENERGY_QP`) and writing output files (cube files, WFN files, density cuts, YAML logs, HDF5 restarts).
*   This module centralizes many of the user-facing output generation tasks and some specialized input tasks.
```
