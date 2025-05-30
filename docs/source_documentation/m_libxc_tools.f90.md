# `m_libxc_tools.f90`

## Overview

The `m_libxc_tools` Fortran module serves as the primary interface between MOLGW and the Libxc library. Libxc is a comprehensive C library providing a vast collection of exchange-correlation (XC) functionals used in Density Functional Theory (DFT). This module defines a Fortran derived type (`dft_xc_info`) to encapsulate information about specific XC functionals and includes the necessary interfaces to call Libxc's C API. These calls allow MOLGW to initialize, query, and evaluate XC energies, potentials, and kernels from the Libxc library. The functionality of this module is conditional on MOLGW being compiled and linked against Libxc (i.e., `NO_LIBXC` preprocessor macro not defined).

## Key Components

*   **Libxc Functional Family Parameters (Integer, Parameters)**:
    *   Named constants like `XC_FAMILY_LDA`, `XC_FAMILY_GGA`, `XC_FAMILY_MGGA`, `XC_FAMILY_HYB_GGA`, etc., which mirror Libxc's internal classification of functional families.

*   **`dft_xc_info` (Type)**:
    *   A Fortran derived data type to store information and the runtime state for one or more components of an XC functional.
    *   **Members**:
        *   `needs_gradient` (Logical): Set to `.TRUE.` if any functional component requires density gradients (i.e., is a GGA, meta-GGA, or hybrid GGA/meta-GGA).
        *   `nxc` (Integer): The number of actual XC functional components defined in this instance (if `dft_xc_info` is used as an array to represent composite functionals).
        *   `id` (Integer(C_INT)): The numerical identifier of the XC functional (or component) as defined in Libxc.
        *   `coeff` (Real(dp)): The coefficient multiplying this functional component in a composite functional (e.g., in B3LYP or hybrid functionals).
        *   `gamma` (Real(C_DOUBLE)): The range-separation parameter &omega; (in Bohr<sup>-1</sup>) for range-separated hybrid functionals or other functionals that use it.
        *   `nspin` (Integer(C_INT)): Number of spin channels (1 for restricted/spin-unpolarized, 2 for unrestricted/spin-polarized).
        *   `family` (Integer): The Libxc family code for this functional (e.g., `XC_FAMILY_LDA`).
        *   `func` (Type(C_PTR), Pointer): A C pointer to the initialized Libxc functional object (`xc_func_type*` in C).

*   **C Function Interfaces (`INTERFACE` block)**:
    *   Provides Fortran interfaces to a wide range of C functions from the Libxc library. Key interfaces include:
        *   `xc_func_alloc()`: Allocates a Libxc functional object.
        *   `xc_func_free(func)`: Frees a Libxc functional object.
        *   `xc_version(major, minor, micro)`: Retrieves the Libxc library version.
        *   `xc_func_init(func, functional_id, nspin)`: Initializes a functional object for a given Libxc `functional_id` and number of spin channels.
        *   `xc_func_end(func)`: Finalizes (deinitializes) a functional object.
        *   `xc_func_get_info(func)`: Gets an info object for an initialized functional.
        *   `xc_func_info_get_name(info)`, `xc_func_info_get_family(info)`, `xc_func_info_get_flags(info)`: Retrieve specific details about a functional from its info object.
        *   `xc_func_set_ext_params_name(func, name, param)`: Sets external parameters (like `_omega`) for certain functionals.
        *   **LDA Evaluation Routines**: `xc_lda_exc` (energy), `xc_lda_exc_vxc` (energy and potential), `xc_lda_fxc` (kernel).
        *   **GGA Evaluation Routines**: `xc_gga_exc` (energy), `xc_gga_exc_vxc` (energy and potential), `xc_gga_vxc` (potential), `xc_gga_fxc` (kernel).
        *   `xc_hyb_cam_coef(func, omega, alpha, beta)`: Retrieves CAM-type hybrid parameters (&omega;, &alpha; for HF exchange, &beta; for HF long-range exchange) from a Libxc functional object.
        *   `xc_functional_get_number(name)`, `xc_functional_get_name(number)`: Convert between Libxc functional names and their numerical IDs.

*   **`init_libxc_info(dft_xc)` (Subroutine)**:
    *   **Purpose**: Initializes an array of `dft_xc_info` structures. For each entry with a non-zero `id`:
        1.  Allocates a Libxc functional object via `xc_func_alloc`.
        2.  Initializes this object using `xc_func_init` with the given Libxc `id` and `nspin`.
        3.  If the functional is a range-separated type (e.g., HSE, wPBEH, ITYH), it calls `xc_func_set_ext_params_name` to set the range-separation parameter `gamma` (&omega;).
        4.  Determines the `family` of the functional and sets the `needs_gradient` flag accordingly.

*   **`copy_libxc_info(dft_xc_in, dft_xc_out)` (Subroutine)**:
    *   Performs a deep copy of an array of `dft_xc_info` structures, ensuring that the Libxc functional pointers (`func`) are also correctly handled (though in this implementation, it seems to be a shallow copy of the pointer itself, which is typical if the lifetime of the pointed-to object is managed elsewhere or if it's intended to point to the same Libxc object).

*   **`destroy_libxc_info(dft_xc)` (Subroutine)**:
    *   Finalizes and deallocates the Libxc functional objects associated with an array of `dft_xc_info` structures by calling `xc_func_end` and `xc_func_free` for each.

## Important Variables/Constants

*   **`XC_FAMILY_*` (Integer, Parameters)**: Constants representing Libxc functional families (LDA, GGA, etc.).
*   **`XC_FLAGS_HAVE_FXC` (Integer, Parameter)**: Libxc flag indicating if a functional provides the fxc kernel (second derivative).
*   **`dft_xc_info%func` (Type(C_PTR))**: The C pointer to the Libxc functional object, which is the handle used for all subsequent calls to Libxc evaluation routines for that functional.

## Usage Examples

This module is primarily used internally by other MOLGW modules, particularly `m_inputparam` (for parsing functional choices) and `m_hamiltonian_twobodies` (for evaluating XC contributions).

Conceptual setup in `m_inputparam`:
```fortran
USE m_libxc_tools
USE m_definitions
IMPLICIT NONE

TYPE(dft_xc_info), ALLOCATABLE :: dft_settings(:)
CHARACTER(LEN=20) :: user_functional_key = "PBE" ! Example from user input

! ... (logic to parse 'user_functional_key' into Libxc IDs) ...
! Assume PBE = XC_GGA_X_PBE (ID 101) + XC_GGA_C_PBE (ID 131)

ALLOCATE(dft_settings(2))
dft_settings(1)%id = 101 ! XC_GGA_X_PBE
dft_settings(1)%nspin = 1 ! Assuming restricted calculation
dft_settings(1)%coeff = 1.0_dp
dft_settings(2)%id = 131 ! XC_GGA_C_PBE
dft_settings(2)%nspin = 1
dft_settings(2)%coeff = 1.0_dp

CALL init_libxc_info(dft_settings) 
! Now dft_settings(1)%func and dft_settings(2)%func are initialized Libxc handles.
! These can be passed to m_hamiltonian_twobodies::dft_exc_vxc_batch.

! ... (at the end of the program) ...
CALL destroy_libxc_info(dft_settings)
```

## Dependencies and Interactions

*   **Libxc C Library**: This module is critically dependent on the external Libxc library. All `xc_*` functions are C functions provided by Libxc. The module's functionality is enabled by the `!defined(NO_LIBXC)` preprocessor directive. The C header files `xc_funcs.h` and `xc_version.h` are included via the C preprocessor.
*   **`m_definitions`**: For `dp` (double precision kind) and C interoperability types (`C_INT`, `C_DOUBLE`, `C_PTR`, `C_CHAR`).
*   **`m_warning`**: For the `die` error handling subroutine.
*   **`m_inputparam`**: Uses `init_libxc_info` to process user-specified DFT functional choices and populate the `dft_xc` array, which is then used globally.
*   **`m_hamiltonian_twobodies`**: The `dft_exc_vxc_batch` subroutine in `m_hamiltonian_twobodies` takes the initialized `dft_xc_info` array as input to evaluate XC energy densities and potential values on the DFT grid by calling the Libxc routines like `xc_lda_exc_vxc` or `xc_gga_exc_vxc`.
```
