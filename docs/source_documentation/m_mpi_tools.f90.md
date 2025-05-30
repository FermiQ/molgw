# `m_mpi_tools.f90`

## Overview

The `m_mpi_tools` Fortran module in MOLGW provides a high-level, object-oriented wrapper for common Message Passing Interface (MPI) operations. It defines a derived type `mpi_communicator` which encapsulates an MPI communicator and provides methods for typical parallel tasks such as barriers, collective reductions (sum, min, max, logical AND), and broadcasts. The use of MPI calls is conditionally compiled based on the `HAVE_MPI` preprocessor macro, allowing the code to run in serial mode if MPI is not enabled. This module simplifies parallel programming constructs in other parts of MOLGW.

## Key Components

*   **`mpi_communicator` (Type)**:
    *   A derived data type that holds information about a specific MPI communicator.
    *   **Members**:
        *   `comm` (Integer): The MPI communicator handle (e.g., `MPI_COMM_WORLD` or a sub-communicator).
        *   `nproc` (Integer): The total number of processes within this communicator.
        *   `rank` (Integer): The rank (0 to `nproc`-1) of the current process within this communicator.
    *   **Type-Bound Procedures (Methods)**:
        *   **`init(comm_in)`**: Initializes the `mpi_communicator` object using an existing MPI communicator handle `comm_in`. It queries and stores `nproc` and `rank`. In serial mode, sets `nproc=1` and `rank=0`.
        *   **`barrier()`**: Synchronizes all processes within the communicator by calling `MPI_BARRIER`.
        *   **`sum(array)` (Generic)**: Performs an in-place global sum reduction on `array` across all processes in the communicator using `MPI_ALLREDUCE` with `MPI_SUM`. Overloaded for `REAL(dp)`, `COMPLEX(dp)`, and `INTEGER(kind=int8)` arrays of any rank.
        *   **`min(array)` (Generic)**: Performs an in-place global minimum reduction using `MPI_ALLREDUCE` with `MPI_MIN`. Overloaded for `REAL(dp)` and `INTEGER` arrays.
        *   **`max(array)` (Generic)**: Performs an in-place global maximum reduction using `MPI_ALLREDUCE` with `MPI_MAX`. Overloaded for `REAL(dp)` and `INTEGER` arrays.
        *   **`and(array)`**: Performs an in-place global logical AND reduction on a `LOGICAL` array using `MPI_ALLREDUCE` with `MPI_LAND`.
        *   **`bcast(rank, array)` (Generic)**: Broadcasts `array` from the process with the specified `rank` to all other processes in the communicator using `MPI_BCAST`. Overloaded for `INTEGER`, `REAL(dp)`, and `COMPLEX(dp)` arrays.

*   **Implementation Subroutines (e.g., `mpic_sum_dp`)**:
    *   These are private module subroutines that implement the actual logic for each type-bound procedure. They contain the conditional MPI calls and handle the `MPI_IN_PLACE` mechanism for reductions. They bypass MPI calls if `mpic%nproc == 1`.

## Important Variables/Constants

This module does not define public constants for general use. Its primary export is the `mpi_communicator` type. Internal MPI constants like `MPI_SUM`, `MPI_DOUBLE_PRECISION`, etc., are used from the standard `mpi` module.

## Usage Examples

```fortran
MODULE my_parallel_routine
  USE m_mpi_tools
  USE m_definitions, ONLY: dp
  IMPLICIT NONE

  SUBROUTINE do_parallel_sum(world_mpi_comm_handle)
    INTEGER :: world_mpi_comm_handle ! Assume this is a valid MPI communicator handle
    TYPE(mpi_communicator) :: world
    REAL(dp) :: local_data, sum_data
    INTEGER :: ierr

    CALL world%init(world_mpi_comm_handle)

    local_data = REAL(world%rank + 1, dp) ! Each process has different data
    sum_data = local_data

    CALL world%sum(sum_data) ! Perform global sum

    IF (world%rank == 0) THEN
      PRINT *, "Global sum of ranks + 1 is: ", sum_data
    END IF

    CALL world%barrier()
    
    ! Example of broadcast
    IF (world%rank == 0) THEN
        local_data = 123.456_dp
    ELSE
        local_data = 0.0_dp
    END IF
    CALL world%bcast(0, local_data) ! Broadcast from rank 0
    ! Now all processes have local_data = 123.456_dp
    
  END SUBROUTINE do_parallel_sum

END MODULE my_parallel_routine

PROGRAM main
  USE mpi            ! For MPI_INIT, MPI_FINALIZE, MPI_COMM_WORLD
  USE my_parallel_routine
  IMPLICIT NONE
  INTEGER :: ierr
  
  CALL MPI_INIT(ierr)
  CALL do_parallel_sum(MPI_COMM_WORLD)
  CALL MPI_FINALIZE(ierr)
  
END PROGRAM main
```

## Dependencies and Interactions

*   **MPI Library**: The module's parallel functionalities are entirely dependent on an MPI library. All MPI-specific calls are enclosed in `#if defined(HAVE_MPI) ... #endif` blocks, allowing the code to compile and run in a serial (non-MPI) environment where these calls are skipped. It `USE mpi` for MPI constants and procedure interfaces.
*   **`m_definitions`**: For precision kinds like `dp` (double precision) and `int8`.
*   **`m_warning`**: For the `die` subroutine (though not directly used in the provided snippet, it's a common utility in MOLGW for error handling).
*   **Other MOLGW Modules**: This module provides a convenient abstraction for MPI communication that can be used by any other module in MOLGW requiring parallel data exchange or synchronization. It simplifies common collective operations and broadcasts.
```
