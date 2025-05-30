# `m_lbfgs.f90`

## Overview

The `m_lbfgs` Fortran module implements the Limited-memory Broyden–Fletcher–Goldfarb–Shanno (L-BFGS) algorithm. L-BFGS is a widely-used quasi-Newton method for solving unconstrained nonlinear optimization problems. It approximates the Hessian matrix using information from a limited number of previous gradient evaluations and steps, making it suitable for large-scale problems where computing or storing the true Hessian is impractical. This module provides routines to initialize the optimizer, execute optimization steps, and manage its state. It includes a line search algorithm (More-Thuente) to satisfy Wolfe conditions.

## Key Components

*   **`lbfgs_state` (Type)**:
    *   A derived data type that encapsulates the state of the L-BFGS optimizer. This allows multiple independent L-BFGS optimizations to be carried out if needed.
    *   **Members**:
        *   `lbfgs_status` (Integer): Current status of the optimizer (e.g., initial call, converged, error).
        *   `ndim` (Integer): Dimensionality of the optimization problem (number of variables).
        *   `history_record` (Integer): The number of past steps and gradient differences to store (parameter M in L-BFGS).
        *   `iter` (Integer): Current iteration count.
        *   `gtol` (Real(dp)): Tolerance for the projected gradient in the line search (related to Wolfe conditions).
        *   `diag(:)` (Real(dp), Allocatable): Diagonal approximation of the inverse Hessian matrix.
        *   `work(:)` (Real(dp), Allocatable): Workspace array used by L-BFGS to store past steps (s<sub>k</sub>) and gradient differences (y<sub>k</sub>), and other intermediate values. Its size is `ndim * (2 * history_record + 1) + 2 * history_record`.
        *   `line_*` variables: A collection of scalar and logical variables used to maintain the state of the More-Thuente line search algorithm (`mcsrch`) between calls. This includes current step length (`line_stp`), step bounds (`line_stpmin`, `line_stpmax`), function and derivative values at trial points, etc.

*   **`lbfgs_init(lbfgs_plan, ndim, history_record, diag_guess)` (Subroutine)**:
    *   **Purpose**: Initializes an `lbfgs_state` object (`lbfgs_plan`).
    *   Sets dimensionality, history size, allocates workspace arrays (`diag`, `work`), and initializes line search parameters and default diagonal Hessian guess.

*   **`lbfgs_destroy(lbfgs_plan)` (Subroutine)**:
    *   **Purpose**: Deallocates the workspace arrays (`diag`, `work`) within an `lbfgs_state` object.

*   **`lbfgs_execute(lbfgs_plan, x, f, gradf)` (Function)**:
    *   **Purpose**: Executes one iteration of the L-BFGS optimization algorithm.
    *   Takes the current variable vector `x`, function value `f`, and gradient `gradf` as input.
    *   Internally calls the `lbfgs` subroutine.
    *   Updates `x` to the new set of variables determined by the L-BFGS step and line search.
    *   Returns an integer status code (`lbfgs_plan%lbfgs_status`) indicating the outcome (e.g., 0 for convergence, 1 for ongoing line search, negative for error).

*   **`lbfgs(...)` (Subroutine, Private)**:
    *   **Purpose**: The core L-BFGS algorithm logic.
    *   If `IFLAG` indicates a new evaluation is needed by the line search, it calls `MCSRCH`.
    *   Updates the L-BFGS history (differences in `x` and `gradient`) and the diagonal Hessian approximation.
    *   Computes the new search direction -H*g using the L-BFGS two-loop recursion.
    *   Calls `MCSRCH` to perform a line search along the new direction to find a step `STP` satisfying Wolfe conditions.
    *   Updates `X` (coordinates) based on the accepted step.
    *   Sets `IFLAG` to manage control flow (e.g., request new function/gradient evaluation if line search is ongoing).

*   **`mcsrch(...)` (Subroutine, Private)**:
    *   **Purpose**: Implements the More-Thuente line search algorithm to find a step length `STP` that satisfies the strong Wolfe conditions (sufficient decrease in function value and sufficient decrease in projected gradient magnitude).
    *   It iteratively calls `mcstep` to refine an interval bracketing a suitable step length.
    *   Communicates with the calling routine (`lbfgs`) via `INFO` flag to indicate convergence, errors, or if a new function/gradient evaluation is needed at the trial `STP`.

*   **`mcstep(...)` (Subroutine, Private)**:
    *   **Purpose**: A helper routine for `mcsrch`. Given an interval [STX, STY] and function/derivative values at these points, it computes a new trial step length `STP` using interpolation (cubic or secant method).
    *   Updates the bracketing interval [STX, STY] based on the function and derivative values at the new `STP`.

## Important Variables/Constants

*   **`lbfgs_plan%history_record` (Integer)**: The number of past updates (s<sub>k</sub>, y<sub>k</sub> pairs) stored and used to approximate the inverse Hessian. Typical values are 5-20.
*   **`lbfgs_plan%gtol` (Real(dp))**: A parameter for the line search, typically related to the curvature condition in Wolfe criteria (e.g., |g<sub>k+1</sub><sup>T</sup>s<sub>k</sub>| &le; gtol * |g<sub>k</sub><sup>T</sup>s<sub>k</sub>|). Default is 0.9.
*   **`FTOL` (Real(dp), Parameter in `mcsrch`)**: Tolerance for the sufficient decrease (Armijo) condition in the line search (f(x<sub>k</sub> + &alpha;<sub>k</sub>p<sub>k</sub>) &le; f(x<sub>k</sub>) + c<sub>1</sub>&alpha;<sub>k</sub>p<sub>k</sub><sup>T</sup>&nabla;f<sub>k</sub>). Default is 1.0e-4.
*   **`XTOL` (Real(dp), Parameter in `mcsrch`)**: Relative tolerance for the width of the bracketing interval in the line search. Default is 1.0e-17.
*   **`STPMIN`, `STPMAX` (Real(dp))**: Minimum and maximum allowable step lengths during the line search.

## Usage Examples

This module is typically used for unconstrained optimization problems like molecular geometry optimization. The `lbfgs_execute` function is called iteratively.

```fortran
! Conceptual usage within a geometry optimization loop:
USE m_lbfgs
USE m_definitions
IMPLICIT NONE

INTEGER, PARAMETER :: N_DIM = 100 ! Number of variables to optimize
INTEGER, PARAMETER :: M_HISTORY = 7 ! L-BFGS history size
TYPE(lbfgs_state) :: opt_state
REAL(DP) :: coords(N_DIM), energy_val, grad_val(N_DIM)
INTEGER :: iter_count, lbfgs_status_flag

! Initialize L-BFGS state
CALL lbfgs_init(opt_state, N_DIM, M_HISTORY)
! Optionally, provide an initial guess for the diagonal Hessian (e.g., 1.0 / approximate force constants)
! CALL lbfgs_init(opt_state, N_DIM, M_HISTORY, diag_guess=0.01_dp) 

! ... Initialize coords with starting geometry ...

DO iter_count = 1, MAX_ITERATIONS
  ! 1. Calculate function value (energy_val) and gradient (grad_val) 
  !    at the current 'coords'. The L-BFGS routine expects the gradient
  !    of the function to be minimized (e.g., -Force for energy minimization).
  CALL calculate_my_function_and_gradient(coords, energy_val, grad_val)

  ! 2. Execute L-BFGS step. Note: grad_val is modified by lbfgs.
  !    The lbfgs_execute wrapper handles the IFLAG logic of the underlying lbfgs routine.
  !    coords array is updated in place by lbfgs_execute.
  lbfgs_status_flag = lbfgs_execute(opt_state, coords, energy_val, grad_val)

  ! 3. Check convergence status
  IF (lbfgs_status_flag == 0) THEN
    PRINT *, "L-BFGS optimization converged."
    EXIT
  ELSE IF (lbfgs_status_flag < 0) THEN
    PRINT *, "L-BFGS optimization error: ", lbfgs_status_flag
    EXIT
  END IF
  ! If lbfgs_status_flag == 1, the L-BFGS routine (specifically mcsrch) 
  ! has updated coords to a new trial point and expects a new function
  ! and gradient evaluation at this new point for the next call.
  ! The lbfgs_execute wrapper manages this internal state.
END DO

CALL lbfgs_destroy(opt_state)

CONTAINS
  SUBROUTINE calculate_my_function_and_gradient(x_vars, f_val, g_val)
    REAL(DP), INTENT(IN)    :: x_vars(N_DIM)
    REAL(DP), INTENT(OUT)   :: f_val
    REAL(DP), INTENT(OUT)   :: g_val(N_DIM)
    ! ... user code to compute f_val and g_val at x_vars ...
  END SUBROUTINE
```

## Dependencies and Interactions

*   **`m_definitions`**: For the `dp` (double precision) kind parameter.
*   **Numerical Algorithms**: The L-BFGS method and the More-Thuente line search are standard numerical optimization algorithms. The implementation details suggest it might be based on or inspired by well-established Fortran codes for these algorithms (e.g., from Nocedal, or MINPACK).
*   **Calling Modules**: This module is designed to be a general-purpose optimizer. In MOLGW, it's notably used by `m_atoms` (specifically the `relax_atoms` subroutine) for geometry optimization, where the "function" is the total energy and the "gradient" is the negative of the atomic forces.
```
