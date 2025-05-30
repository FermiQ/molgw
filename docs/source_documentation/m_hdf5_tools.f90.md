# `m_hdf5_tools.f90`

## Overview

The `m_hdf5_tools` module provides a set of high-level Fortran wrapper subroutines for interacting with HDF5 (Hierarchical Data Format 5) files. It aims to simplify common HDF5 operations, such as opening and closing files, creating and managing groups, and reading and writing datasets and attributes. The module supports various data types including integers, double-precision reals, and strings, for both scalar and multi-dimensional arrays (up to 6D for datasets). The code is adapted from Justin Erwin's `HDF5_utils` library and is conditional on the HDF5 library being available during compilation (`HAVE_HDF5`).

## Key Components

The module uses Fortran generic interfaces to provide a consistent API for different data types and ranks.

*   **Initialization and Finalization**:
    *   **`hdf_init()`**: Initializes the HDF5 Fortran library. Should be called before other HDF5 operations.
    *   **`hdf_finalize()`**: Closes the HDF5 Fortran library. Should be called at the end of all HDF5 operations.

*   **File Operations**:
    *   **`hdf_open_file(file_id, filename, STATUS, ACTION)`**: Opens an existing HDF5 file or creates a new one.
        *   `file_id` (Output, `HID_T`): HDF5 identifier for the opened file.
        *   `filename` (Input, Character): Name of the HDF5 file.
        *   `STATUS` (Input, Optional, Character): File status ('OLD', 'NEW', 'REPLACE'). Default is 'NEW'.
        *   `ACTION` (Input, Optional, Character): File access mode ('READ', 'WRITE', 'READWRITE'). Default is 'READWRITE'.
    *   **`hdf_close_file(file_id)`**: Closes an HDF5 file identified by `file_id`.

*   **Group Operations**:
    *   **`hdf_create_group(loc_id, group_name)`**: Creates a new group with `group_name` at the location specified by `loc_id` (file or group identifier).
    *   **`hdf_open_group(loc_id, group_name, group_id)`**: Opens an existing group.
    *   **`hdf_close_group(group_id)`**: Closes an HDF5 group.

*   **Dataset Query Operations**:
    *   **`hdf_exists(loc_id, obj_name)` (Function)**: Checks if an object (dataset, group, or link) with `obj_name` exists at `loc_id`. Returns `.TRUE.` or `.FALSE.`.
    *   **`hdf_get_rank(loc_id, dset_name, rank)`**: Retrieves the rank (number of dimensions) of the dataset `dset_name`.
    *   **`hdf_get_dims(loc_id, dset_name, dims)`**: Retrieves the dimensions of the dataset `dset_name` into the integer array `dims`.

*   **Dataset I/O Operations**:
    *   **`hdf_create_dataset(loc_id, dset_name, dset_dims, dset_type)`**: Creates a new dataset named `dset_name` at `loc_id` with dimensions `dset_dims` and data type `dset_type` (e.g., 'integer', 'double').
    *   **`hdf_write_dataset(loc_id, dset_name, data)` (Generic Interface)**: Writes data (scalar up to 6D arrays of integer, double, or string) to a new dataset.
    *   **`hdf_read_dataset(loc_id, dset_name, data)` (Generic Interface)**: Reads data from an existing dataset into the provided `data` array/scalar.
    *   **`hdf_update_dataset(loc_id, dset_name, data)` (Generic Interface, e.g., `hdf_update_dataset_integer_0`)**: Updates an existing dataset with new data. Currently implemented for scalar integers.
    *   **`hdf_write_vector_to_dataset(loc_id, dset_name, offset, vector)` (Generic Interface)**: Writes a 1D `vector` into a specified `offset` (hyperslab) within an existing multi-dimensional dataset.
    *   **`hdf_read_vector_from_dataset(loc_id, dset_name, offset, vector)` (Generic Interface)**: Reads a 1D `vector` from a specified `offset` within an existing dataset.

*   **Attribute I/O Operations**:
    *   **`hdf_write_attribute(loc_id, obj_name, attr_name, data)` (Generic Interface)**: Writes an attribute (`data`) with name `attr_name` to an object (`obj_name`) located at `loc_id`. Supports scalar/1D integer, double, and scalar string attributes. If `obj_name` is empty, attaches to `loc_id`.
    *   **`hdf_read_attribute(loc_id, obj_name, attr_name, data)` (Generic Interface)**: Reads an attribute.

*   **Utility**:
    *   **`hdf_set_print_messages(val_print_messages)`**: Toggles verbose diagnostic messages from the module.

## Important Variables/Constants

*   **`HID_T` (Integer, Parameter)**: Publicly exposed if `HAVE_HDF5` is defined. Represents the HDF5 identifier type (typically a kind of integer, often 8-byte for modern HDF5 Fortran interfaces).
*   **`hdf_print_messages` (Logical, Private)**: Module flag controlling verbose output. Default is `.FALSE.`.
*   The module internally uses HDF5 library constants like `H5F_ACC_RDONLY_F`, `H5T_NATIVE_DOUBLE`, `H5S_SCALAR_F`, etc., which are provided by the `USE hdf5` statement.

## Usage Examples

```fortran
USE m_hdf5_tools
USE m_definitions, ONLY: dp ! For real(dp) kind
IMPLICIT NONE

INTEGER(HID_T) :: file_handle, group_handle
REAL(dp) :: array_2d(10, 5), array_read_back(10, 5)
INTEGER :: scalar_int_attr
CHARACTER(LEN=50) :: string_attr

CALL hdf_init()
CALL hdf_set_print_messages(.TRUE.) ! Enable verbose messages

! Create and open a new HDF5 file
CALL hdf_open_file(file_handle, "example.h5", STATUS="REPLACE")

! Create a group
CALL hdf_create_group(file_handle, "MyDataGroup")
CALL hdf_open_group(file_handle, "MyDataGroup", group_handle)

! Prepare some data
array_2d = 1.0_dp 

! Write a 2D double-precision dataset
CALL hdf_write_dataset(group_handle, "DoublePrecisionArray", array_2d)

! Write a scalar integer attribute to the dataset
scalar_int_attr = 123
CALL hdf_write_attribute(group_handle, "DoublePrecisionArray", "Version", scalar_int_attr)

! Write a string attribute to the group
string_attr = "This is a test group"
CALL hdf_write_attribute(group_handle, "", "Description", string_attr) ! obj_name="" attaches to group_handle

! Read back the dataset
CALL hdf_read_dataset(group_handle, "DoublePrecisionArray", array_read_back)

! Check if data matches (simple check)
IF (ABS(SUM(array_2d - array_read_back)) < 1e-9_dp) THEN
  PRINT *, "HDF5 Read/Write successful!"
ELSE
  PRINT *, "HDF5 Read/Write mismatch!"
END IF

CALL hdf_close_group(group_handle)
CALL hdf_close_file(file_handle)

CALL hdf_finalize()
```

## Dependencies and Interactions

*   **HDF5 Library**: This module is a wrapper around the HDF5 Fortran library. Its entire functionality is conditional on the `HAVE_HDF5` preprocessor macro being defined and linking against the HDF5 library. It uses the standard `hdf5` Fortran module.
*   **`m_definitions`**: For the `dp` (double precision) kind parameter.
*   **Other MOLGW Modules**: Any module within MOLGW that requires persistent storage of large or complex data (e.g., wavefunctions, density matrices, ERI tensors if out-of-core, restart information) might use `m_hdf5_tools` for these I/O operations.
```
