
FEISTY library:
 * `FEISTY.F90`: The core routines, incl. wrapper for R.
 * `setup.F90`: FEISTY setups and dependent subroutines.
 * `globals.F90`: Global variables and parameters.
 * `fish.F90`: Initialize fish groups.
 * `spectrum.F90`: Fish size spectrum build up.
 * `input.F90`: Parameter input subroutines for R.
 * `FEISTYtest.F90`: Build up as a program for debugging.

R package required:
 * `Makevars`: This defines how Fortran source files are compiled into object files when building the R package on Mac and Linux.
 * `Makevars.win`: This defines how Fortran source files are compiled into object files when building the R package on Windows.
 * `R_init_feisty.c`: This registers Fortran routines with R, so they can be called from R.

FABM-FEISTY wrappers:
 * `FEISTY_FABM.F90`: The wrapper for the FABM system
 * `CMakeLists.txt`: CMake definitions for FABM
