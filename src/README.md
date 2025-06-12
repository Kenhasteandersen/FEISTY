Contains the fortran code for FEISTY

FEISTY library:
 * `FEISTY.F90`:  The core routines, incl. wrapper for R
 * `setup.F90`:
 * `globals.F90`: global variables
 * `spectrum.F90`: class definition of a size spectrum
 * `fish.F90`: Definition of the fish class
 * `input.F90`:
 * `FEISTYtest.F90`: 

R package required:
 * `Makevars`: This defines how Fortran source files are compiled into object files when building R package on Mac and Linux.
 * `Makevars.win`: This defines how Fortran source files are compiled into object files when building R package on Windows.
 * `R_init_feisty.c`: This registers Fortran routines with R, so they can be called from R.

FABM-FEISTY wrappers:
 * `FEISTY_FABM.F90`: The wrapper for the FABM system
 * `CMakeLists.txt`: CMake definitions for FABM
