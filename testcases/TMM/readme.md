# FABM-FEISTY in the Transport Matrix Method (TMM) through FABM-OS

Configuration files and run scripts for running FEISTY, coupled to ERSEM or PISCES, in the global TMM through FABM-OS. They are adapted from the GOTM test cases in [`../GOTM`](../GOTM/README.md), whose README also describes FABM-FEISTY and how its `fabm.yaml` section is set up. An earlier FEISTY-UVic configuration used an older FEISTY interface and has been retired.

> All three cases have been run for one month on `MITgcm_2.8deg` (Windows, gfortran, FABM-OS with FABM pinned as described in step 2 below). They have not yet been run long enough to check the results.

## Files

| FABM configuration | Run script | Case |
|--------------------|------------|------|
| `fabm_ersem_noben.yaml` | `run_ersem_noben.py` | ERSEM with simple benthic pools and FEISTY's internal benthos (labile deposition redirected into the `fish` instance) |
| `fabm_ersem_ben.yaml` | `run_ersem_ben.py` | ERSEM with the full benthic model; benthic fauna Y2, Y3 and Y4 are the benthic prey |
| `fabm_pisces.yaml` | `run_pisces.py` | PISCES (122:16:1 stoichiometry) with FEISTY's internal benthos |
| `setup.cfg` | | Template build configuration for FABM-OS (compiler paths and the model sources) |
| `gotm_feisty_workflow.Rmd` | | R workflow for working with FEISTY results |

## Installation

### 1. Get the components

| Component | Source |
|-----------|--------|
| FABM-OS | `git clone --recurse-submodules https://github.com/BoldingBruggeman/fabmos.git` (the release archives do not work as well as a clone with submodules) |
| FEISTY | https://github.com/Kenhasteandersen/FEISTY/tree/fabm-feisty |
| ERSEM | https://github.com/pmlmodelling/ersem |
| FABM-PISCES | https://github.com/BoldingBruggeman/fabm-pisces (`fabm_pisces.yaml` names the revision it was made with) |
| TMM data | `MITgcm_2.8deg`, from the table of matrices on https://sites.google.com/view/samarkhatiwala-research-tmm (the table is embedded in the page; download the archive and extract it) |

Create the conda environment (conda-forge):

```bash
conda create -n fabmos -c conda-forge python=3.12 gfortran cmake ninja git pygetm scipy h5py netcdf4 cftime pyyaml xarray setuptools_scm mpi4py msmpi
conda activate fabmos
```

`msmpi` provides `mpiexec` on Windows; on Linux/macOS use `mpich` or `openmpi` instead.

### 2. Fix the FABM version bundled with FABM-OS

FABM-OS always compiles FABM from its own submodule `fabmos/extern/fabm`; the `-DFABM_BASE` CMake option does **not** select another FABM for the Python build.

FABM commits from b5704dbc ("prepare for initial_value without default of 0", fabm-model/fabm PR 179, July 2026) onwards crash when a model is created with FEISTY (`OSError: exception: access violation writing 0x0000000000000000` right after `Initializing fish...`). The crash is in FABM, not FEISTY: a leftover line in `src/fabm_types.F90` appends the optional argument `state_index` without checking that it is present, which FABM's depth-mapping code (used by FEISTY's zooplankton and benthos couplings) triggers. Until FABM is fixed, do one of:

- **Pin FABM** to the last revision tested with FEISTY:
  ```bash
  cd fabmos/extern/fabm
  git checkout f85de87a
  ```
- **Or patch** the bundled FABM: in `src/fabm_types.F90`, inside the block `if (variable%initial_value /= default_missing_value) then ... end if`, delete these two lines (the `state_index` append is already done, with a presence check, just before that block):
  ```fortran
  ! Store a pointer to the variable that should hold the state variable index.
  call variable%state_indices%append(state_index)
  ```
  With this patch the three configurations create correctly on FABM 87fb1b93 (Sep 2026); full runs were done on the pinned revision only.

### 3. Configure the build

Copy `setup.cfg` from this folder to the **root of the FABM-OS clone** (next to `setup.py`) and fill in the placeholders:

- `<ABSOLUTE_PATH>`: where the build files go. Use an absolute path; a relative one is resolved from the directory you run `pip` from.
- `<ENV>`: the fabmos conda environment (`where gfortran`, `where gcc-ar` and `where gcc-ranlib` show the paths).
- `<ERSEM>`, `<FABM-PISCES>`, `<FEISTY>`: the source folders from step 1. `FABM_EXTRA_INSTITUTES` and the `FABM_*_BASE` options register the three models with FABM.

If CMake cannot find `ninja.exe`, add `-DCMAKE_MAKE_PROGRAM="<ENV>/Library/bin/ninja.exe"` to `cmake_opts`.

### 4. Build and install FABM-OS

Build with two parallel jobs: with more, the FEISTY sources can compile before the FABM modules they use (`Cannot open module file 'fabm_types.mod'`).

PowerShell:
```powershell
conda activate fabmos
cd <fabmos>
$env:CMAKE_BUILD_PARALLEL_LEVEL = '2'
pip install .
```

cmd:
```bat
conda activate fabmos
cd <fabmos>
set CMAKE_BUILD_PARALLEL_LEVEL=2
pip install .
```

If two jobs still fail, use `1`. Before rebuilding, delete the build folder (`build_temp`) and also `<fabmos>/build` and `<fabmos>/src/fabmos.egg-info`: leftover libraries there are packaged into the new install and can stop it from loading.

Building with Intel `ifx` instead of gfortran also works, but needs the Visual Studio environment (`vcvars64`), `LIB` pointing to the oneAPI `compiler/<version>/lib` folder, the full path to `ifx.exe` in `CMAKE_Fortran_COMPILER`, and, at run time, the oneAPI `bin` folder added with `os.add_dll_directory` (Python's ctypes ignores `PATH`). gfortran is the simpler route.

See also https://github.com/BoldingBruggeman/fabmos/wiki.

## Running a case

1. **Prepare the TMM directory**: extract `MITgcm_2.8deg`.
2. **Copy the case files**: copy its yaml file and run script to the TMM directory.
3. **Run**, e.g.:
   ```bash
   conda activate fabmos
   cd <PATH_TO>/MITgcm_2.8deg
   mpiexec -n <NCORES> python run_ersem_noben.py
   ```
   Replace `<NCORES>` with the number of processes (e.g. `-n 4`).

Settings in the run scripts:

| Setting | Value |
|---------|-------|
| Time period | 2000-01-01 to 2010-01-01 (change `stop` for a short test) |
| Time step | 12 hours (43200 s) |
| Calendar | 360-day |
| Output | `output_<case>.nc`, monthly means |

One model month took about 65 s on one process for each case, and 25 s for `ersem_noben` with `mpiexec -n 4`; the monthly-mean output is about 110 MB per month.

The output contains the FEISTY fish biomass per functional type (`fish_smpel_totB`, `fish_mesopel_totB`, `fish_lgpel_totB`, `fish_bathypel_totB`, `fish_dem_totB`) and all state variables, including `fish_benthos` for the internal benthos. Change the `out.request()` call in a run script to add variables.

## Differences from the GOTM configurations

- `fish`: `max_depth` is 6000 m (GOTM: 5000) so that the look-up tables reach the deepest cells of the global grid; FEISTY uses the table closest to the bottom depth of each water column. `depth` is only used for the initial set-up and is 100 m (GOTM: 41 m).
- Inputs that the TMM host does not provide are set in the run scripts:
  - ERSEM: atmospheric CO2 of 400 ppm (`mole_fraction_of_carbon_dioxide_in_air`, as in GOTM), `absorption_of_silt` = 0.07 m-1 (as in GOTM) and a constant `bottom_stress` of 0.01 Pa (GOTM computed it).
  - PISCES: `surface_air_pressure` = 101325 Pa, and `vertical_tracer_diffusivity` from the crude `estimate_diffusivity()` defined in the script. Atmospheric CO2 (350 ppm) comes from the `pco2` instance in `fabm_pisces.yaml`.
- The nutrient profiles read from file in GOTM are not used; the nutrients start from the uniform `initialization` values in the yaml files (placeholders: use observed fields or a long spin-up).
- PISCES: the GOTM case relaxed nitrate, phosphate and silicate in the bottom 15 m to observed profiles (nutrient inflow). There is no direct equivalent in TMM and none is added here.
