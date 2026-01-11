# FEISTY-UVic BGC in TMM Test Case

Configuration files and instructions for running FEISTY coupled with the UVic biogeochemistry model in TMM through FABM-OS.

## Files

| File | Description |
|------|-------------|
| `setup.cfg` | Build configuration for FABM-OS (paths, compiler settings, CMake options) |
| `fabm.yaml` | FABM configuration defining FEISTY and UVic components, and their couplings |
| `run.py` | Python script to configure and execute the simulation |


## Prerequisites

- Anaconda/Miniconda
- FABM-OS built from source with FEISTY and UVic models registered
- TMM files
- Ninja build system (may be required; install in fabmos environment: `conda activate fabmos && conda install ninja`)

## Installation

### 1. Download Required Components

| Component | Source |
|-----------|--------|
| FABM-OS | https://github.com/BoldingBruggeman/fabmos/releases |
| UVic BGC | https://github.com/BoldingBruggeman/fabm-uvic |
| FEISTY | https://github.com/Kenhasteandersen/FEISTY |
| TMM Data | https://sites.google.com/view/samarkhatiwala-research-tmm |

For TMM, download `MITgcm_2.8deg` for coarse resolution simulations.

### 2. Register Models in FABM

Place the UVic and FEISTY model folders under:
```
.../fabmos/extern/fabm/src/models/
```

Edit `.../fabmos/extern/fabm/src/CMakeLists.txt` to include the new models:
```cmake
set(DEFAULT_INSTITUTES
   ...
   feisty
   uvic
)
```

### 3. Build FABM-OS

Follow the build instructions at https://github.com/BoldingBruggeman/fabmos/wiki

Customize `setup.cfg` for your system:
```ini
[build_ext]
build_temp=../../build/fabmos
cmake_opts=-G Ninja -DFABM_BASE="<YOUR_PATH>/fabmos/extern/fabm" -DCMAKE_Fortran_COMPILER=gfortran
compiler=gfortran
```

## Running the Simulation

1. **Prepare TMM directory**: Download and extract `MITgcm_2.8deg`

2. **Copy configuration files**: Copy `fabm.yaml` and `run.py` to the TMM directory

3. **Execute the simulation**:
   ```bash
   conda activate fabmos
   cd <PATH_TO>/MITgcm_2.8deg
   mpiexec -n <NCORES> python run.py
   ```
   Replace `<NCORES>` with the number of CPU cores for parallel execution (e.g., `-n 10`).

## Simulation Configuration

Default settings in `run.py`:

| Parameter | Value |
|-----------|-------|
| Time period | 2000-01-01 to 2010-01-01 |
| Time step | 12 hours (43200 s) |
| Calendar | 360-day |
| Output interval | Monthly |
| Output format | NetCDF (`output.nc`) |

Modify `run.py` to adjust time period, timestep, or output variables.

## Output Variables

Results are saved as monthly time-averaged values to `output.nc`. Default requested variables:
- `fish_fft_1_totB` through `fish_fft_5_totB`: FEISTY fish biomass by functional type
- All FABM state variables

To customize output, modify the `out.request()` line in `run.py`.
