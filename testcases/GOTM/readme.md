# FEISTY-UVic BGC in GOTM Test Case

Example configuration files for running FEISTY coupled with UVic biogeochemistry model in GOTM through FABM.

## Files Description

| File | Description |
|------|-------------|
| `fabm.yaml` | FABM configuration file defining biogeochemical models including FEISTY fish model, UVic components, and their couplings |
| `gotm.yaml` | GOTM configuration file containing physical ocean model settings and simulation parameters |
| `setup.zip` | Input data files required for the GOTM simulation (temperature/salinity profiles, meteorological forcing, grid specification, etc.) |

## Prerequisites

- GOTM compiled with FABM support
- FEISTY and UVic models registered in FABM

## Installation

1. **Download required components**:
   - GOTM: https://github.com/gotm-model/code/releases (Source code including submodules)
   - FABM: https://github.com/fabm-model/fabm/releases (Source code including externally maintained biogeochemical models)
   - UVic BGC: https://github.com/BoldingBruggeman/fabm-uvic
   - FEISTY: https://github.com/Kenhasteandersen/FEISTY/tree/fabm-feisty

2. **Register models in FABM**:

   > **Note**: GOTM may contain an outdated version of FABM. Updating to the latest FABM is recommended.

   a. Copy all files from the latest FABM release to `<GOTM>/extern/fabm/` (replace existing files).

   b. Place the UVic and FEISTY folders in `<GOTM>/extern/fabm/src/models/`.

   c. Edit `<GOTM>/extern/fabm/src/CMakeLists.txt` and add the UVic and FEISTY folder names to the `DEFAULT_INSTITUTES` list. For example:
   ```cmake
   set(DEFAULT_INSTITUTES
      ...
      feisty
      uvic
   )
   ```

   d. In `<GOTM>/extern/fabm/src/models/uvic/src/nut_chem.F90`, replace:
   ```fortran
   call self%register_dependency(self%id_atco2, standard_variables%mole_fraction_of_carbon_dioxide_in_air)
   ```
   with:
   ```fortran
   call self%register_dependency(self%id_atco2, 'atco2', '-', 'atmospheric CO2 concentration')
   ```

3. **Compile GOTM**: Follow the instructions at https://github.com/fabm-model/fabm/wiki/GOTM#building-from-scratch

## Running the Test Case

1. **Extract input data**: Unzip `setup.zip` to your GOTM build directory.

2. **Copy configuration files**: Copy `fabm.yaml` and `gotm.yaml` to the build directory, replacing any existing files.

3. **Run the simulation**:
   ```bash
   ./gotm
   ```

   > **Note**: Some commands may require administrator privileges (`sudo` on Linux/macOS).

## Simulation Details

- **Location**: North Sea station (56.43°N, 6.55°E)
- **Depth**: 41 m with 68 vertical layers
- **Time period**: 2003-01-01 to 2023-12-31
- **Time step**: 1800 s (30 minutes)

## Model Components (fabm.yaml)

- **FEISTY**: Fish ecosystem model with size-structured populations, coupled to UVic zooplankton and nutrient cycling
- **UVic biogeochemistry**: Sediment, nutrients/chemistry, phytoplankton, diazotrophs, zooplankton, detritus, light, and solar radiation

## Output

Results are saved daily as mean values. See the `output` section in `gotm.yaml` for configuration.
