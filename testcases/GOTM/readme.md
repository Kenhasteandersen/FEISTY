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
   - FABM with GOTM: https://github.com/fabm-model/fabm/releases
   - UVic BGC: https://github.com/BoldingBruggeman/fabm-uvic
   - FEISTY: https://github.com/Kenhasteandersen/FEISTY

2. **Register models in FABM**:

   Place UVic and FEISTY folders under:
   ```
   .../extern/fabm/src/models/
   ```

   Then edit `.../extern/fabm/src/CMakeLists.txt` and add the UVic and FEISTY folder names from the step above. For example:
   ```cmake
   set(DEFAULT_INSTITUTES
      ...
      feisty
      uvic
   )
   ```

3. **Compile GOTM**: Follow the instructions at https://github.com/fabm-model/fabm/wiki/Building-and-installing

## Running the Test Case

1. **Extract input data**: Unzip `setup.zip` to your GOTM build directory.

2. **Copy configuration files**: Copy `fabm.yaml` and `gotm.yaml` to the build directory, replacing any existing files.

3. **Run the simulation**:
   ```
   ./gotm
   ```

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
