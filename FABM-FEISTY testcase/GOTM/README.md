# GOTM test case: FEISTY at a North Sea station

A 1-D GOTM water column in the central North Sea (56.43°N, 6.55°E), 41 m deep with 68 layers, run from 2003-01-01 to 2023-12-31 with a 300 s time step. FEISTY runs on top of one of three biogeochemical set-ups:

| Case | GOTM file | FABM file | Biogeochemistry | Benthic prey of FEISTY |
|---|---|---|---|---|
| `ersem_noben` | `gotm_ersem_noben.yaml` | `fabm_ersem_noben.yaml` | ERSEM 15.06 with simple benthic pools | FEISTY's internal benthos, fed by ERSEM's detritus deposition |
| `ersem_ben` | `gotm_ersem_ben.yaml` | `fabm_ersem_ben.yaml` | ERSEM 15.06 with the full benthic model | ERSEM's benthic fauna (deposit feeders Y2, suspension feeders Y3, meiofauna Y4) |
| `pisces` | `gotm_pisces.yaml` | `fabm_pisces.yaml` | PISCES | FEISTY's internal benthos, fed by PISCES's detritus deposition |

How FEISTY is set up inside each `fabm_<case>.yaml`, and all its parameters, are described in [`../README.md`](../README.md).

## Files

| File | Contents |
|---|---|
| `gotm_<case>.yaml` | GOTM settings: physics, forcing, the FABM file to use (`fabm: yaml_file:`), the biogeochemical initial conditions and relaxation (`fabm: input:`), and the output file `result_<case>.nc` (monthly means of all variables) |
| `fabm_<case>.yaml` | FABM configuration: the biogeochemical model plus the FEISTY instance `fish`. Every line carries the long name, units and default, and the header explains the coupling |
| `*.dat` | Forcing shared by all cases: grid, temperature and salinity profiles, meteorology, precipitation, short-wave radiation, surface elevation, pressure gradients, and nutrient, oxygen and carbonate profiles |

## Running

GOTM writes its output (and FEISTY its look-up table `matrix.dat`, 78 MB) to the working directory, so run from a copy of this folder:

```bash
gotm gotm_ersem_noben.yaml
gotm gotm_ersem_ben.yaml
gotm gotm_pisces.yaml
```

A full run takes 6-10 minutes. `matrix.dat` is written on the first run and reused afterwards. It is checked only by its array sizes, so delete it after changing FEISTY parameters that affect the vertical distribution or the predator-prey preferences.

The fish need several years to build up from their initial biomass. For results closer to equilibrium, rerun with `restart: load: true` so that the run starts from `restart.nc`, the state at the end of the previous run.

## Switching between the three set-ups

Each set-up is a pair of files, `gotm_<case>.yaml` and `fabm_<case>.yaml`. To switch, pass a different GOTM file on the command line. The forcing files are shared, so nothing else has to change.

The three GOTM files differ only in what the biogeochemical model needs:

- `fabm: yaml_file:` names the FABM file of the case.
- `fabm: input:` sets the initial nutrient, oxygen and carbonate profiles. The variables are named after the biogeochemical model (`N3/n`, `N1/p`, `O3/c`, ... for ERSEM; `no3/c`, `po4/c`, `carbonate/DIC`, ... for PISCES, with scale factors for its mol/L units). ERSEM also needs `absorption_of_silt`. PISCES relaxes nitrate, phosphate and silicate in the bottom 15 m to the observed profiles, which stands in for the lateral nutrient supply that this closed column lacks.
- `output:` names the output file `result_<case>.nc`, and `title:` describes the case.

To run your own variant, copy one pair, edit the FABM file, and point the copied GOTM file to it with `fabm: yaml_file:`. Change the output name too, so that runs do not overwrite each other. Keep the 300 s time step for ERSEM: at 1800 s its detritus goes negative, and the clipping by `repair_state` then adds nitrogen, phosphorus and silicate.

To change fishing, edit `Fsource_<i>`, `Fmax_<i>` and `etaF_<i>` in the FEISTY section of the FABM file. With `Fsource_<i>: 1`, the fishing mortality comes from an input `fish/F_<name>` listed under `fabm: input:` in the GOTM file (see section 3.4 of [`../README.md`](../README.md)).
