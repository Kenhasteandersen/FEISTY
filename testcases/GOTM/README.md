# FABM-FEISTY

FEISTY (Fisheries Size and Functional Type model) as a [FABM](https://fabm.net) biogeochemical module. It adds a size-structured fish community to any FABM host (GOTM, FABM-OS/TMM, ...) and couples it to the zooplankton, benthos, detritus and nutrient pools of the host's biogeochemical model.

Model type in `fabm.yaml`: `feisty/FEISTY_FABM`.

This folder holds the GOTM test cases: `setup.zip` contains the forcing for a North Sea station and three FABM-FEISTY configurations, one with ERSEM and its internal benthos, one with ERSEM's benthic fauna as prey, and one with PISCES. The same three configurations for the transport-matrix model (FABM-OS) are in [`../TMM`](../TMM/readme.md). Earlier example configurations (FEISTY with UVic and with FABM's NPZD example) used an older FEISTY interface and have been retired.

## 1. What it is

**State variables**

- Fish biomass per functional type and size class, `g m-2` (wet weight), depth-integrated, named `<fft_short_name>_size_<k>`. There are five functional types (default names `fft_1`..`fft_5`; the test cases call them `smpel`, `mesopel`, `lgpel`, `bathypel`, `dem`), with `nStages` size classes for the large types. Each type also has a diagnostic `<fft_short_name>_totB` (total biomass).
- Optional internal benthos `benthos` (`g m-2`), used when `internal_benthos: true`.

All fish and benthos biomass is added to the host's total carbon, nitrogen and phosphorus, so conservation checks include it (`gww_gc`, `qnc` and `qpc` define the conversion).

**What it eats**

- Zooplankton, taken from any number of host models. Each is classed small (`s`) or large (`l`) and has its own availability.
- Benthos, either an internal pool fed by the deposition flux of the host's detritus, or any number of external zoobenthos models, each with an availability.
- Fish eat each other according to the FEISTY size-based preference.

**What it returns**

Respiration, feces and carcasses, each as carbon, nitrogen and phosphorus (any subset). Each flux goes to a host model you choose. Because fish take up prey in a fixed stoichiometry, any element in excess of the fish composition is excreted. Benthos excretion goes to the same respiration targets. Deposition that the internal benthos does not use goes back to the host's sediment.

**Vertical structure**

FEISTY works on depth-integrated prey (`mmol m-2`, integrated by FABM over the water column) and on mean temperatures (top 100 m, 500-1500 m, bottom). The fish are distributed vertically by pre-computed profiles that depend on bottom depth, photic depth, shelf depth and fish size, and the fluxes they return are spread over the water column with the same profile. The module therefore needs a water column (cell thickness, depth, bottom depth) from the host.

## 2. Installation

You need a host model with FABM, e.g. [GOTM](https://github.com/gotm-model/code) (which carries FABM in `extern/fabm`), and the FEISTY branch of this repository. FABM is built with CMake; FEISTY is added as an extra "institute".

Configure GOTM (or FABM) with:

```bash
cmake <GOTM_source> -DGOTM_USE_FABM=ON -DFABM_BASE=<GOTM_source>/extern/fabm \
      -DFABM_EXTRA_INSTITUTES="feisty" -DFABM_FEISTY_BASE=<path to this repo>
cmake --build . --config Release
```

Other biogeochemical models are added to the same list, each with a `FABM_<INSTITUTE>_BASE` if it is not bundled with FABM. The test cases need ERSEM and PISCES, e.g. `-DFABM_EXTRA_INSTITUTES="ersem;pisces;feisty" -DFABM_ERSEM_BASE=<ersem> -DFABM_PISCES_BASE=<fabm-pisces>` ([ERSEM](https://github.com/pmlmodelling/ersem), [FABM-PISCES](https://github.com/BoldingBruggeman/fabm-pisces)). Alternatively append `feisty` to `-DFABM_INSTITUTES=...` (the default list plus `feisty`). `FABM_FEISTY_BASE` must point at the directory containing this repository's `CMakeLists.txt`.

In FABM versions that look for models in `src/models/<institute>`, you can instead copy the `feisty` folder into `<fabm>/src/models/` and add `feisty` to `DEFAULT_INSTITUTES` in `<fabm>/src/CMakeLists.txt`. Setting `FABM_FEISTY_BASE` is simpler and leaves FABM untouched.

**FABM version.** FABM commits from b5704dbc (fabm-model/fabm PR 179, July 2026) onwards crash while creating any model that uses FABM's depth mapping, FEISTY included (`access violation` right after `Initializing fish...`). The test cases were run with FABM f85de87a (May 2026). With a newer FABM, check out f85de87a in `extern/fabm` or apply the two-line fix described in [`../TMM/readme.md`](../TMM/readme.md#2-fix-the-fabm-version-bundled-with-fabm-os).

**Running the GOTM test cases.** Unzip `setup.zip` into a folder `setup` and run, from inside that folder:

```bash
gotm gotm_ersem_noben.yaml   # ERSEM, internal benthos fed by deposition
gotm gotm_ersem_ben.yaml     # ERSEM with its benthic fauna Y2-Y4 as prey
gotm gotm_pisces.yaml        # PISCES, internal benthos
```

| File in `setup.zip` | Contents |
|---|---|
| `gotm_<case>.yaml` | GOTM settings for each case; points to `fabm_<case>.yaml` (`fabm: yaml_file:`), writes `result_<case>.nc` (monthly means of all variables) |
| `fabm_<case>.yaml` | FABM configuration for each case; every line carries the registered long name, units and default, and the header explains the coupling |
| `*.dat` | Forcing: grid, temperature and salinity profiles, meteorology, precipitation, short-wave radiation, surface elevation, pressure gradients, and nutrient, oxygen and carbonate profiles |
| `gotm.yaml`, `fabm.yaml` | The original physics set-up of the forcing package (MSI-ERGOM biogeochemistry, no FEISTY) |

The station is in the central North Sea (56.43°N, 6.55°E), 41 m deep with 68 layers, run from 2003-01-01 to 2023-12-31 with a 300 s time step. ERSEM is not mass-conserving at the 1800 s time step of the original set-up: its detritus goes negative and the clipping by `repair_state` adds nutrients. The PISCES case relaxes nitrate, phosphate and silicate in the bottom 15 m to the observed profiles, which stands in for the nutrient inflow that PISCES's sediment losses would otherwise remove from this closed column. A full run takes 6-10 minutes; FEISTY writes its look-up tables to `matrix.dat` (78 MB) on the first run and reuses them afterwards. The fish need several years to build up from their initial biomass; for results near equilibrium, rerun with `restart: load: true` so that the run starts from the state in `restart.nc` at the end of the previous one.

## 3. Setting up the FEISTY section of `fabm.yaml`

An instance is an ordinary FABM instance. The examples below name it `fish`:

```yaml
instances:
  fish:
    model: feisty/FEISTY_FABM
    parameters: {...}
    coupling: {...}
```

Parameters without a default are required.

### 3.1 Physiology

| Parameter | Units | Default | Description |
|---|---|---|---|
| `h` | g^nn yr-1 | 20 | Max. consumption coefficient |
| `nn` | - | -0.25 | Metabolic exponent |
| `gamma` | m2 g^q yr-1 | 70 | Coef. for clearance rate |
| `q` | - | -0.2 | Clearance rate exponent |
| `p` | - | -0.175 | Metabolism exponent |
| `epsAssim` | - | 0.7 | Assimilation efficiency |
| `epsRepro` | - | 0.01 | Reproduction and recruitment efficiency |
| `bgmort` | yr-1 | 0.1 | Fish background mortality |
| `beta` | - | 400 | Beta for size-based predation preference |
| `sigma` | - | 1.3 | Sigma for size-based predation preference |
| `mMin` | g | 0.001 | Minimum fish mass (boundary of the grid) |
| `etamature` | - | 0.25 | Determines the size with 50% maturity |
| `visual` | - | 1.5 | Visual ability (1.5 visual, 1 non-visual predator) |
| `Fmax` | yr-1 | 0 | Maximum fishing mortality |
| `etaF` | - | 0.05 | Determines the size with 50% fishing selectivity |
| `nStages` | - | 9 | Number of size classes of the large functional types |
| `b_ini` | g m-2 | 1e-5 | Initial biomass of each fish size class |

### 3.2 Environment and look-up tables

FEISTY pre-computes vertical distributions and predator-prey preferences (`n_depth` bottom depths from 0 to `max_depth`) and picks the closest one to the local bottom depth.

| Parameter | Units | Default | Description |
|---|---|---|---|
| `depth` | m | 100 | Water column depth used to set up the model |
| `max_depth` | m | 5000 | Maximum depth of the model domain (look-up tables) |
| `n_depth` | - | 50 | Number of depth levels in the look-up tables |
| `photic_depth` | m | 100 | Photic zone depth (limited to the local bottom depth) |
| `shelfdepth` | m | 250 | Continental shelf depth (no vertical migration at or above it) |
| `Tp` | Celsius | 10 | Pelagic layer averaged temperature (used for set-up) |
| `Tb` | Celsius | 8 | Bottom temperature (500-1500 m) (used for set-up) |
| `Tm` | Celsius | `Tb` | Mid-water temperature (used for set-up) |
| `ssigma` | - | 10 | Width of initial vertical distribution |
| `tau` | - | 10 | Increase in width |

At run time the temperatures come from the host's temperature field; `Tp`, `Tb` and `Tm` only affect the tables.

### 3.3 Stoichiometry

| Parameter | Units | Default | Description |
|---|---|---|---|
| `gww_gc` | - | 9 | Gram wet weight per gram carbon |
| `qnc` | N/C | 16/106 | Nitrogen to carbon ratio of fish (mol/mol) |
| `qpc` | P/C | 1/106 | Phosphorus to carbon ratio of fish (mol/mol) |

Prey stoichiometry is variable: prey are converted to wet weight by their limiting element, and fish excrete the excess of the other elements through the respiration targets.

### 3.4 Functional-type names

`fft_short_name_<i>` for `i` = 1..5 (default `fft_<i>`), unit-less. They name the state variables and the flux diagnostics, e.g. `smpel_size_1`.

### 3.5 Zooplankton prey

| Parameter | Default | Description |
|---|---|---|
| `n_zoo` | (required) | Number of zooplankton prey |
| `zoo_size<i>` | (required if `n_zoo` > 1) | `s` (small) or `l` (large) for prey `i` |
| `zoo_<i>_availability` | 1 | Availability of prey `i`, between 0 and 1 (only used if `n_zoo` > 1) |

Couplings (each takes the name of a host model; FEISTY reads its total C, N and P through the standard variables, so you point to the model once):

- `n_zoo` > 1: `small_zooplankton_target<j>` for the j-th prey with `zoo_size` `s`, and `large_zooplankton_target<j>` for the j-th prey with `l`, numbered in the order of the `zoo_size<i>` entries.
- `n_zoo` = 1: `zooplankton_target`. The single prey is split evenly between the small and the large pool. `zoo_size1` and the availability are not read.

### 3.6 Benthos

Set `internal_benthos` (required, `true`/`false`).

**Internal benthos** (`internal_benthos: true`), fed by the host's deposition of detritus:

| Parameter | Units | Default | Description |
|---|---|---|---|
| `benk` | g m-2 | (required) | Benthos carrying capacity |
| `ben_eff` | - | (required) | Detrital fraction going to the benthos |
| `bent_ini` | g m-2 | 5 | Initial benthos biomass |
| `det_c_scale`, `det_n_scale`, `det_p_scale`, `det_s_scale` | mmol (deposited unit)-1 | 1 | Convert deposited C, N, P, Si to mmol |
| `bottom_detritus_elements` | | (required) | Elements of the host's bottom detritus, e.g. `cnps` |

- FEISTY exposes four bottom deposition pools `fish/det_pool_c`, `fish/det_pool_n`, `fish/det_pool_p`, `fish/det_pool_s`. Couple the host's deposition into them (the host model that sinks the matter, e.g. ERSEM `deposition_target1c: fish/det_pool_c`, PISCES `bc: fish/det_pool_c`). FEISTY reads the deposition flux from them (`det_flux_c`, `det_flux_n`, `det_flux_p`, `det_flux_s`, which by default point at `./det_pool_<x>_sms_tot`).
- The part of the deposition the benthos does not use is returned to the host's sediment. For each element listed in `bottom_detritus_elements`, give `det_bot_<x>_target_real` (`x` = `c`, `n`, `p`, `s`), the host model that receives it.
- Benthic respiration goes to the same targets as the fish (`respiration_c_target`, ...).

**External benthos** (`internal_benthos: false`):

| Parameter | Default | Description |
|---|---|---|
| `n_ben` | (required) | Number of zoobenthos prey |
| `ben_<i>_availability` | (required) | Availability of prey `i`, between 0 and 1 |

Couplings: `zoobenthos_target<i>` for each prey. FEISTY reads their total C, N and P.

### 3.7 Return fluxes

| Parameter | Description |
|---|---|
| `respiration_elements` | Elements returned as respiration, e.g. `cnp`. Required. |
| `feces_elements` | Elements returned as feces. Required. |
| `carcass_elements` | Elements returned as carcasses. Required. |

For each element in each of these there is a target coupling, which takes a host model:

| | c | n | p |
|---|---|---|---|
| respiration | `respiration_c_target` | `respiration_n_target` | `respiration_p_target` |
| feces | `feces_c_target` | `feces_n_target` | `feces_p_target` |
| carcasses | `carcasses_c_target` | `carcasses_n_target` | `carcasses_p_target` |

Respiration is normally returned to the inorganic pools (DIC, ammonium, phosphate), feces and carcasses to a POM pool. The target model is resolved through the standard variables `total_carbon`, `total_nitrogen` and `total_phosphorus`. Each flux is added to the target model's pool in the water column using the fish vertical distribution. The flux diagnostics are named `<flux>_fish_<fft_short_name>_size<k>` (`mmol m-2 s-1`).

### 3.8 Units, and when FABM converts them

- **Standard-variable couplings are converted.** Zooplankton and zoobenthos (`total_carbon` etc.), return-flux targets and `det_bot_*_target_real` are all requested through standard variables. FABM converts from the host's units (e.g. ERSEM mg C, mmol N) to FEISTY's mmol, so nothing has to be scaled.
- **Direct variable couplings are not converted.** Coupling a host's deposition flux straight into `fish/det_pool_c` is a plain variable link: the number is passed unchanged. Set `det_c_scale`, `det_n_scale`, `det_p_scale` and `det_s_scale` to convert the deposited unit to mmol. The same holds for `det_flux_*` and `det_bot_<x>` if you set them by hand.

### 3.9 Practical pitfalls

- **Totals that include calcite.** In PISCES the total carbon of `gom` and of `sed` includes calcite. FEISTY therefore does not use `gom` as a target for feces or carcasses (`pom` is used), and the sediment's carbon is coupled with `det_bot_c: sed/c` (the organic pool) instead of the total.
- **Hosts with fixed C:N:P.** PISCES carries nitrogen and phosphorus implicitly in carbon units (122:16:1). Set `qnc = 16/122` and `qpc = 1/122`, use carbon-only `feces_elements` and `carcass_elements` (`c`), and use `bottom_detritus_elements: c`. The nitrogen and phosphorus deposition are derived from the carbon flux: couple `det_flux_n` and `det_flux_p` to `fish/det_pool_c_sms_tot` and set `det_n_scale` = 16/122 and `det_p_scale` = 1/122.
- **`require_initialization`.** FEISTY's state variables have initial values (`b_ini`, `bent_ini`), but the host models used in the test cases leave some of their own state variables without an explicit initial value. These cases set `require_initialization: false` (a top-level FABM setting) so those start from their default. Remove it if you initialize everything explicitly.
- **Check conservation.** `check_conservation: true` at the top level is useful: fish, the benthos and all return fluxes enter the total C, N and P.
- **Names are case-sensitive**, and the numbering in `zoo_size<i>`/`..._target<j>` starts at 1.

### 3.10 Examples

All three are in `setup.zip`. Every line there is annotated with units and defaults, and the file headers explain the choices.

**ERSEM, internal benthos** (`fabm_ersem_noben.yaml`). Z5 (micro) is small and Z4 (meso) is large prey. `det_c_scale` is 1/12.011 because ERSEM deposits carbon in mg C. The deposition not used by the benthos is returned to Q6.

```yaml
   fish:
      model: feisty/FEISTY_FABM
      parameters:
         qnc: 0.1509433962264151        # 16/106
         qpc: 0.009433962264150943      # 1/106
         feces_elements: cnp
         carcass_elements: cnp
         respiration_elements: cnp
         n_zoo: 2
         zoo_size1: s
         zoo_size2: l
         internal_benthos: true
         benk: 80.0
         ben_eff: 0.1
         det_c_scale: 0.0832570         # ERSEM carbon is in mg C: 1/12.011
         bottom_detritus_elements: cnps
      coupling:
         small_zooplankton_target1: Z5
         large_zooplankton_target1: Z4
         det_bot_c_target_real: Q6      # likewise _n_, _p_, _s_
         respiration_c_target: O3
         respiration_n_target: N4
         respiration_p_target: N1
         feces_c_target: R8             # likewise _n_, _p_ and carcasses_*_target
```

The ERSEM particulate matter is coupled into the pools in its own section, e.g. for R4, R6, R8 and P1-P4:

```yaml
         deposition_target1c: fish/det_pool_c
         deposition_target1n: fish/det_pool_n
         deposition_target1p: fish/det_pool_p
```

**ERSEM, benthic fauna as prey** (`fabm_ersem_ben.yaml`). Deposit feeders Y2, suspension feeders Y3 and meiofauna Y4 are the benthic prey.

```yaml
      parameters:
         internal_benthos: false
         n_ben: 3
         ben_1_availability: 1.0
         ben_2_availability: 1.0
         ben_3_availability: 1.0
         # zooplankton, stoichiometry and return fluxes as above
      coupling:
         zoobenthos_target1: Y2
         zoobenthos_target2: Y3
         zoobenthos_target3: Y4
```

**PISCES** (`fabm_pisces.yaml`). Matching 122:16:1 ratios, carbon-only feces and carcasses to `pom`, N and P deposition derived from the carbon flux, and the sediment's organic carbon `sed/c` as the real bottom pool.

```yaml
   fish:
      model: feisty/FEISTY_FABM
      parameters:
         qnc: 0.13114754098360656       # 16/122
         qpc: 0.00819672131147541       # 1/122
         feces_elements: c
         carcass_elements: c
         respiration_elements: cnp
         n_zoo: 2
         zoo_size1: s
         zoo_size2: l
         internal_benthos: true
         benk: 80.0
         ben_eff: 0.1
         det_c_scale: 1.0
         det_n_scale: 0.13114754098360656
         det_p_scale: 0.00819672131147541
         bottom_detritus_elements: c
      coupling:
         small_zooplankton_target1: zoo
         large_zooplankton_target1: mes
         det_bot_c_target_real: sed
         det_flux_n: fish/det_pool_c_sms_tot
         det_flux_p: fish/det_pool_c_sms_tot
         det_bot_c: sed/c
         respiration_c_target: carbonate
         respiration_n_target: nh4
         respiration_p_target: po4
         feces_c_target: pom
         carcasses_c_target: pom
```

In the PISCES `pom` and `gom` sections, `bc: fish/det_pool_c` sends the deposited carbon (already in mmol C m-2) to the benthos.
