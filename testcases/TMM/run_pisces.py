# TMM-FEISTY test case: PISCES with FEISTY's internal benthos.
# FABM configuration: fabm_pisces.yaml (adapted from testcases/GOTM/setup/fabm_pisces.yaml).
# Run from the TMM directory (e.g. MITgcm_2.8deg) after copying fabm_pisces.yaml and this script there:
#   mpiexec -n <NCORES> python run_pisces.py
import cftime
import os

import numpy as np
import fabmos.transport.tmm
import fabmos

calendar = "360_day"  # any valid calendar recognized by cftime, see https://cfconventions.org/cf-conventions/cf-conventions.html#calendar

domain = fabmos.transport.tmm.create_domain(".")

sim = fabmos.transport.tmm.Simulator(domain, calendar=calendar, fabm="fabm_pisces.yaml")

def estimate_diffusivity():
    # Crude estimate of turbulent diffusivity from temperature-only MLD
    # PISCES and iHAMOCC use this to determine turbocline depth
    temp = fabmos.transport.tmm.get_mat_array(
        "GCM/Theta_gcm.mat",
        "Tgcm",
        domain._grid_file,
        times=fabmos.transport.tmm.climatology_times(calendar),
    )
    Kzval = np.where(temp.values < temp.values[:, :1, :, :] - 0.5, 1e-5, 1e-3)
    return fabmos.transport.tmm._wrap_ongrid_array(
        Kzval, domain._grid_file, times=fabmos.transport.tmm.climatology_times(calendar)
    )

# Inputs that GOTM supplied through fabm: input: in testcases/GOTM/setup/gotm_pisces.yaml.
# Atmospheric CO2 is set in fabm_pisces.yaml (instance pco2), not here: PISCES has no mole_fraction_of_carbon_dioxide_in_air dependency.
sim.fabm.get_dependency("surface_air_pressure").set(101325.0)  # surface air pressure [Pa]: the TMM host does not supply it

# PISCES (turbocline depth) needs a turbulent diffusivity, which GOTM computed itself. Crude estimate from temperature, see above.
sim.fabm.get_dependency("vertical_tracer_diffusivity").set(
    estimate_diffusivity(), on_grid=fabmos.input.OnGrid.ALL, climatology=True
)

# NOTE: gotm_pisces.yaml relaxed nitrate, phosphate and silicate in the bottom 15 m to observed profiles (nutrient inflow).
# TMM has no direct equivalent, so nothing is done here; nutrients are only transported and cycled by the model.
# If a nutrient source/relaxation is needed, it has to be designed for the global grid (not invented here).

# Other dependencies (temperature, salinity, light, wind, ...) are supplied by the TMM host; add more here if FABM reports an unset one.

out = sim.output_manager.add_netcdf_file(
    "output_pisces.nc", interval=1, interval_units=fabmos.TimeUnit.MONTHS
)
# FEISTY fish biomass diagnostics (g m-2). FEISTY's internal benthos (fish_benthos) is a state variable and is included below.
out.request(
    'fish_smpel_totB',
    'fish_mesopel_totB',
    'fish_lgpel_totB',
    'fish_bathypel_totB',
    'fish_dem_totB',
    *sim.fabm.state_variables,
    time_average=True,
)
start = cftime.datetime(2000, 1, 1, calendar=calendar)
stop = cftime.datetime(2010, 1, 1, calendar=calendar)
sim.start(start, timestep=12*3600)
while sim.time < stop:
    sim.advance()
sim.finish()
