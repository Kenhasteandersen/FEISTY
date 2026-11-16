# TMM-FEISTY test case: ERSEM with the full benthic model; ERSEM benthic fauna Y2, Y3, Y4 are FEISTY's benthic prey.
# FABM configuration: fabm_ersem_ben.yaml (adapted from testcases/GOTM/setup/fabm_ersem_ben.yaml).
# Run from the TMM directory (e.g. MITgcm_2.8deg) after copying fabm_ersem_ben.yaml and this script there:
#   mpiexec -n <NCORES> python run_ersem_ben.py
import cftime
import os

import numpy as np
import fabmos.transport.tmm
import fabmos

calendar = "360_day"  # any valid calendar recognized by cftime, see https://cfconventions.org/cf-conventions/cf-conventions.html#calendar

domain = fabmos.transport.tmm.create_domain(".")

sim = fabmos.transport.tmm.Simulator(domain, calendar=calendar, fabm="fabm_ersem_ben.yaml")

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

# Inputs that GOTM supplied through fabm: input: in testcases/GOTM/setup/gotm_ersem_ben.yaml.
sim.fabm.get_dependency("mole_fraction_of_carbon_dioxide_in_air").set(400.0)  # atmospheric pCO2 [ppm]
sim.fabm.get_dependency("absorption_of_silt").set(0.07)  # shortwave absorption by suspended silt [m-1], needed by ERSEM light_iop
sim.fabm.get_dependency("bottom_stress").set(0.01)  # bottom stress [Pa]: crude constant, the TMM host does not supply it (needed by ERSEM benthos)

# Other dependencies (temperature, salinity, light, wind, ...) are supplied by the TMM host; add more here if FABM reports an unset one.

out = sim.output_manager.add_netcdf_file(
    "output_ersem_ben.nc", interval=1, interval_units=fabmos.TimeUnit.MONTHS
)
# FEISTY fish biomass diagnostics (g m-2). The benthic prey are ERSEM's Y2, Y3, Y4 (included in the state variables below).
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
