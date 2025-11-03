# NPT Simulation Analysis for the intial 100 angstrom simulation box
from modules import analysis

# 0.1 GPa
analysis.npt_analysis(
    energy_xvg="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/npt_nap_pvp_pot01.xvg",
    temperature_xvg="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/npt_nap_pvp_temp01.xvg",
    pressure_xvg="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/npt_nap_pvp_pressure01.xvg",
    topology_file="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/wet_nap_pvp_npt01.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/wet_nap_pvp_npt01_trim.xtc",
    output_dir="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang",
)

# 0.2 Gpa
analysis.npt_analysis(
    energy_xvg="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/npt_nap_pvp_pot02.xvg",
    temperature_xvg="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/npt_nap_pvp_temp02.xvg",
    pressure_xvg="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/npt_nap_pvp_pressure02.xvg",
    topology_file="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/wet_nap_pvp_npt02.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/wet_nap_pvp_npt02_trim.xtc",
    output_dir="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang",
)
# 0.3 Gpa
analysis.npt_analysis(
    energy_xvg="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/npt_nap_pvp_pot03.xvg",
    temperature_xvg="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/npt_nap_pvp_temp03.xvg",
    pressure_xvg="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/npt_nap_pvp_pressure03.xvg",
    topology_file="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/wet_nap_pvp_npt03.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/wet_nap_pvp_npt03_trim.xtc",
    output_dir="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang",
)
# 0.4 Gpa
analysis.npt_analysis(
    energy_xvg="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/npt_nap_pvp_pot04.xvg",
    temperature_xvg="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/npt_nap_pvp_temp04.xvg",
    pressure_xvg="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/npt_nap_pvp_pressure04.xvg",
    topology_file="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/wet_nap_pvp_npt04.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/wet_nap_pvp_npt04_trim.xtc",
    output_dir="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang",
)
# 0.5 Gpa
analysis.npt_analysis(
    energy_xvg="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/npt_nap_pvp_pot05.xvg",
    temperature_xvg="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/npt_nap_pvp_temp05.xvg",
    pressure_xvg="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/npt_nap_pvp_pressure05.xvg",
    topology_file="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/wet_nap_pvp_npt05.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang/wet_nap_pvp_npt05_trim.xtc",
    output_dir="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang",
)
