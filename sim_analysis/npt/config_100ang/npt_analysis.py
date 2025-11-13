# NPT Simulation Analysis for the intial 100 angstrom simulation box
from modules import analysis

# 0.1 GPa
analysis.npt_analysis(
    energy_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/npt_nap_pvp_pot01.xvg",
    temperature_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/npt_nap_pvp_temp01.xvg",
    pressure_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/npt_nap_pvp_pressure01.xvg",
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt01.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt01_trim.xtc",
    npt_mdp_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/params/01npt.mdp",
    output_dir="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang",
)

# 0.2 GPa
analysis.npt_analysis(
    energy_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/npt_nap_pvp_pot02.xvg",
    temperature_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/npt_nap_pvp_temp02.xvg",
    pressure_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/npt_nap_pvp_pressure02.xvg",
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt02.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt02_trim.xtc",
    npt_mdp_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/params/02npt.mdp",
    output_dir="/storage/chem/phuqdw/breaking_bad/sim_analysis/npt/config_100ang",
)
