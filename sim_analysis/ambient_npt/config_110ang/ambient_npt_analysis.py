# Ambient pressure simulation analysis
from modules import analysis

analysis.npt_analysis(
    energy_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_110ang/npt_ambient_pot.xvg",
    temperature_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_110ang/npt_ambient_temp.xvg",
    pressure_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_110ang/npt_ambient_pressure.xvg",
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_110ang/ambient_npt.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_110ang/ambient_npt_trim.xtc",
    npt_mdp_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_110ang/npt.mdp",
    output_dir="/storage/chem/phuqdw/breaking_bad/sim_analysis/ambient_npt/config_110ang",
)
