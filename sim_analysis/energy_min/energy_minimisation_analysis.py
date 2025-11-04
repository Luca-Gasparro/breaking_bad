# Analysing energy minimisiation for each naproxen PVP configuration

from modules import analysis


# 100 angstrom configuration...
analysis.energy_minimisation_pe(
    polymer_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/em_files/em_pvp_relax_100ang.xvg",
    api_polymer_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/em_files/em_nap_pvp_relax_100ang.xvg",
    wet_api_polymer_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/em_files/em_wet_nap_pvp_relax_100ang.xvg",
    box_size=100,
    output_dir="/storage/chem/phuqdw/breaking_bad/sim_analysis/energy_min",
)

# 110 angstrom configuration
analysis.energy_minimisation_pe(
    polymer_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/em_files/em_pvp_relax_110ang.xvg",
    api_polymer_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/em_files/em_nap_pvp_relax_110ang.xvg",
    wet_api_polymer_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/em_files/em_wet_nap_pvp_relax_110ang.xvg",
    box_size=110,
    output_dir="/storage/chem/phuqdw/breaking_bad/sim_analysis/energy_min",
)

# 130 angstrom configuration
analysis.energy_minimisation_pe(
    polymer_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/em_files/em_pvp_relax_130ang.xvg",
    api_polymer_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/em_files/em_nap_pvp_relax_130ang.xvg",
    wet_api_polymer_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/em_files/em_wet_nap_pvp_relax_130ang.xvg",
    box_size=130,
    output_dir="/storage/chem/phuqdw/breaking_bad/sim_analysis/energy_min",
)

# 140 angstrom configuration
analysis.energy_minimisation_pe(
    polymer_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/em_files/em_pvp_relax_140ang.xvg",
    api_polymer_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/em_files/em_nap_pvp_relax_140ang.xvg",
    wet_api_polymer_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/em_files/em_wet_nap_pvp_relax_140ang.xvg",
    box_size=140,
    output_dir="/storage/chem/phuqdw/breaking_bad/sim_analysis/energy_min",
)

# 160 angstrom configuration
analysis.energy_minimisation_pe(
    polymer_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/em_files/em_pvp_relax_160ang.xvg",
    api_polymer_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/em_files/em_nap_pvp_relax_160ang.xvg",
    wet_api_polymer_xvg="/storage/chem/phuqdw/breaking_bad/gromacs/em_files/em_wet_nap_pvp_relax_160ang.xvg",
    box_size=160,
    output_dir="/storage/chem/phuqdw/breaking_bad/sim_analysis/energy_min",
)
