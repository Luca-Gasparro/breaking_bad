# Analysing energy minimisiation for each naproxen PVP configuration
import analysis

# 100 angstrom configuration...
analysis.energy_minimisation_pe(
    polymer_xvg="em_pvp_relax_100ang.xvg",
    api_polymer_xvg="em_nap_pvp_relax_100ang.xvg",
    wet_api_polymer_xvg="em_wet_nap_pvp_relax_100ang.xvg",
    box_size=100,
)

# 110 angstrom configuration
analysis.energy_minimisation_pe(
    polymer_xvg="em_pvp_relax_110ang.xvg",
    api_polymer_xvg="em_nap_pvp_relax_110ang.xvg",
    wet_api_polymer_xvg="em_wet_nap_pvp_relax_110ang.xvg",
    box_size=110,
)

# 130 angstrom configuration
analysis.energy_minimisation_pe(
    polymer_xvg="em_pvp_relax_130ang.xvg",
    api_polymer_xvg="em_nap_pvp_relax_130ang.xvg",
    wet_api_polymer_xvg="em_wet_nap_pvp_relax_130ang.xvg",
    box_size=130,
)

# 140 angstrom configuration
analysis.energy_minimisation_pe(
    polymer_xvg="em_pvp_relax_140ang.xvg",
    api_polymer_xvg="em_nap_pvp_relax_140ang.xvg",
    wet_api_polymer_xvg="em_wet_nap_pvp_relax_140ang.xvg",
    box_size=140,
)

# 160 angstrom configuration
analysis.energy_minimisation_pe(
    polymer_xvg="em_pvp_relax_160ang.xvg",
    api_polymer_xvg="em_nap_pvp_relax_160ang.xvg",
    wet_api_polymer_xvg="em_wet_nap_pvp_relax_160ang.xvg",
    box_size=160,
)
