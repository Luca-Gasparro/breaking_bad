# RDFS for initial 100 angstrom simulation box

from modules import utility
from modules import rdf_analysis

# 0.1 GPa
utility.api_simulation_extractor(
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt01.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt01_trim.xtc",
    api_residue_name="NAP",
    api_information_filename="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt01.npz",
)

api_bin_01, api_rdf_01 = rdf_analysis.api_api_rdf(
    api_information_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt01.npz",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

poly_bin_01, poly_rdf_01 = rdf_analysis.polymer_polymer_rdf(
    polymer_topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt01.tpr",
    polymer_trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt01_trim.xtc",
    reference_particle_name="N1",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

cross_bin_01, cross_rdf_01 = rdf_analysis.polymer_api_rdf(
    api_information_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt01.npz",
    polymer_topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt01.tpr",
    polymer_trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt01_trim.xtc",
    reference_particle_name="N1",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

rdf_analysis.plot_rdfs(
    api_bin=api_bin_01,
    api_rdf=api_rdf_01,
    poly_bin=poly_bin_01,
    poly_rdf=poly_rdf_01,
    cross_bin=cross_bin_01,
    cross_rdf=cross_rdf_01,
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt01.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt01_trim.xtc",
    output_dir="/storage/chem/phuqdw/breaking_bad/rdf/npt_compression/config_100ang",
)
