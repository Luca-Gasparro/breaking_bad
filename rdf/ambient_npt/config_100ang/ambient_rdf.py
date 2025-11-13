# RDF for initial 100 angstrom simulation box - ambient expansion

from modules import utility
from modules import rdf

utility.api_simulation_extractor(
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_100ang/ambient_npt.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_100ang/ambient_npt_trim.xtc",
    api_residue_name="NAP",
    api_information_filename="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_100ang/api_coms.npz",
)

api_bin, api_rdf = rdf.api_api_rdf(
    api_information_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_100ang/api_coms.npz",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

poly_bin, poly_rdf = rdf.polymer_polymer_rdf(
    polymer_topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_100ang/ambient_npt_poly.tpr",
    polymer_trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_100ang/ambient_npt_poly_trim.xtc",
    reference_particle_name="N1",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

cross_bin, cross_rdf = rdf.polymer_api_rdf(
    api_information_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_100ang/api_coms.npz",
    polymer_topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_100ang/ambient_npt_poly.tpr",
    polymer_trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_100ang/ambient_npt_poly_trim.xtc",
    reference_particle_name="N1",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

rdf.plot_rdfs(
    api_bin=api_bin,
    api_rdf=api_rdf,
    poly_bin=poly_bin,
    poly_rdf=poly_rdf,
    cross_bin=cross_bin,
    cross_rdf=cross_rdf,
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_100ang/ambient_npt.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_100ang/ambient_npt_trim.xtc",
    npt_mdp_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_ambient/config_100ang/npt.mdp",
    output_dir="/storage/chem/phuqdw/breaking_bad/rdf/ambient_npt/config_100ang",
)
