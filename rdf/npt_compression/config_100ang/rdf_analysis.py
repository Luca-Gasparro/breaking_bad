# RDFS for initial 100 angstrom simulation box

from modules import utility
from modules import rdf

# 0.1 GPa
utility.api_simulation_extractor(
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt01.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt01_trim.xtc",
    api_residue_name="NAP",
    api_information_filename="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt01.npz",
)

api_bin_01, api_rdf_01 = rdf.api_api_rdf(
    api_information_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt01.npz",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

poly_bin_01, poly_rdf_01 = rdf.polymer_polymer_rdf(
    polymer_topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt01.tpr",
    polymer_trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt01_trim.xtc",
    reference_particle_name="N1",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

cross_bin_01, cross_rdf_01 = rdf.polymer_api_rdf(
    api_information_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt01.npz",
    polymer_topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt01.tpr",
    polymer_trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt01_trim.xtc",
    reference_particle_name="N1",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

rdf.plot_rdfs(
    api_bin=api_bin_01,
    api_rdf=api_rdf_01,
    poly_bin=poly_bin_01,
    poly_rdf=poly_rdf_01,
    cross_bin=cross_bin_01,
    cross_rdf=cross_rdf_01,
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt01.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt01_trim.xtc",
    npt_mdp_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/params/01npt.mdp",
    output_dir="/storage/chem/phuqdw/breaking_bad/rdf/npt_compression/config_100ang",
)

# 0.2 Gpa
utility.api_simulation_extractor(
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt02.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt02_trim.xtc",
    api_residue_name="NAP",
    api_information_filename="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt02.npz",
)

api_bin_02, api_rdf_02 = rdf.api_api_rdf(
    api_information_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt02.npz",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

poly_bin_02, poly_rdf_02 = rdf.polymer_polymer_rdf(
    polymer_topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt02.tpr",
    polymer_trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt02_trim.xtc",
    reference_particle_name="N1",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

cross_bin_02, cross_rdf_02 = rdf.polymer_api_rdf(
    api_information_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt02.npz",
    polymer_topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt02.tpr",
    polymer_trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt02_trim.xtc",
    reference_particle_name="N1",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

rdf.plot_rdfs(
    api_bin=api_bin_02,
    api_rdf=api_rdf_02,
    poly_bin=poly_bin_02,
    poly_rdf=poly_rdf_02,
    cross_bin=cross_bin_02,
    cross_rdf=cross_rdf_02,
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt02.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt02_trim.xtc",
    npt_mdp_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/params/02npt.mdp",
    output_dir="/storage/chem/phuqdw/breaking_bad/rdf/npt_compression/config_100ang",
)

# 0.3 Gpa
utility.api_simulation_extractor(
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt03.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt03_trim.xtc",
    api_residue_name="NAP",
    api_information_filename="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt03.npz",
)

api_bin_03, api_rdf_03 = rdf.api_api_rdf(
    api_information_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt03.npz",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

poly_bin_03, poly_rdf_03 = rdf.polymer_polymer_rdf(
    polymer_topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt03.tpr",
    polymer_trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt03_trim.xtc",
    reference_particle_name="N1",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

cross_bin_03, cross_rdf_03 = rdf.polymer_api_rdf(
    api_information_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt03.npz",
    polymer_topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt03.tpr",
    polymer_trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt03_trim.xtc",
    reference_particle_name="N1",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

rdf.plot_rdfs(
    api_bin=api_bin_03,
    api_rdf=api_rdf_03,
    poly_bin=poly_bin_03,
    poly_rdf=poly_rdf_03,
    cross_bin=cross_bin_03,
    cross_rdf=cross_rdf_03,
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt03.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt03_trim.xtc",
    npt_mdp_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/params/03npt.mdp",
    output_dir="/storage/chem/phuqdw/breaking_bad/rdf/npt_compression/config_100ang",
)

# 0.4 Gpa
utility.api_simulation_extractor(
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt04.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt04_trim.xtc",
    api_residue_name="NAP",
    api_information_filename="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt04.npz",
)

api_bin_04, api_rdf_04 = rdf.api_api_rdf(
    api_information_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt04.npz",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

poly_bin_04, poly_rdf_04 = rdf.polymer_polymer_rdf(
    polymer_topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt04.tpr",
    polymer_trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt04_trim.xtc",
    reference_particle_name="N1",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

cross_bin_04, cross_rdf_04 = rdf.polymer_api_rdf(
    api_information_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt04.npz",
    polymer_topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt04.tpr",
    polymer_trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt04_trim.xtc",
    reference_particle_name="N1",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

rdf.plot_rdfs(
    api_bin=api_bin_04,
    api_rdf=api_rdf_04,
    poly_bin=poly_bin_04,
    poly_rdf=poly_rdf_04,
    cross_bin=cross_bin_04,
    cross_rdf=cross_rdf_04,
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt04.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt04_trim.xtc",
    npt_mdp_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/params/04npt.mdp",
    output_dir="/storage/chem/phuqdw/breaking_bad/rdf/npt_compression/config_100ang",
)

# 0.5 Gpa
utility.api_simulation_extractor(
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt05.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt05_trim.xtc",
    api_residue_name="NAP",
    api_information_filename="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt05.npz",
)

api_bin_05, api_rdf_05 = rdf.api_api_rdf(
    api_information_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt05.npz",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

poly_bin_05, poly_rdf_05 = rdf.polymer_polymer_rdf(
    polymer_topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt05.tpr",
    polymer_trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt05_trim.xtc",
    reference_particle_name="N1",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

cross_bin_05, cross_rdf_05 = rdf.polymer_api_rdf(
    api_information_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/api_coms_npt05.npz",
    polymer_topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt05.tpr",
    polymer_trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/poly_npt05_trim.xtc",
    reference_particle_name="N1",
    bin_width=0.5,
    instantaneous=True,
    start_time=5000,
)

rdf.plot_rdfs(
    api_bin=api_bin_05,
    api_rdf=api_rdf_05,
    poly_bin=poly_bin_05,
    poly_rdf=poly_rdf_05,
    cross_bin=cross_bin_05,
    cross_rdf=cross_rdf_05,
    topology_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt05.tpr",
    trajectory_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/config_100ang/wet_nap_pvp_npt05_trim.xtc",
    npt_mdp_file="/storage/chem/phuqdw/breaking_bad/gromacs/npt_compression/params/05npt.mdp",
    output_dir="/storage/chem/phuqdw/breaking_bad/rdf/npt_compression/config_100ang",
)
