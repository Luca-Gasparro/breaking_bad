import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utility import api_simulation_extractor, dummy_universe, polymer_atom_extraction


def api_api_rdf(
    topology_file, trajectory_file, api_residue_name, output_filename, start_time
):
    """Computes and plots the RDF for API-API. Uses the centres of mass of the API
    as reference particles."""

    api_simulation_extractor(
        topology_file=topology_file,
        trajectroy_file=trajectory_file,
        api_residue_name=api_residue_name,
        output_filename=output_filename,
        start_time=start_time,
    )
    api_universe = dummy_universe(api_simulation_file=output_filename)

    api_rdf = InterRDF(
        api_universe.atoms,
        api_universe.atoms,
        nbins=75,
        range=(0, 30),
        exclusion_block=(1, 1),
    )
    api_rdf.run()

    return api_rdf


def polymer_polymer_rdf(
    polymer_topology_file, polymer_trajectory_file, atom_name, start_time
):
    """Computes radial distribution function for polymer-polymer using a selected reference atom.
    Frames after the start time are used in the calculation."""

    polymer_traj = mda.Universe(polymer_topology_file, polymer_trajectory_file)
    timestep = polymer_traj.trajectory.dt
    start_frame = int(start_time / timestep)

    atom_selection = polymer_traj.select_atoms(f"name {atom_name}")

    polymer_rdf = InterRDF(
        atom_selection, atom_selection, nbins=75, range=(0, 40), exclusion_block=(1, 1)
    )
    polymer_rdf.run(start=start_frame)

    plt.figure(figsize=(8, 6))
    plt.plot(polymer_rdf.results.bins, polymer_rdf.results.rdf)
    plt.title("poly-poly rdf")
    plt.ylabel("rdf")
    plt.xlabel("angstroms")
    plt.savefig("poly-poly_rdf_test.png", dpi=300)

    return


polymer_polymer_rdf("dry_nvt_polymer.tpr", "dry_nvt_trim_polymers.xtc", "N1", 2000)
