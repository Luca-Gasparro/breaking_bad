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
        trajectory_file=trajectory_file,
        api_residue_name=api_residue_name,
        output_filename=output_filename,
        start_time=start_time,
    )
    api_universe = dummy_universe(api_simulation_file=output_filename)

    # Calculate the API-API RDF
    api_rdf = InterRDF(
        api_universe.atoms,
        api_universe.atoms,
        nbins=75,
        range=(0, 30),
        exclusion_block=(1, 1),
    )
    api_rdf.run()

    # Return information needed to plot - the bins and the RDF
    return api_rdf.results.bins, api_rdf.results.rdf


def polymer_polymer_rdf(
    polymer_topology_file, polymer_trajectory_file, polymer_atom_name, start_time
):
    """Computes radial distribution function for polymer-polymer using a selected reference atom.
    Frames after the start time are used in the calculation."""

    # Extract the reference polymer atom and the timestep
    polymer_atom_selection, timestep = polymer_atom_extraction(
        polymer_topology_file=polymer_topology_file,
        polymer_trajectory_file=polymer_trajectory_file,
        atom_name=polymer_atom_name,
    )
    start_frame = int(start_time / timestep)

    # Calculate the polymer-polymer RDF
    polymer_rdf = InterRDF(
        polymer_atom_selection,
        polymer_atom_selection,
        nbins=75,
        range=(0, 40),
        exclusion_block=(1, 1),
    )
    polymer_rdf.run(start=start_frame)

    # Return information needed to plot - the bins and the RDF
    return polymer_rdf.results.bins, polymer_polymer_rdf.rdf


def api_polymer_rdf(
    topology_file,
    trajectory_file,
    api_residue_name,
    output_filename,
    polymer_topology_file,
    polymer_trajectory_file,
    polymer_atom_name,
    start_time,
):
    """Computes the radial distribution function between the API and
    the polymer. The reference particle for the API is the centre of mass
    while the reference particle for the polymer is the user's choice."""

    api_simulation_extractor(
        topology_file=topology_file,
        trajectory_file=trajectory_file,
        api_residue_name=api_residue_name,
        output_filename=output_filename,
        start_time=start_time,
    )
    api_universe = dummy_universe(api_simulation_file=output_filename)

    polymer_atom_selection, timestep = polymer_atom_selection, timestep = (
        polymer_atom_extraction(
            polymer_topology_file=polymer_topology_file,
            polymer_trajectory_file=polymer_trajectory_file,
            atom_name=polymer_atom_name,
        )
    )

    # This timestep could have been calculated with information from the
    # API simulation file if we included it but there is no need as they
    # will both be the same as the polymer information just comes from
    # the main simulation information.
    start_frame = int(start_time / timestep)

    api_polymer_rdf = InterRDF(
        api_universe.atoms, polymer_atom_selection, nbins=120, range=(0, 40)
    )
    api_polymer_rdf.run(start=start_frame)

    return api_polymer_rdf.results.bins, api_polymer_rdf.results.rdf


def rdf_plotter(rdf_bins, rdf_values, rdf_type):
    """Minimal plotting function for radial distribution functions."""
    plt.figure(figsize=(8, 6))
    plt.xlabel("Distance (Angstrom)", fontsize=15)
    plt.ylabel("RDF Value", fontsize=15)
    plt.tick_params(labelsize=15)
    plt.grid()
    plt.plot(rdf_bins, rdf_values)
    plt.savefig(f"{rdf_type}_rdf.png", dpi=300)
    return
