from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt
import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utility import api_simulation_extractor, dummy_universe, polymer_atom_extraction


def api_api_rdf(
    topology_file,
    trajectory_file,
    api_residue_name,
    simulation_information_filename,
    start_time,
    frame_strides,
):
    """Computes and plots the RDF for API-API. Uses the centres of mass of the API
    as reference particles."""

    # Extract API COMs and create a dummy universe for them
    api_simulation_extractor(
        topology_file=topology_file,
        trajectory_file=trajectory_file,
        api_residue_name=api_residue_name,
        output_filename=simulation_information_filename,
        start_time=start_time,
    )
    api_universe = dummy_universe(api_simulation_file=simulation_information_filename)

    rdf_values = []
    bins = None

    # Calculate the API-API RDF
    for stride in frame_strides:
        api_rdf = InterRDF(
            api_universe.atoms,
            api_universe.atoms,
            nbins=100,
            range=(0, 42.5),
            exclusion_block=(1, 1),
        )
        api_rdf.run(step=stride)
        rdf_values.append(api_rdf.results.rdf)
        if bins is None:
            bins = api_rdf.results.bins

    # Return information needed to plot - the bins and the RDF
    return np.array(bins), np.array(rdf_values)


def polymer_polymer_rdf(
    polymer_topology_file,
    polymer_trajectory_file,
    polymer_atom_name,
    start_time,
    frame_strides,
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

    rdf_values = []
    bins = None

    # Calculate the polymer-polymer RDF
    for stride in frame_strides:
        polymer_rdf = InterRDF(
            polymer_atom_selection,
            polymer_atom_selection,
            nbins=150,
            range=(0, 42.5),
            exclusion_block=(1, 1),
        )
        polymer_rdf.run(start=start_frame, step=stride)
        rdf_values.append(polymer_rdf.results.rdf)
        if bins is None:
            bins = polymer_rdf.results.bins

    # Return information needed to plot - the bins and the RDF
    return np.array(bins), np.array(rdf_values)


def api_polymer_rdf(
    topology_file,
    trajectory_file,
    api_residue_name,
    output_filename,
    polymer_topology_file,
    polymer_trajectory_file,
    polymer_atom_name,
    start_time,
    frame_strides,
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

    rdf_values = []
    bins = None

    for stride in frame_strides:
        api_polymer_rdf = InterRDF(
            api_universe.atoms, polymer_atom_selection, nbins=120, range=(0, 40)
        )
        api_polymer_rdf.run(start=start_frame, step=stride)
        rdf_values.append(api_polymer_rdf.results.rdf)
        if bins is None:
            bins = api_polymer_rdf.results.bins

    return np.array(bins), np.array(rdf_values)


def rdf_plotter(rdf_bins, rdf_values, rdf_type, frame_strides):
    """Plotting RDF curves for different frame strides"""
    plt.figure(figsize=(8, 6))
    plt.xlabel("Distance (Angstrom)", fontsize=15)
    plt.ylabel("RDF Value", fontsize=15)
    plt.tick_params(labelsize=15)
    plt.grid()
    for i, stride in enumerate(frame_strides):
        plt.plot(rdf_bins, rdf_values[i], label=f"stride {stride}")
    plt.legend()
    plt.savefig(f"{rdf_type}_rdf.png", dpi=300)
    return


frame_strides = np.array([1, 2, 4, 6, 8, 10])
api_bins, api_rdf = api_api_rdf(
    "dry_nvt.tpr",
    "dry_nvt_trim_whole.xtc",
    "NAP",
    "test_api",
    2000,
    frame_strides=frame_strides,
)
rdf_plotter(api_bins, api_rdf, rdf_type="api-api", frame_strides=frame_strides)

poly_bins, poly_rdf = polymer_polymer_rdf(
    "dry_nvt_polymer.tpr",
    "dry_nvt_trim_polymers.xtc",
    "N1",
    2000,
    frame_strides=frame_strides,
)
rdf_plotter(poly_bins, poly_rdf, rdf_type="poly-poly", frame_strides=frame_strides)

api_poly_bins, api_poly_rdf = api_polymer_rdf(
    "dry_nvt.tpr",
    "dry_nvt_trim_whole.xtc",
    "NAP",
    "test_api",
    "dry_nvt_polymer.tpr",
    "dry_nvt_trim_polymers.xtc",
    "N1",
    2000,
    frame_strides=frame_strides,
)

rdf_plotter(api_poly_bins, api_poly_rdf, "api-poly", frame_strides=frame_strides)
