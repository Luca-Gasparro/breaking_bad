from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt
import sys
import os
import numpy as np
import MDAnalysis as mda

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utility import (
    api_simulation_extractor,
    load_simulation,
    dummy_universe,
    polymer_atom_extraction,
)


def api_api_rdf(
    topology_file,
    trajectory_file,
    api_residue_name,
    simulation_information_filename,
    start_time,
    frame_strides,
    nbins,
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

    # Obtaining a safe interval for calculating the rdf
    api_universe_load = api_universe.universe
    half_box_lengths = []
    for ts in api_universe_load.trajectory:
        lx, ly, lz = ts.dimensions[:3]
        half_box_lengths.append(min(lx, ly, lz) / 2.0)
    # Safest range to use: minimum of half box lengths
    max_range = min(half_box_lengths)

    rdf_values = []
    bins = None

    # Calculate the API-API RDF
    for stride in frame_strides:
        api_rdf = InterRDF(
            api_universe.atoms,
            api_universe.atoms,
            nbins=nbins,
            range=(0, max_range),
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
    polymer_trajectory, polymer_atom_selection, timestep = polymer_atom_extraction(
        polymer_topology_file=polymer_topology_file,
        polymer_trajectory_file=polymer_trajectory_file,
        atom_name=polymer_atom_name,
    )
    start_frame = int(start_time / timestep)

    # Safe calculation range
    polymer_trajectory_uni = polymer_trajectory.universe
    half_box_lengths = []
    for ts in polymer_trajectory_uni.trajectory:
        lx, ly, lz = ts.dimensions[:3]
        half_box_lengths.append(min(lx, ly, lz) / 2.0)
    max_range = min(half_box_lengths)

    rdf_values = []
    bins = None

    # Calculate the polymer-polymer RDF
    for stride in frame_strides:
        polymer_rdf = InterRDF(
            polymer_atom_selection,
            polymer_atom_selection,
            nbins=150,
            range=(0, max_range),
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

    # Extract API COM data starting at start_time
    api_simulation_extractor(
        topology_file=topology_file,
        trajectory_file=trajectory_file,
        api_residue_name=api_residue_name,
        output_filename=output_filename,
        start_time=start_time,
    )
    api_universe = dummy_universe(api_simulation_file=output_filename)

    # Determine safe RDF range
    half_box_lengths = []
    for ts in api_universe.trajectory:
        lx, ly, lz = ts.dimensions[:3]
        half_box_lengths.append(min(lx, ly, lz) / 2.0)
    max_range = min(half_box_lengths)

    # Load polymer universe and trim to start_time
    polymer_universe = mda.Universe(polymer_topology_file, polymer_trajectory_file)
    timestep = polymer_universe.trajectory.dt
    start_frame = int(start_time / timestep)

    polymer_universe.trajectory[start_frame:]  # trim polymer frames to start_time
    polymer_atom_selection = polymer_universe.select_atoms(f"name {polymer_atom_name}")

    rdf_values = []
    bins = None

    for stride in frame_strides:
        rdf_calc = InterRDF(
            api_universe.atoms,
            polymer_atom_selection,
            nbins=75,
            range=(0, max_range),
        )
        rdf_calc.run(
            start=0, step=stride
        )  # start=0 because both universes are now aligned
        rdf_values.append(rdf_calc.results.rdf)
        if bins is None:
            bins = rdf_calc.results.bins

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


def api_api_rdf_manual(simulation_information_filename, bin_width, instantaneous):
    """Compute API-API rdf manually from saved centres of mass."""

    api_coms, box_lengths = load_simulation(simulation_information_filename)
    number_of_frames, number_of_molecules, _ = api_coms.shape

    # Safest cut-off: half the minimum box length for all frames
    half_box_lengths = np.min(box_lengths, axis=1) / 2.0
    r_max = np.min(half_box_lengths)
    bins = np.arange(0, r_max + bin_width, bin_width)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    shell_volumes = (4 / 3) * np.pi * (bins[1:] ** 3 - bins[:-1] ** 3)

    # Unique number of pairs
    n_pairs = number_of_molecules * (number_of_molecules - 1) / 2

    # Per frame normalisation
    if instantaneous:
        gr_frames = []
        for i in range(number_of_frames):
            frame = api_coms[i]
            box = box_lengths[i]
            # Calculate pairwise displacement vectors
            delta = frame[np.newaxis, :, :] - frame[:, np.newaxis, :]
            # Apply the minimum image convention to ensure we get the correct displacement vector
            delta -= np.rint(delta / box) * box
            # Calculate distances and avoid double counting as well as self-self ditances
            distances = np.sqrt(np.sum(delta**2, axis=-1))
            dists = distances[np.triu_indices(number_of_molecules, k=1)]
            hist, _ = np.histogram(dists, bins=bins)

            vol = np.prod(box)
            expected = n_pairs * (shell_volumes / vol)
            gr_frames.append(hist / expected)

        g_r = np.mean(gr_frames, axis=0)

    # Average density normalisation
    else:
        rdf_accumulation = np.zeros(len(bin_centers))
        for i in range(number_of_frames):
            frame = api_coms[i]
            box = box_lengths[i]
            # Calculate pairwise displacement vectors
            delta = frame[np.newaxis, :, :] - frame[:, np.newaxis, :]
            # Apply the minimum image convention to ensure we get the correct displacement vector
            delta -= np.rint(delta / box) * box
            # Calculate distances and avoid double counting as well as self-self ditances
            distances = np.sqrt(np.sum(delta**2, axis=-1))
            dists = distances[np.triu_indices(number_of_molecules, k=1)]
            hist, _ = np.histogram(dists, bins=bins)
            rdf_accumulation += hist

        mean_vol = np.mean(np.prod(box_lengths, axis=1))
        expected = n_pairs * number_of_frames * (shell_volumes / mean_vol)
        g_r = rdf_accumulation / expected

    return bin_centers, g_r


def polymer_polymer_rdf_manual(
    polymer_topology_file,
    polymer_trajectory_file,
    reference_particle_name,
    bin_width,
    instantaneous,
    start_time,
):
    """Compute polymer-polymer rdf manually using a chosen reference atom."""

    # Obtain reference particle positions and box lengths from time of equilbration
    polymer_trajectory, polymer_atom_selection, timestep = polymer_atom_extraction(
        polymer_topology_file=polymer_topology_file,
        polymer_trajectory_file=polymer_trajectory_file,
        atom_name=reference_particle_name,
    )
    start_frame = int(start_time / timestep)

    polymer_atom_positions = []
    box_lengths = []

    for ts in polymer_trajectory.trajectory[start_frame:]:
        polymer_atom_positions.append(polymer_atom_selection.positions.copy())
        box_lengths.append(ts.dimensions[:3].copy())

    polymer_atom_positions = np.array(polymer_atom_positions)
    box_lengths = np.array(box_lengths)
    number_of_frames = len(polymer_atom_positions)

    # Set up the bins using a safe cut off: half the minimum box length for all frames
    half_box_lengths = np.min(box_lengths, axis=1) / 2.0
    r_max = np.min(half_box_lengths)
    bins = np.arange(0, r_max + bin_width, bin_width)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    shell_volumes = (4 / 3) * np.pi * (bins[1:] ** 3 - bins[:-1] ** 3)

    # Unique number of pairs
    number_of_atoms = len(polymer_atom_selection)
    n_pairs = number_of_atoms * (number_of_atoms - 1) / 2

    # Per frame normalisation
    if instantaneous:
        gr_frames = []
        for i in range(number_of_frames):
            frame = polymer_atom_positions[i]
            box = box_lengths[i]
            # Calculate pairwise displacement vectors
            delta = frame[np.newaxis, :, :] - frame[:, np.newaxis, :]
            # Apply the minimum image convention to ensure we get the correct displacement vector
            delta -= np.rint(delta / box) * box
            # Calculate distances and avoid double counting as well as self-self ditances
            distances = np.sqrt(np.sum(delta**2, axis=-1))
            dists = distances[np.triu_indices(number_of_atoms, k=1)]
            hist, _ = np.histogram(dists, bins=bins)

            vol = np.prod(box)
            expected = n_pairs * (shell_volumes / vol)
            gr_frames.append(hist / expected)

        g_r = np.mean(gr_frames, axis=0)

    # Average density normalisation
    else:
        rdf_accumulation = np.zeros(len(bin_centers))
        for i in range(number_of_frames):
            frame = polymer_atom_positions[i]
            box = box_lengths[i]
            # Calculate pairwise displacement vectors
            delta = frame[np.newaxis, :, :] - frame[:, np.newaxis, :]
            # Apply the minimum image convention to ensure we get the correct displacement vector
            delta -= np.rint(delta / box) * box
            # Calculate distances and avoid double counting as well as self-self ditances
            distances = np.sqrt(np.sum(delta**2, axis=-1))
            dists = distances[np.triu_indices(number_of_atoms, k=1)]
            hist, _ = np.histogram(dists, bins=bins)
            rdf_accumulation += hist

        mean_vol = np.mean(np.prod(box_lengths, axis=1))
        expected = n_pairs * number_of_frames * (shell_volumes / mean_vol)
        g_r = rdf_accumulation / expected

    return bin_centers, g_r


def polymer_api_rdf_manual(
    polymer_topology_file,
    polymer_trajectory_file,
    reference_particle_name,
    simulation_information_filename,
    bin_width,
    instantaneous,
    start_time,
):
    """Compute polymer–API RDF manually between chosen polymer reference atoms and API COMs.
    Combines features from both RDF functions."""

    # Polymer atom extraction
    polymer_trajectory, polymer_atom_selection, timestep = polymer_atom_extraction(
        polymer_topology_file=polymer_topology_file,
        polymer_trajectory_file=polymer_trajectory_file,
        atom_name=reference_particle_name,
    )
    start_frame = int(start_time / timestep)

    polymer_atom_positions = []
    box_lengths = []

    for ts in polymer_trajectory.trajectory[start_frame:]:
        polymer_atom_positions.append(polymer_atom_selection.positions.copy())
        box_lengths.append(ts.dimensions[:3].copy())

    polymer_atom_positions = np.array(polymer_atom_positions)
    box_lengths = np.array(box_lengths)
    number_of_frames, number_of_polymer_atoms, _ = polymer_atom_positions.shape

    # API centres of mass
    api_coms, api_box_lengths = load_simulation(simulation_information_filename)
    _, number_of_api, _ = api_coms.shape

    # Check if number of frames are consistent
    if api_coms.shape[0] != number_of_frames:
        raise ValueError("Mismatch in frames between polymer trajectory and API COMs.")

    # Setting up the RDF bins
    half_box_lengths = np.min(box_lengths, axis=1) / 2.0
    r_max = np.min(half_box_lengths)
    bins = np.arange(0, r_max + bin_width, bin_width)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    shell_volumes = (4 / 3) * np.pi * (bins[1:] ** 3 - bins[:-1] ** 3)

    # Total number of cross pairs per frame
    n_pairs = number_of_polymer_atoms * number_of_api

    # Compute RDF
    # Instantaneous normalisation
    if instantaneous:
        gr_frames = []
        for i in range(number_of_frames):
            poly_frame = polymer_atom_positions[i]  # (N_poly, 3)
            api_frame = api_coms[i]  # (N_api, 3)
            box = box_lengths[i]

            # Pairwise displacements: shape (N_poly, N_api, 3) - also applying
            # minimum image convention
            delta = poly_frame[:, np.newaxis, :] - api_frame[np.newaxis, :, :]
            delta -= np.rint(delta / box) * box

            distances = np.sqrt(np.sum(delta**2, axis=-1)).ravel()

            hist, _ = np.histogram(distances, bins=bins)

            vol = np.prod(box)
            expected = n_pairs * (shell_volumes / vol)
            gr_frames.append(hist / expected)

        g_r = np.mean(gr_frames, axis=0)

    # Average density normalisation
    else:
        rdf_accumulation = np.zeros(len(bin_centers))
        for i in range(number_of_frames):
            poly_frame = polymer_atom_positions[i]
            api_frame = api_coms[i]
            box = box_lengths[i]

            delta = poly_frame[:, np.newaxis, :] - api_frame[np.newaxis, :, :]
            delta -= np.rint(delta / box) * box
            distances = np.sqrt(np.sum(delta**2, axis=-1)).ravel()

            hist, _ = np.histogram(distances, bins=bins)
            rdf_accumulation += hist

        mean_vol = np.mean(np.prod(box_lengths, axis=1))
        expected = n_pairs * number_of_frames * (shell_volumes / mean_vol)
        g_r = rdf_accumulation / expected

    return bin_centers, g_r
