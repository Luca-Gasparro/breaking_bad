import matplotlib.pyplot as plt
import sys
import os
import numpy as np


sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utility import (
    load_simulation,
    polymer_atom_extraction,
)


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


def api_api_rdf_manual(
    simulation_information_filename, bin_width, instantaneous, start_time
):
    """Compute API-API rdf manually from saved centres of mass."""

    api_coms, box_lengths, api_timestep = load_simulation(
        simulation_information_filename
    )
    _, number_of_molecules, _ = api_coms.shape
    api_start_frame = int(start_time / api_timestep)

    # Restrict to frames after start_time
    api_coms = api_coms[api_start_frame:]
    box_lengths = box_lengths[api_start_frame:]
    number_of_frames = len(api_coms)

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

    # API centres of mass after start frame
    api_coms, _, _ = load_simulation(simulation_information_filename)
    api_coms = api_coms[start_frame:]
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
