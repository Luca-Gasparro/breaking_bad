from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt
import sys
import os
import numpy as np

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
        lx, ly, lz = ts.dimension[:3]
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

    api_simulation_extractor(
        topology_file=topology_file,
        trajectory_file=trajectory_file,
        api_residue_name=api_residue_name,
        output_filename=output_filename,
        start_time=start_time,
    )
    api_universe = dummy_universe(api_simulation_file=output_filename)

    # Obtaining a safe interval for calculating the rdf
    api_universe_load = api_universe.universe
    half_box_lengths = []
    for ts in api_universe_load.trajectory:
        lx, ly, lz = ts.dimensions[:3]
        half_box_lengths.append(min(lx, ly, lz) / 2.0)
    # Safest range to use: minimum of half box lengths
    max_range = min(half_box_lengths)

    _, polymer_atom_selection, timestep = polymer_atom_selection, timestep = (
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
            api_universe.atoms, polymer_atom_selection, nbins=120, range=(0, max_range)
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


def api_api_rdf_manual(simulation_information_filename, bin_width, instantaneous):
    """Compute API-API rdf manually from saved centres of mass."""

    api_coms, box_lengths = load_simulation(simulation_information_filename)
    number_of_frames, number_of_molecules, _ = api_coms.shape

    # Setting up bins for rdf
    mean_box = np.mean(box_lengths, axis=0)
    r_max = 0.5 * np.mean(mean_box)
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

    else:
        # Average density normalisation
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


bin_man_inst, g_r_man_inst = api_api_rdf_manual("test_api.npz", 0.5, True)
bin_man, g_r_man = api_api_rdf_manual("test_api.npz", 0.5, False)
bin_mda, g_r_mda = api_api_rdf(
    "dry_nvt.tpr", "dry_nvt_trim_whole.xtc", "NAP", "test_api.npz", 2000, [1], 75
)
plt.figure(figsize=(8, 6))
plt.plot(bin_man, g_r_man, label="manual")
plt.plot(bin_mda, g_r_mda[0], label="mda")
plt.plot(bin_man_inst, g_r_man_inst, label="inst")
plt.xlabel("distace / ang")
plt.ylabel("rdf")
plt.title("manual api-api rdf test")
plt.legend()
plt.savefig("manual_api-api_rdf_test.png", dpi=300)


bin_mda_npt, g_r_mda_npt = api_api_rdf(
    "dry_npt.tpr", "dry_npt_trim_whole.xtc", "NAP", "test_api_npt.npz", 2000, [1], 75
)
bin_man_npt, g_r_ma_npt = api_api_rdf_manual("test_api_npt.npz", 0.5, False)
bin_man_inst_npt, g_r_man_inst_npt = api_api_rdf_manual("test_api_npt.npz", 0.5, True)
plt.figure(figsize=(8, 6))
plt.plot(bin_man, g_r_man, label="average density norm")
plt.plot(bin_mda_npt, g_r_mda_npt[0], label="mda")
plt.plot(bin_man_inst_npt, g_r_man_inst_npt, label="inst")
plt.xlabel("distace / ang")
plt.ylabel("rdf")
plt.title("manual api-api rdf test npt")
plt.legend()
plt.savefig("manual_api-api_rdf_test_npt.png", dpi=300)
