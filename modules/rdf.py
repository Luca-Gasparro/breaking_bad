import matplotlib.pyplot as plt
import sys
import os
import numpy as np


sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


from .utility import (
    load_simulation,
    polymer_atom_extraction,
)

from .analysis import average_box_length, parse_mdp


def unified_rdf(
    mode,
    bin_width,
    instantaneous,
    start_time,
    api_information_file,
    polymer_topology_file,
    polymer_trajectory_file,
    reference_particle_name,
):
    """Unified radial distribution function. The three modes are:
    1) mode == "api-api" which uses api_information_file,
    2) mode == "polymer-polymer" which uses  use polymer_topology_file, polymer_trajectory_file, reference_particle_name, and
    3) mode == "polymer-api" which uses all the above.

    The function is capaable of producing either an instantaneous density normalisation or an average density normalisation
    """

    if mode == "api-api":
        # Load API coms
        api_coms, box_lenths, api_timestep = load_simulation(
            api_simulation_file=api_information_file
        )
        api_start_frame = int(start_time / api_timestep)
        # Get positions of the API
        positions_A = api_coms[api_start_frame:]
        box_lengths = box_lenths[api_start_frame:]
        # A-A RDF
        positions_B = None
        number_of_frames = len(positions_A)

    elif mode == "polymer-polymer":
        # Extract information about the polymer
        polymer_trajectory, polymer_atom_selection, polymer_timestep = (
            polymer_atom_extraction(
                polymer_topology_file=polymer_topology_file,
                polymer_trajectory_file=polymer_trajectory_file,
                atom_name=reference_particle_name,
            )
        )
        polymer_start_frame = int(start_time / polymer_timestep)

        # Polymer atom positions
        polymer_atom_positions = []
        box_lengths_list = []
        for ts in polymer_trajectory.trajectory[polymer_start_frame:]:
            polymer_atom_positions.append(polymer_atom_selection.positions.copy())
            box_lengths_list.append(ts.dimensions[:3].copy())
        positions_A = np.array(polymer_atom_positions)
        box_lengths = np.array(box_lengths_list)
        # A-A RDF
        positions_B = None
        number_of_frames = len(positions_A)

    elif mode == "polymer-api":
        # polymer extraction
        polymer_trajectory, polymer_atom_selection, polymer_timestep = (
            polymer_atom_extraction(
                polymer_topology_file=polymer_topology_file,
                polymer_trajectory_file=polymer_trajectory_file,
                atom_name=reference_particle_name,
            )
        )
        polymer_start_frame = int(start_time / polymer_timestep)

        polymer_atom_positions = []
        box_lengths_list = []
        for ts in polymer_trajectory.trajectory[polymer_start_frame:]:
            polymer_atom_positions.append(polymer_atom_selection.positions.copy())
            box_lengths_list.append(ts.dimensions[:3].copy())

        # Polymer atoms positions
        positions_A = np.array(polymer_atom_positions)
        box_lengths = np.array(box_lengths_list)

        # API COMs (use same slicing approach as your polymer_api_rdf_manual: slice by polymer start_frame)
        api_coms, _, api_timestep = load_simulation(api_information_file)
        api_coms = api_coms[polymer_start_frame:]
        positions_B = api_coms

        # Ensuring number of frames are consistent
        if positions_B.shape[0] != positions_A.shape[0]:
            raise ValueError(
                "Mismatch in frames between polymer trajectory and API COMs."
            )
        number_of_frames = len(positions_A)

    else:
        raise ValueError(
            "mode must be one of: 'api-api', 'polymer-polymer', 'polymer-api'"
        )

    # RDF calculation
    n_A = positions_A.shape[1]
    use_cross = positions_B is not None

    if use_cross:
        n_B = positions_B.shape[1]
        n_pairs = n_A * n_B
    else:
        n_pairs = n_A * (n_A - 1) / 2

    # Safe maximum radius: half minimum box dimension across frames
    half_box_lengths = np.min(box_lengths, axis=1) / 2.0
    r_max = np.min(half_box_lengths)
    bins = np.arange(0, r_max + bin_width, bin_width)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    shell_volumes = (4 / 3) * np.pi * (bins[1:] ** 3 - bins[:-1] ** 3)

    if instantaneous:
        gr_frames = []
        for i in range(number_of_frames):
            A = positions_A[i]
            B = A if not use_cross else positions_B[i]
            box = box_lengths[i]

            # <-- YOUR exact delta style preserved -->
            delta = A[np.newaxis, :, :] - B[:, np.newaxis, :]
            delta -= np.rint(delta / box) * box

            distances = np.sqrt(np.sum(delta**2, axis=-1))

            if not use_cross:
                distances = distances[np.triu_indices(n_A, k=1)]
            else:
                distances = distances.ravel()

            hist, _ = np.histogram(distances, bins=bins)
            vol = np.prod(box)
            expected = n_pairs * (shell_volumes / vol)
            gr_frames.append(hist / expected)

        g_r = np.mean(gr_frames, axis=0)

    else:
        rdf_accumulation = np.zeros(len(bin_centers))
        for i in range(number_of_frames):
            A = positions_A[i]
            B = A if not use_cross else positions_B[i]
            box = box_lengths[i]

            delta = A[np.newaxis, :, :] - B[:, np.newaxis, :]
            delta -= np.rint(delta / box) * box

            distances = np.sqrt(np.sum(delta**2, axis=-1))

            if not use_cross:
                distances = distances[np.triu_indices(n_A, k=1)]
            else:
                distances = distances.ravel()

            hist, _ = np.histogram(distances, bins=bins)
            rdf_accumulation += hist

        mean_vol = np.mean(np.prod(box_lengths, axis=1))
        expected = n_pairs * number_of_frames * (shell_volumes / mean_vol)
        g_r = rdf_accumulation / expected

    return bin_centers, g_r


def api_api_rdf(
    api_information_file,
    bin_width,
    instantaneous,
    start_time,
):
    return unified_rdf(
        mode="api-api",
        bin_width=bin_width,
        instantaneous=instantaneous,
        start_time=start_time,
        api_information_file=api_information_file,
        polymer_topology_file=None,
        polymer_trajectory_file=None,
        reference_particle_name=None,
    )


def polymer_polymer_rdf(
    polymer_topology_file,
    polymer_trajectory_file,
    reference_particle_name,
    bin_width,
    instantaneous,
    start_time,
):
    return unified_rdf(
        mode="polymer-polymer",
        bin_width=bin_width,
        instantaneous=instantaneous,
        start_time=start_time,
        api_information_file=None,
        polymer_topology_file=polymer_topology_file,
        polymer_trajectory_file=polymer_trajectory_file,
        reference_particle_name=reference_particle_name,
    )


def polymer_api_rdf(
    api_information_file,
    polymer_topology_file,
    polymer_trajectory_file,
    reference_particle_name,
    bin_width,
    instantaneous,
    start_time,
):
    return unified_rdf(
        mode="polymer-api",
        bin_width=bin_width,
        instantaneous=instantaneous,
        start_time=start_time,
        api_information_file=api_information_file,
        polymer_topology_file=polymer_topology_file,
        polymer_trajectory_file=polymer_trajectory_file,
        reference_particle_name=reference_particle_name,
    )


def plot_rdfs(
    api_bin,
    api_rdf,
    poly_bin,
    poly_rdf,
    cross_bin,
    cross_rdf,
    topology_file,
    trajectory_file,
    npt_mdp_file,
    output_dir=None,
):
    """Plots the API-API, Polymer-Polymer and Polymer-API radial distirbution function for a given
    simulation"""
    # If no output folder is provided, use the current working directory
    if output_dir is None:
        output_dir = os.getcwd()
    os.makedirs(output_dir, exist_ok=True)

    # Obtain mean box in the length (cubic box so component doesn't matter)
    mean_box_length = average_box_length(
        topology_file=topology_file, trajectory_file=trajectory_file
    )

    # Extract pressure for saving and plotting purposes
    npt_params = parse_mdp(mdp_file=npt_mdp_file)
    pressure = npt_params["ref-p"]

    # Construct full path for the output figure
    output_path = os.path.join(output_dir, f"npt_{pressure}_bars.png")

    # Plotting
    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 12,
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 12,
            "figure.titlesize": 12,
        }
    )

    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(6, 9), sharex=False, sharey=False
    )
    fig.suptitle(f"RDFs ({pressure} bars, {mean_box_length:.2f} $\\AA$)")

    # API–API
    ax1.plot(api_bin, api_rdf, label="API-API RDF")
    ax1.set_ylabel("g(r)")
    ax1.set_xlabel("r (Å)")
    ax1.grid(True, linestyle=":", linewidth=0.5)
    ax1.legend()

    # Polymer–Polymer
    ax2.plot(poly_bin, poly_rdf, label="Polymer-Polymer RDF")
    ax2.set_ylabel("g(r)")
    ax2.set_xlabel("r (Å)")
    ax2.grid(True, linestyle=":", linewidth=0.5)
    ax2.legend()

    # Polymer–API
    ax3.plot(cross_bin, cross_rdf, label="Polymer-API RDF")
    ax3.set_ylabel("g(r)")
    ax3.set_xlabel("r (Å)")
    ax3.grid(True, linestyle=":", linewidth=0.5)
    ax3.legend()

    plt.tight_layout()

    plt.savefig(output_path, dpi=300)

    return
