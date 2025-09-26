# Functions related to calculating the diffusion coefficient

import os
import re
import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def traj_organiser_300k(directory, is_dry):
    """Organises trajectory files numerically by the number in their name."""

    selector = "dry" if is_dry else "wet"
    files = []

    for filename in os.listdir(directory):
        if filename.endswith(".xtc") and selector in filename:
            # Match e.g. dry_traj202950.xtc or segment202950.xtc
            # This is because I am lazy sometimes with what I name these cut down trajectories
            # I will need to rename the trajectories something better soon to indicate that
            # they are at 300 K.
            match = re.search(r"(?:traj|segment)(\d+)\.xtc", filename)
            if match:
                files.append((int(match.group(1)), os.path.join(directory, filename)))

    # Sort by the extracted number
    files.sort(key=lambda x: x[0])

    return [f[1] for f in files]


def msd_calculator(
    topology_file, sorted_trajectory_array, api_residue_name, msd_cache_file_name
):
    """Calculates MSDs from given trajectory array, storing them in an `.npz` file.
    Need to come up with a better way to relate this to temperature but first focus on
    diffusion."""

    # Calculation is lengthy, so first check if it is already done. If already done
    # skip it
    if os.path.exists(msd_cache_file_name):
        print(f"Loading cached MSD data from {msd_cache_file_name}")
        data = np.load(msd_cache_file_name, allow_pickle=True)
        return data["msd_array"], data["lagtimes_array"]

    msds = []
    lagtimes = []

    # Calculating MSD and lagtime for each of the given trajectories
    # Use Fast Fourier Transforms to increase calculation speed

    for traj in sorted_trajectory_array:
        u = mda.Universe(topology_file, traj)
        mean_square_displacement = msd.EinsteinMSD(
            u, select=f"resname {api_residue_name}", msd_tpe="xyz", fft=True
        )
        mean_square_displacement.run()

        # MDAnalysis works in Angstrom - we need to go from squared Angstroms to squared cm
        mean_square_displacement_cm2 = (
            mean_square_displacement.results.timeseries / 1e16
        )

        # Calculate the lagtimes - these are in ps
        number_of_frames = mean_square_displacement.n_frames
        timestep = u.trajectory.dt
        lagtime = np.arange(number_of_frames) * timestep

        msds.append(mean_square_displacement_cm2)
        lagtimes.append(lagtime)

    # Save results into `.npz` file - this might be redudant so
    # after getting a convergence plot tidy this up
    np.savez(
        msd_cache_file_name,
        msd_array=np.array(msds, dtype=object),
        lagtimes_array=np.array(lagtimes, dtype=object),
    )
    print(f"Saved MSD data to {msd_cache_file_name}")

    return np.array(msds, dtype=object), np.array(lagtimes, dtype=object)


def msd_300k_plotter(msd_array, lagtime_array, is_dry):
    """Plots the MSD curves relating to the temperature 300 K"""
    wet_label = "Dry" if is_dry else "Wet"
    for i, (lagtimes, msd_vals) in enumerate(zip(lagtime_array, msd_array), start=1):
        plt.figure(figsize=(10, 6))
        plt.plot(lagtimes, msd_vals)
        plt.xlabel("Time (ps)", fontsize=15)
        plt.ylabel(r"MSD (cm$^2$)", fontsize=15)
        plt.title(f"MSD {wet_label}")
        plt.savefig(f"msd_300k_{i}.png", dpi=300)
        plt.close()
    return


def diffusion_coefficients_300k(msd_array, lagtime_array, start_array_ps, end_array_ps):
    """Calculates diffusion coefficients at 300 K for each given trajectory. Need to give
    the fitting windows determined by the start and end array. The end times cannot be greater than
    half the trajectory length due to fading statstics. This function should be used after
    inspecting the relevant MSD plots."""

    D = []

    # Final lagtime of each trajectory will be the associated plotting value
    lagtime_end = []

    # Calculation of the diffusion coefficient using Einstein's relation
    lagtime_end = []
    for msd_vals, lagtimes, start, end in zip(
        msd_array, lagtime_array, start_array_ps, end_array_ps
    ):
        # Convert the lagtimes to seconds
        lagtimes_seconds = lagtimes / 1e12
        start_index = np.searchsorted(lagtimes_seconds, start / 1e12)
        end_index = np.searchsorted(lagtimes_seconds, end / 1e12)
        if end_index > len(lagtimes_seconds):
            end_index = len(lagtimes_seconds)

        # Get the slope of the linear region denoted by the start and end array
        slope, _, _, _, _ = linregress(
            lagtimes_seconds[start_index:end_index], msd_vals[start_index:end_index]
        )
        diff_coeff = slope / 6
        D.append(diff_coeff)
        lagtime_end.append(lagtimes[-1])

    return np.array(D), np.array(lagtime_end)


def plot_diff_time_300k_inset(diff_array, lagtime_end_array, is_dry):
    """Plots diffusion against simulation time at 300 K, with an inset for last points."""
    wet_label = "Dry" if is_dry else "Wet"
    fig, ax = plt.subplots(figsize=(10, 6))

    # Main plot
    ax.plot(lagtime_end_array, diff_array, marker="o")
    ax.set_xlabel("Time (ps)", fontsize=20)
    ax.set_ylabel("Diffusion coefficient", fontsize=20)
    ax.tick_params(labelsize=20)

    inset_ax = inset_axes(
        ax,
        width=5,
        height=3,  # inches, not percentages
        loc="upper right",
        bbox_to_anchor=(0.95, 0.9),  # (x=1 fixed to right edge, y<1 shifts down)
        bbox_transform=ax.transAxes,
        borderpad=0,
    )
    n_points = max(5, len(lagtime_end_array) // 10)  # at least 5 points or 10% of data
    inset_ax.plot(lagtime_end_array[-n_points:], diff_array[-n_points:], marker="o")
    inset_ax.set_title("Zoom (last points)", fontsize=10)

    # Optional: tighter limits around inset data
    inset_ax.set_xlim(
        min(lagtime_end_array[-n_points:]), max(lagtime_end_array[-n_points:])
    )
    inset_ax.set_ylim(min(diff_array[-n_points:]), max(diff_array[-n_points:]))

    plt.savefig(f"diffusion_vs_time_{wet_label.lower()}.png", bbox_inches="tight")
    return


trajs = traj_organiser_300k("/storage/chem/phuqdw/breaking-bad/diffusion", True)

msd_300k, lagtimes_300k = msd_calculator(
    "dry_cooling_ramp.tpr", trajs, "NAP", "conv_test2.npz"
)

msd_300k_plotter(msd_300k, lagtimes_300k, is_dry=True)

start_fit = [
    300,
    300,
    500,
    1000,
    1000,
    1000,
    1000,
    1000,
    1500,
    2000,
    2000,
    2000,
    2000,
    2000,
    2000,
    2000,
    2000,
    2500,
    5000,
    2500,
    2500,
    6000,
    6000,
    6000,
    6000,
    6000,
    7000,
    7000,
    7500,
    7500,
]

end_fit = [
    500,
    750,
    1500,
    2000,
    2500,
    3000,
    3000,
    4000,
    4000,
    5000,
    5000,
    5000,
    6000,
    6000,
    7000,
    8000,
    8000,
    9000,
    9500,
    10000,
    10000,
    11000,
    11500,
    12000,
    12500,
    13000,
    13500,
    14000,
    14500,
    15000,
]

D_300, lagtime_end = diffusion_coefficients_300k(
    msd_300k, lagtimes_300k, start_fit, end_fit
)
plot_diff_time_300k_inset(D_300, lagtime_end, True)
