# Functions related to calculating the diffusion coefficient

import os
import re
import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd
import matplotlib.pyplot as plt


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

    return np.array(msd, dtype=object), np.array(lagtimes, dtype=object)


def msd_300k_plotter(msd_array, lagtime_array, is_dry):
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


trajs = traj_organiser_300k("/storage/chem/phuqdw/breaking-bad/diffusion", True)

msd_300k, lagtimes_300k = msd_calculator(
    "dry_cooling_ramp.tpr", trajs, "NAP", "conv_test.npz"
)

msd_300k_plotter(msd_300k, lagtimes_300k, is_dry=True)
