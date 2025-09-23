# Functions related to calculating the diffusion coefficient

import os
import re


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
