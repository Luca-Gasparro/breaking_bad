# Utility functions
import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader


def api_com_extractor(
    topology_file, trajectroy_file, api_residue_name, output_filename, start_time
):
    """Extract centres-of-mass of the API from the time of equilibration starting."""

    # Make sure the output has a .npz extension
    if not output_filename.endswith(".npz"):
        output_filename += ".npz"

    # Check for existence of filename. If it exists, skip the recalculation
    if os.path.exists(output_filename):
        print(f"API COM file {output_filename} already exists, skipping extraction.")
        return output_filename

    # Use MDAnalysis to load in the trajectories
    traj = mda.Universe(topology_file, trajectroy_file)
    timestep = traj.trajectory.dt
    start_frame = int(start_time / timestep)
    number_of_frames = len(traj.trajectory) - start_frame

    # Select residues of API
    api = traj.select_atoms(f"resname {api_residue_name}")
    api_resids = sorted(set(api.resids))
    number_of_api = len(api_resids)

    # COM extraction
    api_coms = np.zeros((number_of_frames, number_of_api, 3))
    for i, _ in enumerate(traj.trajectory[start_frame:]):
        for j, residue in enumerate(api_resids):
            mol = traj.select_atoms(f"resname {api_residue_name} and resid {residue}")
            api_coms[i, j] = mol.center_of_mass()

    # Save the results to a file to avoid repeat computation
    np.savez_compressed(output_filename, api_coms=api_coms)
    return


def load_api_coms(api_com_file):
    """Loads API COM file"""

    # Load the dta stored in the COM file
    com_file_data = np.load(api_com_file, allow_pickle=True)
    return com_file_data["api_coms"]


def dummy_universe(api_com_array):
    """Creates a minimal dummy Universe for the API COMS. This
    allows radial distribution functions to be calculated with MDAnalysis"""

    # Extract the number of atoms
    _, number_of_atoms, _ = api_com_array.shape

    # Make the emtpy universe and attach the trajectory with
    # the MemoryReader function
    dummy_universe = mda.Universe.empty(
        n_atoms=number_of_atoms,
        n_residues=number_of_atoms,
        atom_resindex=np.arange(number_of_atoms),
        trajectory=True,
    )
    dummy_universe.trajectory = MemoryReader(api_com_array)

    return dummy_universe
