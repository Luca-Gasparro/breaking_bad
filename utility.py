# Utility functions
import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader


def api_simulation_extractor(
    topology_file, trajectory_file, api_residue_name, output_filename, start_time
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
    traj = mda.Universe(topology_file, trajectory_file)
    timestep = traj.trajectory.dt
    start_frame = int(start_time / timestep)
    number_of_frames = len(traj.trajectory) - start_frame

    # Select residues of API
    api = traj.select_atoms(f"resname {api_residue_name}")
    api_resids = sorted(set(api.resids))
    number_of_api = len(api_resids)

    # COM and box length per frame extraction
    api_coms = np.zeros((number_of_frames, number_of_api, 3))
    box_lengths = np.zeros((number_of_frames, 3))
    for i, ts in enumerate(traj.trajectory[start_frame:]):
        box_lengths[i] = ts.dimensions[:3]
        for j, residue in enumerate(api_resids):
            mol = traj.select_atoms(f"resname {api_residue_name} and resid {residue}")
            api_coms[i, j] = mol.center_of_mass()

    # Save the results to a file to avoid repeat computation
    np.savez_compressed(output_filename, api_coms=api_coms, box_lengths=box_lengths)
    return


def load_simulation(api_simulation_file):
    """Loads API COM file"""

    # Load the dta stored in the COM file
    simulation_file_data = np.load(api_simulation_file, allow_pickle=True)
    return simulation_file_data["api_coms"], simulation_file_data["box_lengths"]


def dummy_universe(api_simulation_file):
    """Creates a minimal dummy Universe for the API COMS. This
    allows radial distribution functions to be calculated with MDAnalysis"""

    if not api_simulation_file.endswith(".npz"):
        api_simulation_file += ".npz"

    api_com_array, box_lengths = load_simulation(
        api_simulation_file=api_simulation_file
    )
    # Extract the number of frames and atoms
    number_of_frames, number_of_atoms, _ = api_com_array.shape

    # Construct the box dimensions for each frame
    # First three entries are the Lx, Ly, Lz lengths and the last
    # three are the angles - all 90 degrees.
    box_dims = np.zeros((number_of_frames, 6))
    box_dims[:, :3] = box_lengths
    box_dims[:, 3:] = 90.0

    # Make the emtpy universe and attach the trajectory and box dimensions
    # with the MemoryReader function
    dummy_universe = mda.Universe.empty(
        n_atoms=number_of_atoms,
        n_residues=number_of_atoms,
        atom_resindex=np.arange(number_of_atoms),
        trajectory=True,
    )
    dummy_universe.trajectory = MemoryReader(api_com_array, dimensions=box_dims)

    return dummy_universe


def polymer_atom_extraction(polymer_topology_file, polymer_trajectory_file, atom_name):
    """Extracts a representative atom that will be used as a reference particle
    for the polymer radial distribution funciton. The topology and trajectory files are based purely on
    the combined polymer index. This is done to make atom selection easier."""

    # Use MDAnalysis to load in the trajectory
    polymer_traj = mda.Universe(polymer_topology_file, polymer_trajectory_file)
    polymer_timestep = polymer_traj.trajectory.dt

    # Atom selection
    polymer_atom_selection = polymer_traj.select_atoms(f"name {atom_name}")

    return polymer_traj, polymer_atom_selection, polymer_timestep
