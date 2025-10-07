# Utility functions
import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader


def api_simulation_extractor(
    topology_file, trajectory_file, api_residue_name, output_filename
):
    """Extract centres-of-mass of the API."""

    # Make sure the output has a .npz extension
    if not output_filename.endswith(".npz"):
        output_filename += ".npz"

    # Check for existence of filename. If it exists, skip the recalculation
    if os.path.exists(output_filename):
        print(f"API COM file {output_filename} already exists, skipping extraction.")
        return output_filename

    # Use MDAnalysis to load in the trajectories
    traj = mda.Universe(topology_file, trajectory_file)

    # Extract the number of frames and the timestep
    number_of_frames = len(traj.trajectory)
    timestep = traj.trajectory.dt

    # Select residues of API
    api = traj.select_atoms(f"resname {api_residue_name}")
    api_resids = sorted(set(api.resids))
    number_of_api = len(api_resids)

    # COM and box length per frame extraction
    api_coms = np.zeros((number_of_frames, number_of_api, 3))
    box_lengths = np.zeros((number_of_frames, 3))
    for i, ts in enumerate(traj.trajectory):
        box_lengths[i] = ts.dimensions[:3]
        for j, residue in enumerate(api_resids):
            mol = traj.select_atoms(f"resname {api_residue_name} and resid {residue}")
            api_coms[i, j] = mol.center_of_mass()

    # Save the results to a file to avoid repeat computation
    np.savez_compressed(
        output_filename, api_coms=api_coms, box_lengths=box_lengths, timestep=timestep
    )
    return


def load_simulation(api_simulation_file):
    """Loads API COM file"""

    # Load the dta stored in the COM file
    api_simulation_file_data = np.load(api_simulation_file, allow_pickle=True)
    return (
        api_simulation_file_data["api_coms"],
        api_simulation_file_data["box_lengths"],
        api_simulation_file_data["timestep"],
    )


def dummy_universe(api_simulation_file):
    """Creates a minimal dummy Universe for the API COMS. This
    allows radial distribution functions to be calculated with MDAnalysis"""

    if not api_simulation_file.endswith(".npz"):
        api_simulation_file += ".npz"

    api_com_array, box_lengths, _ = load_simulation(
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


def combined_universe(
    polymer_topology_file,
    polymer_trajectory_file,
    api_simulation_file,
    polymer_atom_name,
    api_atom_name,
):
    polymer_u = mda.Universe(polymer_topology_file, polymer_trajectory_file)
    polymer_coords = np.array([ts.positions.copy() for ts in polymer_u.trajectory])
    box_lengths = np.array([ts.dimensions.copy() for ts in polymer_u.trajectory])
    n_frames, n_poly, _ = polymer_coords.shape

    api_com_array, api_box_lengths, api_dt = load_simulation(api_simulation_file)
    n_frames_api, n_api, _ = api_com_array.shape

    if n_frames != n_frames_api:
        raise ValueError(f"Frame mismatch: polymer {n_frames}, API {n_frames_api}")
    if not np.allclose(box_lengths[:, :3], api_box_lengths, rtol=1e-5):
        raise ValueError(
            "Box lengths do not match between polymer and API trajectories"
        )

    combined_coords = np.concatenate([polymer_coords, api_com_array], axis=1)
    n_total = n_poly + n_api

    combined_u = mda.Universe.empty(
        n_atoms=n_total,
        n_residues=n_total,
        atom_resindex=np.arange(n_total),
        trajectory=True,
    )

    # Give unique names so selections work
    atom_names = [f"{polymer_atom_name}{i}" for i in range(n_poly)] + [
        f"{api_atom_name}{j}" for j in range(n_api)
    ]
    combined_u.add_TopologyAttr("name", atom_names)

    combined_u.trajectory = MemoryReader(
        combined_coords, dimensions=box_lengths, dt=api_dt
    )
    return combined_u
