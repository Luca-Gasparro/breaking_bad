# Utility functions
import os
import numpy as np
import MDAnalysis as mda


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
