# Module contaning code to analyse quantities obtained via MD simulations
import numpy as np
import matplotlib.pyplot as plt
import os
import MDAnalysis as mda


def energy_minimisation_pe(
    polymer_xvg, api_polymer_xvg, wet_api_polymer_xvg, box_size, output_dir=None
):
    """Produces three subplots which display how the potential energy of the system
    evolves via the steepest descent method during the energy minimisation stage."""

    # If no output folder is provided, use the current working directory
    if output_dir is None:
        output_dir = os.getcwd()
    os.makedirs(output_dir, exist_ok=True)  # ensure folder exists

    # Construct full path for the output figure
    output_path = os.path.join(output_dir, f"energy_minimisation_{box_size}_ang.png")

    # Extract the data from each .xvg file
    poly_step, poly_potential = np.loadtxt(
        polymer_xvg, comments=("#", "@"), unpack=True
    )
    api_poly_step, api_poly_potential = np.loadtxt(
        api_polymer_xvg, comments=("#", "@"), unpack=True
    )
    wet_api_poly_step, wet_api_poly_potential = np.loadtxt(
        wet_api_polymer_xvg, comments=("#", "@"), unpack=True
    )

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
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(6, 8), sharex=False)
    fig.suptitle(f"Energy minimisation {box_size} $\\AA$")
    ax1.plot(poly_step, poly_potential, label="Polymer")
    ax1.legend()
    ax2.plot(api_poly_step, api_poly_potential, label="API-Polymer")
    ax2.legend()
    ax3.plot(wet_api_poly_step, wet_api_poly_potential, label="API-Polymer (hydrated)")
    ax3.legend()
    for ax in (ax1, ax2, ax3):
        ax.set_ylabel("Potential Energy (kJ/mol)")
        ax.grid(True, linestyle=":", linewidth=0.5)
    ax3.set_xlabel("Minimisation Step")
    plt.tight_layout()
    plt.savefig(
        output_path,
        dpi=300,
    )
    return


def average_box_length(topology_file, trajectory_file):
    """Extract the average box length over a trajectory. This will
    be used when plotting quantities that come from NPT simulations."""

    u = mda.Universe(topology_file, trajectory_file)
    x_lengths = np.array([ts.dimensions[0] for ts in u.trajectory])
    mean_x = x_lengths.mean()

    return mean_x


def npt_analysis(
    energy_xvg,
    temperature_xvg,
    pressure_xvg,
    topology_file,
    trajectory_file,
    output_dir=None,
):
    """Calculates the potential energy, the temperature and the pressure evolution in an NPT
    simulation."""

    # If no output folder is provided, use the current working directory
    if output_dir is None:
        output_dir = os.getcwd()
    os.makedirs(output_dir, exist_ok=True)

    # Obtain mean box in the length (cubic box so component doesn't matter)
    mean_box_length = average_box_length(
        topology_file=topology_file, trajectory_file=trajectory_file
    )
    # Construct full path for the output figure

    output_path = os.path.join(output_dir, f"npt_{mean_box_length}_ang.png")

    # Extract the data from each .xvg file - potential energy, temperature and pressure
    npt_time, npt_potential = np.loadtxt(energy_xvg, comments=("#", "@"), unpack=True)
    _, npt_temperature = np.loadtxt(temperature_xvg, comments=("#", "@"), unpack=True)
    _, npt_pressure = np.loadtxt(pressure_xvg, comments=("#", "@"), unpack=True)

    # Plotting
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
        3, 1, figsize=(6, 8), sharex=False, sharey=False
    )
    fig.suptitle(f"NPT analysis {mean_box_length} $\\AA$")
    ax1.plot(npt_time, npt_potential, label="Potential energy")
    ax1.set_ylabel("Energy (kJ/mol)")
    ax1.legend()
    ax2.plot(npt_time, npt_temperature, label="Temperature")
    ax2.set_ylabel("Temperature (K)")
    ax2.legend()
    ax3.plot(npt_time, npt_pressure, label="Pressure")
    ax3.set_ylabel("Pressure (bar)")
    ax3.set_xlabel("Time (ps)")
    ax3.legend()
    for ax in (ax1, ax2, ax3):
        ax.grid(True, linestyle=":", linewidth=0.5)
    plt.tight_layout()
    plt.savefig(
        output_path,
        dpi=300,
    )

    return
