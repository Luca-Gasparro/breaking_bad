# Module contaning code to analyse quantities obtained via MD simulations
import numpy as np
import matplotlib.pyplot as plt


def energy_minimisation_pe(polymer_xvg, api_polymer_xvg, wet_api_polymer_xvg, box_size):
    """Produces three subplots which display how the potential energy of the system
    evolves via the steepest descent method during the energy minimisation stage."""

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
    fig.suptitle(f"Energy minimisation {box_size} $\AA$")
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
        f"/storage/chem/phuqdw/breaking_bad/sim_analysis/energy_min/energy_minimisation_{box_size}_ang.png",
        dpi=300,
    )
    return
