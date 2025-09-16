import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt


def rdf_api(topology, trajectory, api_residue_name):
    """Computes the average radial distribution of the API using MDAnalysis"""

    # Select the API atoms
    u = mda.Universe(topology, trajectory)
    api = u.select_atoms(f"resname {api_residue_name}")

    # Compute the average radial distriution function between API molecules and plot
    rdf_api = InterRDF(api, api, nbins=50, range=(0.0, 10.0))
    rdf_api.run()
    plt.plot(rdf_api.results.bins, rdf_api.results.rdf)
    plt.savefig("api_rdf.png", dpi=300)

    return


rdf_api("dry_nvt.tpr", "dry_nvt_trim.xtc", "NAP")
