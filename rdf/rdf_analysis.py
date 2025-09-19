import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utility import api_simulation_extractor, dummy_universe, polymer_atom_extraction


def api_api_rdf(
    topology_file, trajectory_file, api_residue_name, output_filename, start_time
):
    """Computes and plots the RDF for API-API. Uses the centres of mass of the API
    as reference particles."""

    api_simulation_extractor(
        topology_file=topology_file,
        trajectroy_file=trajectory_file,
        api_residue_name=api_residue_name,
        output_filename=output_filename,
        start_time=start_time,
    )
    api_universe = dummy_universe(api_simulation_file=output_filename)

    api_rdf = InterRDF(
        api_universe.atoms,
        api_universe.atoms,
        nbins=75,
        range=(0, 30),
        exclusion_block=(1, 1),
    )
    api_rdf.run()

    return api_rdf


nitrogen_polymer_coords = polymer_atom_extraction(
    "dry_nvt_polymer.tpr", "dry_nvt_trim_polymers.xtc", "N1", start_time=2000
)
print(nitrogen_polymer_coords)
