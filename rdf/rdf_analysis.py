import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utility import api_com_extractor, load_api_coms, dummy_universe

api_com_extractor(
    topology_file="dry_nvt.tpr",
    trajectroy_file="dry_nvt_trim_whole.xtc",
    api_residue_name="NAP",
    output_filename="test",
    start_time=2000,
)
api_com_array, box_lengths = load_api_coms("test.npz")

fake_uni = dummy_universe(api_com_file="test.npz")
print(fake_uni.dimensions)

rdf = InterRDF(
    fake_uni.atoms, fake_uni.atoms, nbins=200, range=(0, 40), exclusion_block=(1, 1)
)
rdf.run()
plt.plot(rdf.results.bins, rdf.results.rdf)
plt.savefig("test.png", dpi=300)
