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
api_com_array = load_api_coms("test.npz")

fake_uni = dummy_universe(api_com_array=api_com_array)
