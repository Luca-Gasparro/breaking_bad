import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import matplotlib.pyplot as plt
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utility import api_com_extractor

api_com_extractor(
    topology_file="dry_nvt.tpr",
    trajectroy_file="dry_nvt_trim.xtc",
    api_residue_name="NAP",
    output_filename="test",
    start_time=2000,
)
