import MDAnalysis
import sys

modulename = "MDAnalysis"
if modulename not in sys.modules:
    print("You have not imported the {} module".format(modulename))
