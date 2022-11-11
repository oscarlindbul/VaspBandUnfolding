import yaml
from yaml import CLoader
import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser("Calculates the ipr of phonon modes in data file")
parser.add_argument("data_file", metavar="data_file", help="File containing phonon modes (yaml file)")
args = parser.parse_args()

with open(args.data_file, "r") as f:
	data = yaml.load(f, Loader=CLoader)

phonon_data = data["phonon"][0]
print(phonon_data["band"][0].keys())
nlattice = len(data["points"])
nphonon = len(phonon_data["band"])
mode_freqs = np.array([d["frequency"] for d in phonon_data["band"]])
mode_eigenvectors = [0]*nphonon
for i in range(nphonon):
	mode_eigenvectors[i] = np.reshape(np.array([np.array(d)[:,0] for d in phonon_data["band"][i]["eigenvector"]]), (nlattice, 3))
lattice_points = np.array([d["coordinates"] for d in data["points"]])
lattice_vecs = [np.array(d) for d in data["lattice"]]

np.savez("phonon_data.npz", atoms=lattice_points, freqs = mode_freqs, eigs = mode_eigenvectors, lattice = lattice_vecs)
