import re
import numpy as np
from argparse import ArgumentParser as Parser

parser = Parser("Extract band info from BSEFATBANDS file")
parser.add_argument("fatbandfile", metavar="FATBANDS file")
args = parser.parse_args()

band_pattern = re.compile("\s*(\d+)\s*BSE eigenvalue\s*([\d.]+)\s*IP-eigenvalue:\s*([\d.]+)\s*")
with open(args.fatbandfile, "r") as reader:
    lines = reader.readlines()

eigenvals = []
e_h_pairs = []
c = 0
i = 0
while i < len(lines):
    line = lines[i]
    m = band_pattern.match(line)
    if m is not None:
        c += 1
        eh_pairs = []
        eigs = tuple(map(float,m.groups([1,2])))
        eigenvals.append(eigs)
        i += 1
        for j in range(i,len(lines)):
            vals = lines[j].split()
            if len(vals) < 7:
                i = j-1
                break
            else:
                kpoint = np.array(map(float,vals[:3]))
                E_h = float(vals[3])
                E_e = float(vals[4])
                coupling = float(vals[5])
                band_from = int(vals[6])
                band_to = int(vals[7])
                eh_pairs.append((kpoint, E_h, E_e, coupling, band_from, band_to))
        e_h_pairs.append(eh_pairs)
    i += 1

print("Saving...")
np.savez("fatband_data.npz", eigen=eigenvals, eh_pairs=e_h_pairs)

