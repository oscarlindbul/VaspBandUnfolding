#!/usr/bin/python
import numpy as np
import ase.build,ase.io
from ase import Atoms
import argparse

parser = argparse.ArgumentParser("Transfer (applicable) positions of cell 1 to cell 2")
parser.add_argument("cell1", metavar="Cell 1")
parser.add_argument("cell2", metavar="Cell 2")
parser.add_argument("-o", dest="output", type=str, default="POSCAR_out", help="Name of output file")
parser.add_argument("-neg", dest="allow_neg", default=False, action="store_true", help="Allow negative values in resulting POSCAR")
args = parser.parse_args()

from_cell = ase.io.read(args.cell1, format="vasp")
to_cell = ase.io.read(args.cell2, format="vasp")

for cell in [from_cell, to_cell]:
    for atom in cell:
        rel_pos = np.linalg.solve(cell.get_cell()[:], atom.position)
        for i,p in enumerate(rel_pos):
            if p > 0.5:
                rel_pos[i] -= 1
        atom.position = np.matmul(cell.get_cell()[:], rel_pos)

for atom1 in from_cell:
    for atom2 in to_cell:
        diff = 0
        is_pair = True
        for i in range(3):
            diff += (atom1.position[i] - atom2.position[i])**2
            if np.sqrt(diff) > 0.8:
                is_pair = False
                break
        if is_pair:
            atom2.position = atom1.position
            break
if not args.allow_neg:
    for cell in [from_cell, to_cell]:
        for atom in cell:
            rel_pos = np.linalg.solve(cell.get_cell()[:], atom.position)
            for i,p in enumerate(rel_pos):
                if p < 0:
                    rel_pos[i] += 1
            atom.position = np.matmul(cell.get_cell()[:], rel_pos)

ase.io.write(args.output, to_cell, format="vasp") 
        

