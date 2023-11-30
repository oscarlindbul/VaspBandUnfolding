#!/usr/bin/python
import numpy as np
import ase.build,ase.io
from ase import Atoms
import argparse

parser = argparse.ArgumentParser("Transfer (applicable) positions of cell 1 to cell 2")
parser.add_argument("cell1", metavar="From_Cell")
parser.add_argument("cell2", metavar="To_Cell")
parser.add_argument("-reffrom", dest="ref_1", nargs=3, type=float, default=None, help="Origin for cell 1")
parser.add_argument("-refto", dest="ref_2", nargs=3, type=float, default=None, help="Origin for cell 2")
parser.add_argument("-rel", dest="relative", action="store_true", default=False, help="Input positions in relative coordinates")
parser.add_argument("-o", dest="output", type=str, default="POSCAR_out", help="Name of output file")
parser.add_argument("-neg", dest="allow_neg", default=False, action="store_true", help="Allow negative values in resulting POSCAR")
parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="Add more printout statements")
args = parser.parse_args()


from_cell = ase.io.read(args.cell1, format="vasp")
to_cell = ase.io.read(args.cell2, format="vasp")

vaccuum_shift = np.array([0.2,0.2,0.2])
for cell,arg in ((from_cell,args.ref_1),(to_cell, args.ref_2)):
    if arg is not None:
        if args.relative:
            pos = cell.cell.cartesian_positions(arg)
        else:
            pos = np.array(arg)
        middle = cell.cell.cartesian_positions((0.5,0.5,0.5))
        cell.translate(vaccuum_shift - pos)
        cell.wrap(center=(0,0,0))

ase.io.write("to_cell.vasp", to_cell, format="vasp") 
ase.io.write("from_cell.vasp", from_cell, format="vasp") 

# removes atoms outside of to_cell from from_cell (if from_cell is bigger)
defects_outofbounds = set()
for i,from_atom in enumerate(from_cell):
    scaled_pos = to_cell.cell.scaled_positions(from_atom.position)
    if np.any(scaled_pos < -0.500) or np.any(scaled_pos > 0.500):
        defects_outofbounds.add(i)
if args.verbose:
    print("{} atoms ignored from from_cell".format(len(defects_outofbounds)))
for j in reversed(sorted(list(defects_outofbounds))):
    from_cell.pop(j)

smaller_cell = from_cell if from_cell.get_volume() < to_cell.get_volume() else to_cell

defect_vacancies = set()
defect_addons = set()
tolerance = 0.5 # A
matched = [False for i in range(len(from_cell))]
for j,to_atom in enumerate(to_cell):
    if args.verbose:
        print("Scanning atom {}/{}, {}/{} matched".format(j+1,to_cell.positions.shape[0], sum(matched), len(from_cell)), end="\r")
    scaled_pos = from_cell.cell.scaled_positions(to_atom.position)

    # ignore atom outside of from_cell (if to_cell is bigger)
    if np.any(scaled_pos < -0.501) or np.any(scaled_pos > 0.501):
        if args.verbose:
            print("Ignoring {} at {}                  ".format(to_atom.symbol, scaled_pos))
        continue
    min_mic_distance = [[np.inf,np.inf,np.inf], np.inf]
    min_ind = 0 
    for i,from_atom in enumerate(from_cell):
        if from_atom.symbol != to_atom.symbol or matched[i]:
            continue
        #scaled_pos = to_cell.cell.scaled_positions(from_atom.position)
        #if np.any(-0.55 >= scaled_pos) or np.any(scaled_pos >= 0.55):
        #    matched[i] = True
        #    continue
        mic_distance = ase.geometry.get_distances(to_atom.position, from_atom.position, cell=smaller_cell.cell, pbc=(True,True,True))
        #if mic_distance[1] < tolerance:
        #    matched[i] = True
        if min_mic_distance[1] > mic_distance[1]:
            min_mic_distance = mic_distance
            min_ind = i
    if min_mic_distance[1] < tolerance:
        if matched[min_ind]: # double match due to periodicity
            print("Error???")
            exit()
        matched[min_ind] = True
        to_atom.position = from_cell[min_ind].position
        if args.verbose:
            print("Found atom {}/{}, {}/{} matched               ".format(j+1,to_cell.positions.shape[0], sum(matched),len(from_cell)))
    else:
        if args.verbose:
            scaled_pos = from_cell.cell.scaled_positions(to_atom.position)
            print("Removed {} at {} due to dist {}".format(to_atom.symbol, scaled_pos, min_mic_distance[1]))
        defect_vacancies.add(j)

for i,from_atom in enumerate(from_cell):
    if not matched[i]:
        position = to_cell.cell.scaled_positions(from_atom.position)
        if np.all(-0.51 <= position) and np.all(position <= 0.51):
            if args.verbose:
                print("Added atom {} at {}".format(from_atom.symbol, position))
            defect_addons.add(from_atom)

for j in reversed(sorted(list(defect_vacancies))):
    to_cell.pop(j)
for atom in defect_addons:
    to_cell.append(atom)


if args.ref_2 is not None:
    to_cell.translate(pos - vaccuum_shift)
to_cell.wrap()
to_cell.center()
to_cell = to_cell[to_cell.numbers.argsort()]

ase.io.write(args.output, to_cell, format="vasp") 
exit()

for cell in [from_cell, to_cell]:
    for atom in cell:
        rel_pos = np.linalg.solve(cell.get_cell()[:], atom.position)
        for i,p in enumerate(rel_pos):
            if p > 0.5:
                rel_pos[i] -= 1
        atom.position = np.matmul(cell.get_cell()[:], rel_pos)

for atom2 in to_cell:
    for atom1 in from_cell:
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
        

