import numpy as np
import ase.io, ase.build
from ase import Atoms
import argparse

parser = argparse.ArgumentParser("Create supercell from basic cell")
parser.add_argument("structure", metavar="structure", help="Primitive cell of structure")
parser.add_argument("defect_info", metavar="defect instructions", help="Defect implantation instructions (removal/addition, species, position)")
parser.add_argument("-o", dest="output", type=str, default="POSCAR_out", help="Name of output file")
args = parser.parse_args()

struct = ase.io.read(args.structure, format="vasp")

cell = struct.get_cell()
#read defect instructions
removal = []
new_atoms = None
with open(args.defect_info, "r") as reader:
    for line in reader.readlines():
        infos = line.split()
        instruction = infos[0]
        species = infos[1]
        pos = list(map(float,infos[2:5]))
        if instruction == "-":
            removal.append((species, pos))
        elif instruction == "+":
            new_atom = Atoms(infos[1], positions=[pos], cell=cell, pbc=[1, 1, 1])
            if new_atoms is None:
                new_atoms = new_atom
            else:
                new_atoms += new_atom

for vacancy in removal:
    del struct[[atom.index for atom in struct if atom.symbol == vacancy[0] and np.linalg.norm(np.array(atom.position) - np.array(vacancy[1])) < 0.5]]

defect_struct = ase.build.sort(struct + new_atoms)

ase.io.write(args.output, defect_struct, format="vasp")


