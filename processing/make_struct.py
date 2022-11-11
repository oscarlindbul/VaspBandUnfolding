import numpy as np
import ase.build,ase.io
import argparse

parser = argparse.ArgumentParser("Create supercell from basic cell")
parser.add_argument("structure", metavar="structure", help="Primitive cell of structure")
parser.add_argument("size", metavar="size", type=int, nargs=3, help="Scaling factors to be used for the structure")
parser.add_argument("-o", dest="output", type=str, default="POSCAR_out", help="Name of output file")
args = parser.parse_args()

prim_struct = ase.io.read(args.structure,format="vasp")
struct = prim_struct
print(args.size)
for i,repeat in enumerate(args.size):
    base_struct = struct
    for j in range(repeat-1):
        struct = ase.build.stack(struct, base_struct, axis=i)

struct = ase.build.sort(struct)
ase.io.write(args.output, struct, format="vasp")
