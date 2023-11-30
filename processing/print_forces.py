from argparse import ArgumentParser as Parser
import numpy as np
import xml.etree.ElementTree as ET

parser = Parser("Read vasp xml file to create xsf file with positions and forces")
parser.add_argument("xml_file", metavar="xml file", help="xml file to read")
parser.add_argument("-o", dest="output", default="forces.xsf", help="Name of output file")
args = parser.parse_args()

tree = ET.parse(args.xml_file)
root = tree.getroot()

cell_basis = np.zeros((3,3))
basis_root = root.find("./structure/crystal/varray[@name='basis']")
for i,v in enumerate(basis_root):
    cell_basis[i,:] = list(map(float,v.text.split()))

atom_root = root.find("./atominfo")
atom_N = int(atom_root.find("./atoms").text)
atom_type_roots = atom_root.findall("./array[@name='atoms']/set/rc/c")
atom_types = [atom.text.strip() for atom in atom_type_roots[::2]]
position_tags = root.findall("./structure[@name='finalpos']/varray[@name='positions']/v")
atom_pos = np.array([list(map(float,pos.text.split())) for pos in position_tags])
for i in range(atom_pos.shape[0]):
    atom_pos[i] = np.matmul(atom_pos[i,:],cell_basis)

force_tags = root.findall("./calculation/varray[@name='forces']/v")
for force in force_tags:
    print(force.text)
forces = np.array([list(map(float,force.text.split())) for force in force_tags])

with open(args.output, "w") as xsf:
    xsf.write("ATOMS\n")
    for i in range(atom_N):
        xsf.write("{}\t{}  {}  {}  {}  {}  {}\n".format(atom_types[i], *np.concatenate((atom_pos[i,:], forces[i,:]))))

