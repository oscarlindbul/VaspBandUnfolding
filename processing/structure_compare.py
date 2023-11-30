import argparse
import httk
import math

parser = argparse.ArgumentParser("Get a measure of relaxation difference between structures")
parser.add_argument("poscar1", metavar="poscar1", help="Position file 1")
parser.add_argument("poscar2", metavar="poscar2", help="Position file 2")
parser.add_argument("-max", dest="print_max", default=False, action="store_true", help="Prints the maximum dislocation")
parser.add_argument("-Q", dest="print_Q", default=False, action="store_true", help="Calculations the reduced coordinate perturbation between the structures")
args = parser.parse_args()

start_data = httk.iface.vasp_if.poscar_to_strs(args.poscar1, included_decimals=10)
end_data = httk.iface.vasp_if.poscar_to_strs(args.poscar2, included_decimals=10)
start_coords = httk.core.vectors.FracVector.create(start_data[3], simplify=True)
end_coords = httk.core.vectors.FracVector.create(end_data[3], simplify=True)
cell = httk.core.vectors.FracVector.create(start_data[0], simplify=True)
total_distance = 0
max_distance = 0
ipr = 0
Q = 0
element_list = [ Z for n,Z in zip(start_data[5],start_data[6]) for k in range(n) ]
if args.print_Q:
    from pymatgen.core.periodic_table import Element
    element_masses = [ Element.from_Z(Z).atomic_mass for Z in element_list ]
for i in range(len(start_coords)):
    distance_vector = (start_coords[i]-end_coords[i]).normalize_half()*cell
    distance = distance_vector.lengthsqr().to_float()
    if args.print_Q:
        Q += distance*element_masses[i]
    total_distance += distance
    ipr += distance**2
    if distance > max_distance: max_distance = distance
    if args.print_max:
        print("Max dislocation: {}", max_distance)
Q = math.sqrt(Q)
total_atoms = 0
for i in range(len(start_data[5])):
    total_atoms += start_data[5][i]
average_relaxation = total_distance/total_atoms
ipr = ipr/total_distance**2
print("Participation: {} A".format(1/ipr))
print("Position difference (RMS): {} A:".format(math.sqrt(total_distance)))
print("Position difference (abs max): {} A". format(math.sqrt(max_distance)))
if args.print_Q:
    print("Reduced coordinate deviation (Q): {} amu^0.5 A".format(Q)) 
