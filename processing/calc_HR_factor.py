import numpy as np
import matplotlib.pyplot as plt
import re, sys, glob, os
from argparse import ArgumentParser as Parser
from pymatgen.core.periodic_table import Element
from scipy import constants
import httk

parser = Parser()
parser.add_argument("outcars", metavar="outcar files", type=str, nargs="+", help="Outcar files of static perturbation runs")
parser.add_argument("gref", metavar="ground reference", type=str, help="Outcar of reference calculation defined as Q=0")
parser.add_argument("exref", metavar="ex reference", type=str, help="Outcar of excited calculation defined as Q=Q_max")
parser.add_argument("-plot", dest="plot", action="store_true", default=False, help="Plots intermediate results for debugging and visualization")
parser.add_argument("-precision", dest="precision", default=6, type=int, help="Number of decimals to include in FracVector")
args = parser.parse_args()


name_pattern = re.compile("POSCAR = (.*)")
valenz_pattern = re.compile("Ionic Valenz")
zval_pattern = re.compile("ZVAL\s+=\s+(.*)")
mass_pattern = re.compile("Mass of Ions in am")
ion_type_n_pattern = re.compile("ions per type \=\s+(.*)")
ion_type_pattern = re.compile("VRHFIN\s*=\s*([A-Za-z]+)")
energy_pattern = re.compile("energy without entropy\s*=\s*([-e.\d]+).*")
#electron_pattern = re.compile("NELECT\s*=\s*([\d.]+)")
#volume_pattern = re.compile("volume of cell\s+:\s+([\d.]+)")
#spin_comp_pattern = re.compile("spin component (\d)")
#kpoint_pattern = re.compile("k-point\s+(\d+)")
#band_pattern = re.compile("(\d+)\s+([-.\d]+)\s+([\d.])")

cell_vec_begin_pattern = re.compile("direct lattice vectors")
cell_pos_begin_pattern = re.compile("position of ions in fractional coordinates")
three_val_pattern = re.compile("(-?[e\d.]+)\s*(-?[e\d.]+)\s*(-?[e\d.]+)")

class OutcarInfo:
    def __init__(self, energy, ion_count, ion_masses, ion_types, lattice_vectors, lattice_pos, lattice_frac, species):
        self.energy = energy
        self.ion_n = ion_count
        self.ion_m = ion_masses
        self.ion_types = ion_types
        self.Avecs = lattice_vectors
        self.positions = lattice_pos
        self.frac_positions = lattice_frac
        self.species_of_ions = species


def read_outcar(outcar):
    with open(outcar, "r") as outcar_file:
        lines = outcar_file.readlines()

    ion_types = []
    ion_n = {}
    ion_m = {}
    A_vecs = np.zeros((3,3))
    lattice_frac = []
    lattice_element = []
    energy = 0
    line_ind = 0
    while line_ind < len(lines):
        # check ion types
        m = ion_type_pattern.search(lines[line_ind])
        if m is not None:
            element = m.group(1)
            ion_types.append(element)
            if element not in ion_n:
                ion_n[element] = 0

        m = ion_type_n_pattern.search(lines[line_ind])
        if m is not None:
            element_n = list(map(int, m.group(1).split()))
            tot = 0
            for i,N in enumerate(element_n):
                ion_n[ion_types[i]] += N
                tot += N
            lattice_element = [0]*tot
            i = 0
            for k, element in enumerate(ion_types):
                for n in range(element_n[k]):
                    lattice_element[i] = element
                    i += 1
            lattice_element = np.array(lattice_element)

        # check mass of ions
        m = mass_pattern.search(lines[line_ind])
        if m is not None:
            line_ind += 1
            masses = list(map(float, lines[line_ind].split()[2:]))
            for i,element in enumerate(ion_types):
                ion_m[element] = masses[i]
        
        # check energy
        m = energy_pattern.search(lines[line_ind])
        if m is not None:
            energy = float(m.group(1))

        # read lattice
        m = cell_vec_begin_pattern.search(lines[line_ind])
        if m is not None:
            for i in range(3):
                line_ind += 1
                # OBS: each row is one vector
                vals = three_val_pattern.search(lines[line_ind])
                A_vecs[i,:] = np.array(list(map(float, vals.group(1,2,3))))

        m = cell_pos_begin_pattern.search(lines[line_ind])
        if m is not None:
            line_ind += 1
            vals = three_val_pattern.search(lines[line_ind])
            vals = vals.group(1,2,3)
            while len(vals) == 3:
                lattice_frac.append(list(map(float,vals)))
                line_ind += 1
                vals = lines[line_ind].split()
            lattice_frac = np.array(lattice_frac) 

        line_ind += 1
    ion_types = list(ion_n.keys())
    lattice = np.matmul(lattice_frac, A_vecs)

    return OutcarInfo(energy, ion_n, ion_m, ion_types, A_vecs, lattice, lattice_frac, lattice_element)

g_outcar = read_outcar(args.gref)
g_contcar = glob.glob(os.path.dirname(args.gref) + "/CONTCAR*")[0]
ex_outcar = read_outcar(args.exref)
ex_contcar = glob.glob(os.path.dirname(args.exref) + "/CONTCAR*")[0]

g_httk_data = httk.iface.vasp_if.poscar_to_strs(g_contcar, included_decimals=args.precision)
g_coords = httk.core.vectors.FracVector.create(g_httk_data[3], simplify=True)
cell = httk.core.vectors.FracVector.create(g_httk_data[0], simplify=True)
ex_httk_data = httk.iface.vasp_if.poscar_to_strs(ex_contcar, included_decimals=args.precision)
ex_coords = httk.core.vectors.FracVector.create(ex_httk_data[3], simplify=True)
element_list = [ Z for n,Z in zip(g_httk_data[5],g_httk_data[6]) for k in range(n) ]
element_masses = [ Element.from_Z(Z).atomic_mass for Z in element_list ]

ref_Q = 0
for i in range(len(g_coords)):
    distance_vector = (ex_coords[i]-g_coords[i]).normalize_half()*cell
    distance = distance_vector.lengthsqr().to_float()
    ref_Q += distance*element_masses[i]
ref_Q = np.sqrt(ref_Q)


assert np.all(np.array([ t1 == t2 for t1,t2 in zip(g_outcar.ion_types,ex_outcar.ion_types)]))
assert np.all(np.array([ g_outcar.ion_n[e] == ex_outcar.ion_n[e] for e in g_outcar.ion_types]))

run_infos = [read_outcar(outcar) for outcar in args.outcars]
contcars = [glob.glob(os.path.dirname(outcar) + "/CONTCAR*")[0] for outcar in args.outcars]
run_Qs = [0]
run_Ediff = [0]
for i, info in enumerate(run_infos):
    try:
        print(g_outcar.ion_n.items(), info.ion_n.items())
        assert np.all(np.array([ t1 == t2 for t1,t2 in zip(g_outcar.ion_types,info.ion_types)]))
        assert np.all(np.array([ g_outcar.ion_n[e] == info.ion_n[e] for e in g_outcar.ion_types]))
    except AssertionError as e:
        print(repr(e))
        print("The structures may be incompatible, proceed anyway?")
        answer = "y"
        #if sys.version_info[0] < 3:
        #    answer = raw_input("y/n?")
        #else:
        #    answer = input("y/n?")
        if "y" not in answer: 
            exit()
    httk_data = httk.iface.vasp_if.poscar_to_strs(contcars[i], included_decimals=args.precision)
    httk_coords = httk.core.vectors.FracVector.create(httk_data[3], simplify=True)
    Q = 0
    for i in range(len(httk_coords)):
        distance_vector = (httk_coords[i]-g_coords[i]).normalize_half()*cell
        distance = distance_vector.lengthsqr().to_float()
        Q += distance*element_masses[i]
    Q = np.sqrt(Q)
    Ediff = info.energy - g_outcar.energy
    run_Qs.append(Q)
    run_Ediff.append(Ediff)

# should be a 0.5 * omega**2 * Q^2 function
deg = 2
fitting = np.polyfit(run_Qs, run_Ediff, deg=deg)
if args.plot:
    plt.scatter(run_Qs, run_Ediff)
    poly = lambda x,deg: sum([fitting[i]*x**(deg-i) for i in range(0,deg+1)])
    x = np.linspace(-ref_Q, ref_Q, 300)
    vals = poly(x, deg)
    plt.plot(x, vals)
    plt.xlabel("Q (amu^-0.5 A)")
    plt.ylabel("E (eV)")
    plt.show()

hbar_eV = constants.value("reduced Planck constant in eV s")
unit_Q_to_kgM2 = np.sqrt(1e-20*constants.value("atomic mass constant"))
omega = np.sqrt(2*fitting[0])
omega_eV = omega * hbar_eV * np.sqrt(constants.e / unit_Q_to_kgM2**2) # conversion from (eV/(amu * A^2))^0.5 to eV via hbar
omega_Hz = omega_eV / hbar_eV

HR = omega_Hz*(ref_Q*unit_Q_to_kgM2)**2/(2*constants.hbar)
print("Total reduced coordinate (Q): {} amu^0.5 A".format(ref_Q))
print("One-phonon frequency (w): {} eV".format(omega_eV))
print("Huang-Rhys factor (S): {}".format(HR))
print("Quantum efficiency: {}%".format(np.exp(-HR)*100))
