import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
import re, sys, os
from helperfuncs import open_file_with_check
from argparse import ArgumentParser as Parser

parser = Parser("Script for calculating formation energies of give input files")
sub_parsers = parser.add_subparsers(help="Either create new formation energy file or combine several existing ones", dest="command")

make_parser = sub_parsers.add_parser("make", help="Create formation energy file")
make_parser.add_argument("outcars", metavar="OUTCARs", nargs="+", type=str, help="List of outcar files for runs")
make_parser.add_argument("ref_outcar", metavar="Reference OUTCAR", type=str, help="Unmodified neutral bulk system")
gap_manual_group = make_parser.add_argument_group("Manual setting")
gap_auto_group = make_parser.add_argument_group("Semi-automatic setting")
gap_manual_group.add_argument("-bandgap", dest="bandgap", type=float, default=None, help="Value of bandgap to assume for formation energy range")
gap_manual_group.add_argument("-vbm", dest="vbm", type=float, default=None, help="VBM of the system")
gap_auto_group.add_argument("-eigfile", dest="eigfile", type=str, default=None, help="EIGENVAL file to use for bandgap estimation (such as that of reference system)")
make_parser.add_argument("-o", dest="output", type=str, default="form_E.txt", help="Name of output data file")
make_parser.add_argument("-plot", dest="plot", default=False, action="store_true", help="Should we plot the result?")
make_parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="Extra output for debugging purposes")
make_parser.add_argument("--only-all", dest="only_all", default=False, action="store_true", help="Setting to skip the calculation of the minimum and (if plotting) instead plot all states")
make_parser.add_argument("-label", dest="label", default="None", help="Header for minimum value column")
make_parser.add_argument("-FNV", dest="FNV_pos", default=None, nargs=3, type=float, help="Sets the correction to be calculated using the FNV approach, provided relative position of defect")
make_parser.add_argument("-mu", dest="mu_pots", nargs="+", default=None, help="List of chemical potentials (species, value)")

combine_parser = sub_parsers.add_parser("combine", help="Combine formation energy files into new one")
combine_parser.add_argument("energy_files", metavar="Energy file list", nargs="+", type=str, help="List of energy files to combine")
combine_parser.add_argument("-o", dest="output", default="form_E_combined.txt", help="Name of combined file")
combine_parser.add_argument("-plot", dest="plot", default=False, action="store_true", help="Should we plot the result?")
combine_parser.add_argument("--only-all", dest="only_all", default=False, action="store_true", help="Setting to skip the calculation of the minimum and (if plotting) instead plot all states")

plot_parser = sub_parsers.add_parser("plot", help="Plot the given energy file")
plot_parser.add_argument("file", metavar="file", type=str, help="File to plot the energy of")
plot_parser.add_argument("--all", dest="plot_all", action="store_true", default=False, help="Plot all containing formation energies")

args = parser.parse_args()

if args.command == "plot":
    data = np.loadtxt(args.file)
    with open(args.file, "r") as reader:
        headers = reader.readline().split()[1:]
    x = data[:,0]
    vals = data[:,1:-1]
    min_vals = data[:,-1]
    if args.plot_all:
        for d in range(vals.shape[1]):
            plt.plot(x, vals[:,d], label=headers[d+1])
    else:
        plt.plot(x, min_vals, label=headers[-1])
    plt.legend()
    plt.xlabel("Ef - VBM")
    plt.ylabel("Energy (eV)")
    plt.show()
    exit()
    
elif args.command == "combine":
    data_files = []
    headers = []
    for i in range(len(args.energy_files)):
        data_files.append(np.loadtxt(args.energy_files[i]))
        with open(args.energy_files[i], "r") as reader:
            headers.append(reader.readline().split()[-1])
    # combine data points assuming linear behavior between points
    x_vals = [ data[:,0] for data in data_files]
    main_vals = [ data[:,-1] for data in data_files]
    x_inds = [0]*len(x_vals)
    new_x = [0]*sum([len(vals) for vals in x_vals])
    new_vals = [[None for t_len in range(len(new_x))] for d in main_vals ]
    new_t_ind = 0
    # expand time values (include all and reduce later)
    while np.any(np.array(x_inds) < np.array([len(d) for d in x_vals])):
        min_t = np.inf
        t_ind = 0
        for i in range(len(x_inds)):
            if x_inds[i] < len(x_vals[i]):
                if x_vals[i][x_inds[i]] <= min_t:
                    min_t = x_vals[i][x_inds[i]]
                    t_ind = i
        new_x[new_t_ind] = min_t
        new_vals[t_ind][new_t_ind] = main_vals[t_ind][x_inds[t_ind]]
        x_inds[t_ind] += 1
        new_t_ind += 1
    
    # expand energy values with linear interpolation
    for d_ind in range(len(main_vals)):
        d_list = new_vals[d_ind]
        # find initial values for interpolation
        start_ind = -1
        end_ind = -1
        for i in range(len(new_x)): 
            if d_list[i] is not None:
                if start_ind < 0:
                    start_ind = i
                elif end_ind < 0:
                    end_ind = i
                    break
        interpol = lambda x,ind2,ind1: (d_list[ind2] - d_list[ind1])/(new_x[ind2]-new_x[ind1]) * x + (d_list[ind2] - (d_list[ind2] - d_list[ind1])/(new_x[ind2] - new_x[ind1])*new_x[ind2])
        # fill value lists with interpolation
        for i in range(len(new_x)):
            if d_list[i] is None:
                d_list[i] = interpol(new_x[i], end_ind, start_ind)
            else:
                if i != end_ind:
                    tmp = end_ind
                    end_ind = i
                    start_ind = tmp
    # remove duplicates
    start_ind = 0
    final_x = []
    final_vals = [[] for d in main_vals]
    for i in range(len(new_x)):
        diff = abs(new_x[i] - new_x[start_ind])
        if diff < 1e-5:
            handle_dup = True
            continue
        if diff > 1e-5 or i == len(new_x)-1:
            if handle_dup:
                handle_dup = False
                x = np.mean(np.array(new_x[start_ind:i+1]))
                final_x.append(x)
                for d in range(len(final_vals)):
                    avg_val = np.mean(np.array(new_vals[d][start_ind:i+1]))
                    final_vals[d].append(new_vals[d][i])
                start_ind = i
            else:
                final_x.append(new_x[i])
                for d in range(len(final_vals)):
                    final_vals[d].append(new_vals[d][i])
    final_vals = np.array(final_vals).T
    final_x = np.array(final_x).reshape((len(final_x), 1))

    if not args.only_all:
        min_vals = np.min(final_vals, axis=1).reshape((final_vals.shape[0], 1))

    if args.plot:
        if args.only_all:
            for i in range(final_vals.shape[1]):
                plt.plot(final_x, final_vals[:,i], label=headers[i])
        else:
            plt.plot(final_x, min_vals, label="Hull")
        plt.show()
    
    if args.only_all:
        save_data = np.concatenate([final_x, final_vals], axis=1)
    else:
        save_data = np.concatenate([final_x, final_vals, min_vals], axis=1)
    np.savetxt(args.output, save_data, header="Ef " + " ".join(headers) + " Total")

    exit()

# otherwise "command" must be "make"

my_path = os.path.abspath(os.path.dirname(__file__))
file_name = os.path.join(my_path, "pot_dict.save.npy")

if args.mu_pots is None:
    if sys.version_info[0] < 3:
        potential_vals = np.load(file_name, allow_pickle=True).item()
    else:
        potential_vals = np.load(file_name).item()
else:
    potential_vals = {}
    for i in range(0, len(args.mu_pots), 2):
        species, val = args.mu_pots[i:i+2]
        potential_vals[species] = float(val)

# outcar patterns
name_pattern = re.compile("POSCAR = (.*)")
valenz_pattern = re.compile("Ionic Valenz")
zval_pattern = re.compile("ZVAL\s+=\s+(.*)")
ion_type_n_pattern = re.compile("ions per type \=\s+(.*)")
ion_type_pattern = re.compile("VRHFIN\s*=\s*([A-Za-z]+)")
energy_pattern = re.compile("energy without entropy\s*=\s*([-e.\d]+).*")
electron_pattern = re.compile("NELECT\s*=\s*([\d.]+)")
volume_pattern = re.compile("volume of cell\s+:\s+([\d.]+)")
spin_comp_pattern = re.compile("spin component (\d)")
kpoint_pattern = re.compile("k-point\s+(\d+)")
band_pattern = re.compile("(\d+)\s+([-.\d]+)\s+([\d.])")

# physical parameters
madelung = 1.64132 # 4H-SiC
epsilon = 9.66 # 4H-SiC
LZ_coeff = 0.65
eps_0 = constants.epsilon_0
e = constants.e

ref_energy = 0
ref_ion_n = {}
ref_ion_types = []
with open_file_with_check(args.ref_outcar) as outcar_file:
    lines = outcar_file.readlines()
line_ind = 0
while line_ind < len(lines): 
    # check ion types
    m = ion_type_pattern.search(lines[line_ind])
    if m is not None:
        element = m.group(1)
        ref_ion_types.append(element)
        if element not in ref_ion_n:
            ref_ion_n[element] = 0
    m = ion_type_n_pattern.search(lines[line_ind])
    if m is not None:
        element_n = list(map(int, m.group(1).split()))
        for i,N in enumerate(element_n):
            ref_ion_n[ref_ion_types[i]] += N
    
    # check energy
    m = energy_pattern.search(lines[line_ind])
    if m is not None:
        ref_energy = float(m.group(1))

    line_ind += 1

# bandgap estimation
bandgap = 0
VBM = 0
if args.bandgap is not None and args.vbm is not None:
    VBM = args.bandgap
    bandgap = args.bandgap
elif args.eigfile is not None:
    pass
else:
    import warnings
    warnings.warn("No info given for bandgap, will look in reference outcar")

    bandgap_info = [[],[]]
    with open_file_with_check(args.ref_outcar) as outcar_file:
        lines = outcar_file.readlines()
    line_ind = 0
    spin_comp = None
    search_kpoint = False
    while line_ind < len(lines):
        m = spin_comp_pattern.search(lines[line_ind])
        if m is not None:
            spin_comp = int(m.group(1))
            search_kpoint = True
        if search_kpoint:
            m = kpoint_pattern.search(lines[line_ind])
            if m is not None:
                kpoint = int(m.group(1))
                if len(bandgap_info[spin_comp-1]) < kpoint:
                    bandgap_info[spin_comp-1].append(0)
                line_ind += 2
                occ_E = 0
                while True:
                    m = band_pattern.search(lines[line_ind])
                    if m is None:
                        search_kpoint = False
                        break
                    band_id,band_E,band_occ = m.group(1,2,3)
                    if float(band_occ) > 0:
                        occ_E = float(band_E)
                    else:
                        bandgap_info[spin_comp-1][kpoint-1] = (float(band_E) - occ_E, occ_E)
                        search_kpoint = False
                        break
                    line_ind += 1

        line_ind += 1
    n = sum([sum([ 1 for i in range(len(bandgap_info[j]))]) for j in range(len(bandgap_info))])
    bandgap = sum([sum([gap for gap,vbm in l]) for l in bandgap_info]) / n
    VBM = sum([sum([vbm for gap,vbm, in l]) for l in bandgap_info]) / n

    
E_F = np.linspace(0, bandgap, 500)
form_E = np.zeros((len(E_F), len(args.outcars)))
system_names = []
for outcar_ind,outcar in enumerate(args.outcars):
    print(outcar)
    with open_file_with_check(outcar) as outcar_file:
        lines = outcar_file.readlines()
    line_ind = 0
    zvals = {}
    ion_types = []
    ion_n = {}
    energy = 0
    volume = 0
    while line_ind < len(lines):
        # check name
        m = name_pattern.search(lines[line_ind])
        if m is not None:
            system_names.append(m.group(1))

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
            for i,N in enumerate(element_n):
                ion_n[ion_types[i]] += N
        
        # check valences
        m = valenz_pattern.search(lines[line_ind])
        if m is not None:
            line_ind += 1
            m = zval_pattern.search(lines[line_ind])
            zvals_n = list(map(int,map(float,m.group(1).split())))
            for i,number in enumerate(zvals_n):
                element = ion_types[i]
                zvals[element] = zvals_n[i]
        
        # check energy
        m = energy_pattern.search(lines[line_ind])
        if m is not None:
            energy = float(m.group(1))

        # check electrons
        m = electron_pattern.search(lines[line_ind])
        if m is not None:
            outcar_nelect = int(float(m.group(1)))

        # check volume
        m = volume_pattern.search(lines[line_ind])
        if m is not None:
            volume = float(m.group(1))

        line_ind += 1

    nelect = np.sum([ion_n[element]*zvals[element] for element in ion_n.keys()])
    charge = nelect - outcar_nelect # effective charge
    effective_L = np.cbrt(volume) # effective volume
    
    # calculate differences in atomic makeup
    component_diffs = {}
    all_elements = set(ion_n.keys() + ref_ion_n.keys())
    for element in all_elements:
        now = ion_n[element] if element in ion_n else 0
        before = ref_ion_n[element] if element in ref_ion_n else 0
        component_diffs[element] = now - before
    pot_E = [potential_vals[element]*change for element,change in component_diffs.items()] # chemical potential energy differences
    pot_E_val = sum(pot_E)

    correction = 0
    if args.FNV_pos is not None:
        pass    
    else:
        correction = LZ_coeff*(charge*e)**2*madelung/(2*4*np.pi*eps_0*epsilon*(effective_L*1e-10)) * (1/e) # J to eV
    form_E[:,outcar_ind] = energy - ref_energy - pot_E_val + charge*(E_F + VBM) + correction
    
    if args.verbose:
        print("For system {}:".format(system_names[outcar_ind]))
        print("Charge = {}".format(charge))
        print("L = {}; V = {}".format(effective_L, volume))
        print("Energy = {}; Reference energy = {}".format(energy, ref_energy))
        print("Total potential energies = {}; Net result = {}".format(zip(pot_E, list(component_diffs.keys())), -pot_E_val))
        print("Estimated bandgap = {}; with VBM = {}".format(bandgap, VBM))
        print("correction = {}".format(correction))

# calculate "hull"
min_form = np.min(form_E, axis=1)

if args.plot:
    if args.only_all:
        for i in range(form_E.shape[1]):
            plt.plot(E_F, form_E[:,i], label=system_names[i])
        plt.legend()
    else:
        plt.plot(E_F, min_form, label=args.label)
    plt.show()

if not args.only_all:
    save_data = np.concatenate([E_F.reshape((E_F.shape[0], 1)), form_E, min_form.reshape((E_F.shape[0], 1))], axis=1)
else:
    save_data = np.concatenate([E_F.reshape(form_E.shape), form_E], axis=1)
np.savetxt(args.output, save_data, header="Ef " + "None "*(save_data.shape[1]-2) + args.label)

