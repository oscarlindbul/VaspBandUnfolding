import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
import re, sys, os
from helperfuncs import open_file_with_check
import helperfuncs as helper
from calculate_charge_correction import get_LZ_correction, get_FNV_correction
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
make_parser.add_argument("-FNV", dest="FNV_pos", default=None, nargs=3, help="Sets the correction to be calculated using the FNV approach, provided relative position of defect")
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
    x = data[:, 0]
    vals = data[:, 1:-1]
    min_vals = data[:, -1]
    if args.plot_all:
        for d in range(vals.shape[1]):
            plt.plot(x, vals[:, d], label=headers[d+1])
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
    if sys.version_info[0] >= 3:
        potential_vals = np.load(file_name, allow_pickle=True).item()
    else:
        potential_vals = np.load(file_name).item()
else:
    potential_vals = {}
    for i in range(0, len(args.mu_pots), 2):
        species, val = args.mu_pots[i:i+2]
        potential_vals[species] = float(val)

# physical parameters
madelung = 1.64132  # 4H-SiC
epsilon = 9.66      # 4H-SiC
LZ_coeff = 0.65
eps_0 = constants.epsilon_0
e = constants.e

ref_energy = helper.get_final_energy(args.ref_outcar)
ref_ion_n = helper.get_ion_content(args.ref_outcar)
line_ind = 0

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
    bands, occups = helper.get_bands(args.ref_outcar)
    bandgap = [0, 0]
    for i in range(2):
        bands_s = bands[i, :, :]
        occups_s = occups[i, :, :]
        CBM = np.min(bands_s[occups_s == 0])
        VBM = np.max(bands_s[occups_s > 0])
        bandgap[i] = CBM - VBM
    bandgap = sum(bandgap) / 2

E_F = np.linspace(0, bandgap, 500)
form_E = np.zeros((len(E_F), len(args.outcars)))
system_names = []
for outcar_ind, outcar in enumerate(args.outcars):
    ion_n = helper.get_ion_content(outcar)
    energy = helper.get_final_energy(outcar)
    volume = helper.get_volume(outcar)
    charge = helper.get_charge_OUTCAR(outcar)   # effective charge
    effective_L = np.cbrt(volume)               # effective volume

    # calculate differences in atomic makeup
    component_diffs = {}
    all_elements = set(list(ion_n.keys()) + list(ref_ion_n.keys()))
    for element in all_elements:
        component_diffs[element] = ion_n[element] - ref_ion_n[element]
    # chemical potential energy differences
    pot_E = [potential_vals[element]*change for element, change in component_diffs.items()]
    pot_E_val = sum(pot_E)

    correction = 0
    if args.FNV_pos is not None:
        def_locpot = os.path.dirname(os.path.realpath(outcar)) + "/LOCPOT"
        ref_locpot = os.path.dirname(os.path.realpath(args.ref_outcar)) + "/LOCPOT"
        _, correction = get_FNV_correction(outcar, epsilon, def_locpot, ref_locpot, args.FNV_pos)
    else:
        # J to eV
        correction = get_LZ_correction(outcar, epsilon, madelung)
        #correction = LZ_coeff*(charge*e)**2*madelung/(2*4*np.pi*eps_0*epsilon*(effective_L*1e-10)) * (1/e)
    form_E[:, outcar_ind] = energy - ref_energy - pot_E_val + charge*(E_F + VBM) + correction

    if args.verbose:
        print("For system {}:".format(outcar))
        print("Charge = {}".format(charge))
        print("L = {}; V = {}".format(effective_L, volume))
        print("Energy = {}; Reference energy = {}".format(energy, ref_energy))
        print("Total potential energies = {}; Net result = {}".format(list(zip(pot_E, list(component_diffs.keys()))), -pot_E_val))
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

