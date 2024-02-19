import numpy as np
import matplotlib.pyplot as plt
import os, re
from scipy import constants
import helperfuncs as helper
import ase.io
from argparse import ArgumentParser as Parser
import subprocess

LZ_coeff = 0.65
eps_0 = constants.epsilon_0
e = constants.e
ev2ry = e/(constants.h*constants.c*constants.Rydberg)

def get_LZ_correction(outcar, epsilon, madelung):
    volume = helper.get_volume(outcar)
    effective_L = np.cbrt(volume)
    charge = helper.get_charge_OUTCAR(outcar)

    correction = LZ_coeff*(charge*e)**2*madelung/(2*4*np.pi*eps_0*epsilon*(effective_L*1e-10)) * (1/e)
    return correction


def get_FNV_correction(outcar, epsilon, pot_def, pot_ref, defect_pos, averaging=3, verbose=True):
    script_path = os.path.realpath(os.path.dirname(__file__))
    pot_def = os.path.realpath(pot_def)
    pot_ref = os.path.realpath(pot_ref)
    cwd = os.getcwd()
    cmd_list = ["{}/sxdefectalign".format(script_path)]

    encut = helper.get_encut(outcar)
    charge = helper.get_charge_OUTCAR(outcar)
    if charge == 0:
        return "Defect correction (eV): 0 (incl. screening & alignment)", 0

    os.chdir(os.path.dirname(pot_def))
    
    cmd_list.append("--ecut {}".format(encut*ev2ry))
    cmd_list.append("--eps {}".format(epsilon))
    cmd_list.append("--charge {}".format(-charge))
    cmd_list.append("--center {}".format(",".join(defect_pos)))
    cmd_list.append("--relative")
    cmd_list.append("--vdef{}".format(pot_def))
    cmd_list.append("--vref{}".format(pot_ref))
    cmd_list.append("--average {}".format(averaging))
    cmd_list.append("--vasp")
    cmd_list.append("--format=matrix")

    if verbose:
        print("Running:")
        print(" ".join(cmd_list))
    result = subprocess.run(cmd_list, stdout=subprocess.PIPE, text=True)
    result_str = result.stdout

    if verbose:
        print(result_str)
    os.rename("vline-eV-a0.dat", "potential_fit_0.dat")

    # process initial results for shift
    data = np.loadtxt("potential_fit_0.dat")
    x = data[:, 0]
    model = data[:, 1]*np.sign(-charge)
    difference = data[:, 3]

    model_deriv = (model[1:] - model[:-1])/(x[1:] - x[:-1])
    extrema = []
    for i in range(len(model_deriv)-1):
        if np.sign(model_deriv[i+1]) != np.sign(model_deriv[i]):
            extrema.append(i)
    extrema += [0, len(x)-1]

    extreme_vals = [model[i] for i in extrema]
    if verbose:
        print([x[m] for m in extrema])
        print(extreme_vals)
    min_ex_ind = np.argmin(extreme_vals)
    min_ind = extrema[min_ex_ind]
    test_range = 1
    limits = [max(0, np.argmin(np.abs(x-(x[min_ind]-test_range)))), min(np.argmin(np.abs(x - (x[min_ind] + test_range))), len(x)-1)]

    if verbose:
        print("Plateu identified at {}".format(x[min_ind]))

    av_shift = np.mean(difference[limits[0]:limits[1]])
    std_shift = np.std(difference[limits[0]:limits[1]])

    if np.abs(std_shift/av_shift) > 0.5:
        plt.plot(x, model)
        plt.plot(x, difference)
        print("Warning! Standard deviation of shift region high! Make sure potential files are correct.")
        plt.show()

    cmd_list.append("-C {}".format(av_shift))
    result = subprocess.run(cmd_list, stdout=subprocess.PIPE, text=True)
    result_str += "\n" + result.stdout
    if verbose:
        with open("charge_correction.txt", "w") as writer:
            writer.write(result_str)
        print(result_str)

    # extract result
    correction = float(re.findall(r"Defect correction \(eV\): ([\-\d.]+)", result_str)[-1])
    os.chdir(cwd)

    return result_str.strip(), correction


if __name__ == "__main__":
    parser = Parser("Calculates the charge correction of the given VASP system")
    sub_parsers = parser.add_subparsers(help="Choice of correction schemes", dest="scheme")

    LZ_parser = sub_parsers.add_parser("LZ", help="Use the LZ correction")
    LZ_parser.add_argument("outcar", metavar="OUTCAR", help="OUTCAR file of the run")
    LZ_parser.add_argument("eps", metavar="epsilon", type=float, help="Dielectric constant (full)")
    LZ_parser.add_argument("madelung", metavar="madelung", type=float, help="Madelung constant of the lattice (LZ method)")

    FNV_parser = sub_parsers.add_parser("FNV", help="Use the FNV correction")
    FNV_parser.add_argument("outcar", metavar="OUTCAR", help="OUTCAR file of the run")
    FNV_parser.add_argument("eps", metavar="epsilon", type=float, help="Dielectric constant (full)")
    FNV_parser.add_argument("ref_pot", metavar="Reference LOCPOT", help="LOCPOT file of reference system")
    FNV_parser.add_argument("def_pot", metavar="Defect LOCPOT", help="LOCPOT file of defect system")
    FNV_parser.add_argument("-p", dest="def_pos", nargs=3, default=None, help="Relative position of defect in supercell")

    args = parser.parse_args()

    if args.scheme.upper() == "LZ":
        correction = get_LZ_correction(args.outcar, args.eps, args.madelung)
        print("LZ correction: {} eV".format(correction))
    elif args.scheme.upper() == "FNV":
        if args.def_pos is None:
            struct = ase.io.read(args.def_pot, format="vasp")
            ref_struct = ase.io.read(args.ref_pot, format="vasp")
            def_pos = helper.find_defect_pos(struct, ref_struct)
        else:
            def_pos = args.def_pos
        correction_info, _ = get_FNV_correction(args.outcar, args.eps, args.def_pot, args.ref_pot, def_pos)
        print(correction_info)
    else:
        raise Exception("Unknown correction scheme provided")
