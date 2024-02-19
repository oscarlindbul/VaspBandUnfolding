import os, sys, shutil
import numpy as np
from scipy import constants
import subprocess

import argparse as arg
import traceback

import glob
import bz2
import re

from helperfuncs import iterate_ground_states, convert_charge, open_file_with_check, search_for

parser = arg.ArgumentParser(description="Calculates the FNV correction of ground state formation energies")

#parser.add_argument("defect_name", metavar="name", type=str, nargs="+", help="name of defect to search for")
parser.add_argument("-c", dest="charge", type=int, required=True, help="Option to specify charge")
#parser.add_argument("-s", dest="spin", nargs="+", type=int, default=None, help="Option to specify spin runs to search for (default all)")
parser.add_argument("-eps", dest="eps", type=float, required=True, help="The dielectric constant of bulk material")
#parser.add_argument("-R", dest="run_path", type=str, default="../Runs", help="Option to specify data folder")
parser.add_argument("-bulk", dest="bulk_path", type=str, required=True, help="The search path for bulk simulation data")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-posfile", dest="pos_file", type=str, default=None, help="Option to specify files containing defect positions in cell")
group.add_argument("-pos", dest="pos", type=float,nargs=3, default=None, help="Option to specify defect position manually")
parser.add_argument("-Cpath", dest="C_path", type=str, default=None, help="Give correction data through file with row format (state C averaging)")
parser.add_argument("-C", dest="C", type=float, default=None, help="Option to specify defect align shift")
parser.add_argument("-average", dest="average", type=float, default=None, help="Option to average potentials")
parser.add_argument("--no-log", dest="log", default=True, action="store_false", help="Option to not save FNV data (in case of testing etc.)")
parser.add_argument("-beta", dest="beta", default=None, type=float, help="Beta parameter to use in FNV model (determines localization of defect)")

input = parser.parse_args(sys.argv[1:])

write_log = input.log and (input.C_path is not None or input.C is not None)
log_file_written = False


##### Load bulk LOCPOT data
bulk_path = input.bulk_path
pot_file = search_for("{}/LOCPOT*".format(bulk_path))
state_zipped_files = set()
bulk_zipped_files = set()
print("Decompressing bulk potential data...")
if ".bz2" in pot_file:
    print("... {}".format(pot_file))
    zip_process = subprocess.Popen(["bunzip2", pot_file])
    zip_process.wait()
    bulk_zipped_files.add(pot_file.replace(".bz2", ""))
print("Done")
try:
    bulk_pot_file = search_for("{}/LOCPOT*".format(bulk_path))

    ################ Load defect data

    #### Load defect position in cell
    if input.pos is not None:
        defect_pos = tuple(input.pos)
    if input.pos_file is not None:
        print("Reading cell position data...")
        with open(input.pos_file, "r") as reader:
            lines = reader.readlines()
            for line in lines:
                defect_pos = [float(p) for p in line.split()[1:]]
        print("Done")
            
    print("Reading cutoff energy data from OUTCAR...")
    outcar_file = search_for("./OUTCAR*")
    with open_file_with_check(outcar_file, "r") as reader:
        for line in reader.readlines():
            search = re.search("energy-cutoff\s+:\s+([.\d]+)", line)
            if search is not None:
                E_cutoff = float(search.group(1))
                break
        
    n_elec = -input.charge

    #### Find cutoff E
#               incar_file = glob.glob("{}/INCAR.ground*".format(ground_sub_folder))
    #if len(incar_file) == 0: # if incar not in data folder, check run folder
    #    incar_file = glob.glob("{}/INCAR.ground*".format(ground_folder))
    #incar_file = incar_file[0]
    #incar_reader = open_file_with_check(incar_file)
    #cutoff_pattern = re.compile("^ENCUT\s*=\s*([-\d.]+)")
    #for line in incar_reader.readlines():
    #    match = cutoff_pattern.match(line)
    #    if match is not None:
    #        E_cutoff = float(match.group(1))
    #incar_reader.close()
    print("Done")

    #### Find real spin value for later categorizing
    #analyze_file = glob.glob("{}/analyse.txt*".format(ground_sub_folder))
    #if len(analyze_file) == 0:
    #    raise Exception("No analyze file found for groundstate ({},{}) in {}".format(charge, spin, ground_sub_folder))
    #analyze_file = analyze_file[0]
    #with open_file_with_check(analyze_file) as analyze_reader:
    #    HOB_pattern = re.compile("spin (\d) HOB\s*:\s*(\d+)")
    #    HOB = [0,0]
    #    for line in analyze_reader.readlines():
    #        HOB_match = HOB_pattern.match(line)
    #        if HOB_match is not None:
    #            spin_ind = int(HOB_match.group(1))-1
    #            HOB_val = int(HOB_match.group(2))
    #            HOB[spin_ind] = HOB_val
    #real_spin = (HOB[0] - HOB[1])/2.0

    # Conversion to Ry
    Ry_E = constants.Rydberg*constants.h*constants.c
    e = constants.e
    E_cutoff *= (e/Ry_E) # energy in Ry

    #### Load eps
    eps = input.eps

    #### Load defect state LOCPOT data
    pot_files = search_for("./LOCPOT*")
    zipped_files = []
    print("Decompressing ground state potential data... ")
    for pot_file in pot_files:
        if ".bz2" in pot_file:
            print("... {}".format(pot_file))
            zip_process = subprocess.Popen(["bunzip2", pot_file])
            zip_process.wait()
            state_zipped_files.add(pot_file.replace(".bz2", ""))
    print("Done")
    def_pot_file = search_for("./LOCPOT*")

    print("Calculating FNV corrections and plots...")
    script_path = os.path.realpath(os.path.dirname(__file__))
    defect_align_command = [script_path + "/sxdefectalign"]
    defect_align_command.extend(["--ecut", str(E_cutoff)])
    defect_align_command.extend(["--charge", str(n_elec)])
    defect_align_command.extend(["--eps", str(eps)])
    if input.beta is not None:
        defect_align_command.extend(["--beta", str(input.beta)])
    defect_align_command.extend(["--relative"])
    defect_align_command.extend(["--center", str(defect_pos).strip("() ").replace(" ", "")])
    defect_align_command.extend(["--vdef", def_pot_file])
    defect_align_command.extend(["--vref", bulk_pot_file])
    
    if input.average is not None:
        averaging = input.average
        defect_align_command.extend(["--average", str(averaging)])
    if input.C is not None:
        correction = input.C
        defect_align_command.extend(["-C", str(correction)])
    defect_align_command.append("--vasp")

    ex_string = " ".join(defect_align_command)
    print("Executing:\n" + ex_string)
    FNV_process = subprocess.Popen(ex_string, stdout=subprocess.PIPE, shell=True)
    output, stderror = FNV_process.communicate()
    output = str(output)
    print(output)
    print(stderror)
    if stderror is not None:
        raise Exception(stderror)
    print("Done")

    print("Saving FNV data and plots")

    ##### Parse output data
    E_iso_pattern = re.compile("Isolated energy\s*:\s*([\-.\d]+)")
    E_period_pattern = re.compile("Periodic energy\s*:\s*([\-.\d]+)")
    E_diff_pattern = re.compile("Difference \(eV\)\s*:\s*([\-.\d]+)")
    used_eps_pattern = re.compile("Calculation performed with epsilon = ([.\d]+)")
    align_correction_pattern = re.compile("Defect correction \(eV\): ([\-.\d]+)\s*\(incl\. screening & alignment\)")

    E_iso_match = E_iso_pattern.search(output)
    E_period_match = E_period_pattern.search(output)
    E_diff_match = E_diff_pattern.search(output)
    used_eps_match = used_eps_pattern.search(output)
    align_correction_match = align_correction_pattern.search(output)

    if E_iso_match is None or E_period_match is None \
        or E_diff_match is None or used_eps_match is None \
        or align_correction_match is None:
        raise Exception("Output format of sxdefectalign script not according to standard. Is version not 2.2?")
    else:
        E_iso = float(E_iso_match.group(1))
        E_period = float(E_period_match.group(1))
        E_diff_match = float(E_diff_match.group(1))
        used_eps = float(used_eps_match.group(1))
        align_correction = float(align_correction_match.group(1))

    ###### Save FNV data

    save_folder = "."
    if not os.path.isdir(save_folder):
        os.mkdir(save_folder)
    created_plot_name = "vline-eV-a0.dat"
    if input.average is None:
        if input.C is None:
            category_name = "FNVplot_dat"
        else:
            category_name = "corrected.FNVplot_dat"
    else:
        if input.C is None:
            category_name = "avg_{}.FNVplot_dat".format(str(input.average))
        else:
            category_name = "avg_{}.corrected.FNVplot_dat".format(str(input.average))
    os.rename(created_plot_name, category_name)
    shutil.move(category_name, os.path.join(save_folder, category_name))
    unnecessary_files = ["vline-eV-a1.dat", "vline-eV-a2.dat", "vAtoms.dat", "fftwisdom.dat"]
    for file_path in unnecessary_files:
        if os.path.isfile(file_path):
            os.remove(file_path)

    if write_log:
        log_file = "{}/FNV_dat".format(save_folder)
        log_modifier = "w" if not log_file_written else "a"
        with open(log_file, log_modifier) as log_writer:
            log_writer.write("Potential Correction = {}\n".format(correction))
            log_writer.write("Potential Averaging = {} (Bohr)\n".format(averaging))
            log_writer.write("Isolation Energy = {} eV\n".format(E_iso))
            log_writer.write("Periodic Energy = {} eV\n".format(E_period))
            log_writer.write("Energy Difference = {} eV\n".format(E_diff_match))
            log_writer.write("Used EPS = {}\n".format(used_eps))
            log_writer.write("Lattice Correction = {} eV\n".format(E_iso - E_period))
            log_writer.write("Total Correction = {} eV\n".format(align_correction))
            log_writer.write("Proposed Total Correction = {} eV\n\n".format(E_iso - E_period + correction))
            log_file_written = True
    print("Done")

    print("Recompressing ground state potential data...")
    for file in state_zipped_files:
        print("... {}".format(file))
        zip_process = subprocess.Popen(["bzip2", file])
        zip_process.wait()
    state_zipped_files.clear()
    print("Done")
except Exception as e:
    traceback.print_exc()
finally:
    print("Recompressing bulk potential data...")
    for file in bulk_zipped_files.union(state_zipped_files):
        print("... {}".format(file))
        zip_process = subprocess.Popen(["bzip2", file])
        zip_process.wait()
    print("Done")
