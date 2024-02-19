import glob
import re, os
import bz2, gzip
import numpy as np

#from ase.build import minimize_rotation_and_translation as optimize_pos
from ase import geometry

from collections import defaultdict

def convert_charge(charge, convention="file"):
    """
    Converts the charge convetion "+", "0", "-" to corresponding integers
    """
    string_convention1 = re.compile("[0+\-]+")
    string_convention2 = re.compile("c[np]?\d")
    if isinstance(charge, str) or isinstance(charge, unicode):
        if string_convention1.match(charge) is not None:
            c = 0
            for i in range(len(charge)):
                if charge[i] == '+':
                    c += 1
                elif charge[i] == '-':
                    c -= 1
                elif charge[i] == '0':
                    c = 0
                    break
            return c
        elif string_convention2.match(charge) is not None:
            if "cp" in charge:
                return int(charge[2])
            elif "cn" in charge:
                return -int(charge[2])
            else:
                return 0
        else:
            Exception("Unknown string charge convention")
    elif isinstance(charge, int):
        if convention == "file":
            if charge < 0:
                return "cn" + str(abs(charge))
            elif charge > 0:
                return "cp" + str(charge)
            else:
                return "c0"
        elif convention == "compact":
            if charge == 0:
                return "0"
            elif charge > 0:
                return "+"*charge
            else:
                return "-"*(-charge)
        else:
            Exception("Unknown integer charge convention")


def iterate_ground_states(run_folder, charges=None, spins=None):
    run_sub_folder = glob.glob("{}/ht.run.*".format(run_folder))[0]

    if charges is None or spins is None:
        # find all present states
        charge_set = set()
        spin_set = set()
        charge_set.add(0)   # unpolarized ground state always present
        spin_set.add(0)

        # match with format cp1-s2 for example and extract 1 and 2
        charged_pattern = re.compile(".*\\.(c(?:n|p)?\d)\\..*")
        polarized_pattern = re.compile(".*\\.(c(?:n|p)?\d)-s(\\d)\\..*")
        for root, state_names, files in os.walk(run_sub_folder):
            for state_name in state_names:
                m = charged_pattern.match(state_name)
                if m is not None:
                    charge_set.add(convert_charge(m.group(1)))
                m = polarized_pattern.match(state_name)
                if m is not None:
                    charge_set.add(convert_charge(m.group(1)))
                    spin_set.add(int(m.group(2)))
        if charges is None:
            charges = sorted(list(charge_set))
        if spins is None:
            spins = sorted(list(spin_set))

    for charge in charges:
        charge_id = convert_charge(charge)
        for spin in spins:
            if charge == 0 and spin == 0:
                state_folder = run_folder
                data_folder = run_sub_folder
            else:
                ground_id = charge_id
                if not spin == 0:
                    ground_id += "-s" + str(spin)
                state_folders = glob.glob("{}/*.{}.*".format(run_sub_folder, ground_id))
                if len(state_folders) > 1:
                    raise Exception("duplicate ground states found in run folder")
                elif len(state_folders) == 0:
                    continue
                state_folder = state_folders[0]
                data_folder = glob.glob("{}/ht.run.*".format(state_folder))[0]

            yield data_folder, state_folder, charge, spin


def iterate_excited_states(run_folder, charge, spin):
    run_sub_folder = glob.glob("{}/ht.run.*".format(run_folder))[0]
    if spin == 0:
        ground_ID = convert_charge(charge, "file")
    else:
        ground_ID = convert_charge(charge, "file") + "-s" + str(spin)

    # find all present states
    ex_type_set = set()

    # match with format cp1-s2-k1 for example and extract k and 1
    pattern = re.compile(".*\\.{}-([klv])(\\d)\\..*".format(ground_ID))
    for root, state_names, files in os.walk(run_sub_folder):
        for state_name in state_names:
            m = pattern.match(state_name)
            if m is not None:
                ex_type_set.add((m.group(1), int(float(m.group(2)))))

    ex_types = sorted(list(ex_type_set))

    for ex_type in ex_types:
        ex_ID = ground_ID + "-{}{}".format(ex_type[0], ex_type[1])
        state_folder = glob.glob("{}/*.{}.*".format(run_sub_folder, ex_ID))[0]
        data_folder = glob.glob("{}/ht.run*".format(state_folder))[0]

        yield data_folder, state_folder, ex_type[0], ex_type[1]


def open_file_with_check(filename, option="r"):
    if ".bz2" in filename:
        return bz2.BZ2File(filename, option)
    elif ".gz" in filename:
        return gzip.GzipFile(filename, mode=option)
    else:
        return open(filename, option)


def search_for(f_name):
    results = glob.glob(f_name)
    if len(results) > 1:
        print("Found several results:")
        for i in range(len(results)):
            print("{}) {}".format(i+1, results[i]))
        prompt = int(input("Which would be preferred?"))
        return results[prompt-1]
    elif len(results) == 1:
        return results[0]
    else:
        return None


def write_lists_to_csv(name, headers, sep, *lists):
    file = open(name, "w")
    file.write(sep.join(headers) + sep + "\n")
    i = 0
    keep_writing = True
    while keep_writing:
        keep_writing = False
        line = ""
        for values in lists:
            if i < len(values):
                keep_writing = True
                line += str(values[i])
            line += sep
        line += "\n"
        i += 1
        if keep_writing:
            file.write(line)
    file.close()


def get_charge_OUTCAR(outcar):
    zval_pattern = re.compile(r"Ionic Valenz\s*\n\s*ZVAL\s*\=\s*(.*)")
    natom_pattern = re.compile(r"ions per type\s*\=\s*([\d\s]+)")
    nelect_pattern = re.compile(r"NELECT\s*=\s*([\d.]+)")
    with open_file_with_check(outcar) as outcar_reader:
        file_data = outcar_reader.read()
    zval_match = zval_pattern.search(file_data)
    natom_match = natom_pattern.search(file_data)
    nelect_match = nelect_pattern.search(file_data)
    if zval_match is None or natom_match is None:
        raise Exception("Valence or atom data not found in OUTCAR")

    zvals = list(map(int,map(float, zval_match.group(1).split())))
    natoms = list(map(int, natom_match.group(1).split()))
    nelect = float(nelect_match.group(1))
    electrons = sum([i*j for i, j in zip(zvals, natoms)])
    charge = electrons - nelect
    return charge


def get_final_energy(outcar):
    energy_pattern = re.compile(r"energy without entropy\s*\=\s*([-e.\d]+).*")
    with open_file_with_check(outcar) as reader:
        matches = energy_pattern.findall(reader.read())
    if len(matches) <= 0:
        raise Exception("Format of OUTCAR ({}) incorrect! No energy found".format(outcar))
    energy = float(matches[-1])
    return energy


def get_ion_content(outcar):
    ion_type_n_pattern = re.compile(r"ions per type \=\s+(.*)")
    ion_type_pattern = re.compile(r"VRHFIN\s*=\s*([A-Za-z]+)")

    ion_n = defaultdict(lambda: 0)
    with open_file_with_check(outcar) as reader:
        outcar_data = reader.read()
    type_matches = ion_type_pattern.findall(outcar_data)
    n_match = ion_type_n_pattern.search(outcar_data)

    if len(type_matches) <= 0 or n_match is None:
        raise Exception("Format of OUTCAR incorrect! Ion content not found")

    ion_amounts = list(map(int, n_match.group(1).split()))
    for i, ion_match in enumerate(type_matches):
        ion_n[ion_match] += ion_amounts[i]

    return ion_n


def get_volume(outcar):
    volume_pattern = re.compile(r"volume of cell\s+:\s+([\d.]+)")

    with open_file_with_check(outcar) as reader:
        matches = volume_pattern.findall(reader.read())
    if len(matches) <= 0:
        raise Exception("Format of OUTCAR incorrect! No volume data found")

    volume = float(matches[-1])
    return volume


def get_encut(outcar):
    encut_pattern = re.compile(r"ENCUT\s*\=\s*([\d.]+)")

    with open_file_with_check(outcar) as reader:
        matches = encut_pattern.findall(reader.read())
    if len(matches) <= 0:
        raise Exception("Format of OUTCAR incorrect! No volume data found")

    encut = float(matches[-1])
    return encut


def get_bands(outcar):
    eigval_pattern = re.compile(r"(spin component \d\n[a-zA-Z\d.\-\s\n\:]*)-------")
    with open_file_with_check(outcar) as outcar_file:
        data = outcar_file.read()
    nband_match = re.search(r"NBANDS\s*\=\s*(\d+)", data)
    kpts_match = re.search(r"NKPTS\s*\=\s*(\d+)", data)
    ispin_match = re.search(r"ISPIN\s*\=\s*(\d)", data)
    if nband_match is None or kpts_match is None or ispin_match is None:
        raise Exception("Failed to read nbands, number of kpoints or spin settings in OUTCAR")
    nbands = int(nband_match.group(1))
    nkpts = int(kpts_match.group(1))
    ispin = int(ispin_match.group(1))
    eig_reading = eigval_pattern.findall(data)
    if len(eig_reading) < 1:
        raise Exception("Format of OUTCAR incorrect! Error finding eigenvalues. Do you print out eigenvalues in OUTCAR?")

    spin_pattern = re.compile(r"spin component (\d)([a-zA-Z\d\-.\s\:\n]*?)(?=spin component|------)")
    k_pattern = re.compile(r"k-point\s*(\d+)\s*:.*\n.*\n([a-zA-Z\d\-.\s\n\:]*)")
    spin_readings = spin_pattern.findall(eig_reading[-1])

    bands = np.zeros((ispin, nkpts, nbands))
    occups = np.zeros((ispin, nkpts, nbands))
    if len(spin_readings) < 1:
        raise Exception("Format of OUTCAR incorrect! Error finding index of spin eigenvalues")

    for s_reading in spin_readings:
        s_ind = int(s_reading[0]) - 1
        k_readings = k_pattern.findall(s_reading[1])
        if len(k_readings) < 1:
            raise Exception("Format of OUTCAR incorrect! Error reading k-point eigenvalues")

        for k_reading in k_readings:
            k_ind = int(k_reading[0]) - 1
            eigenvals = k_reading[1].split("\n")
            for band_ind in range(nbands):
                eigs = eigenvals[band_ind]
                vals = list(map(float, eigs.split()))
                bands[s_ind, k_ind, band_ind] = float(vals[1])
                occups[s_ind, k_ind, band_ind] = float(vals[2])

    return bands, occups


def find_defect_pos(structure, ref_structure, threshold = 0.7):
    # TODO: Structure-based approach did not work...
    vacancy_weight = np.cbrt(structure.get_volume())
    substitution_weight = vacancy_weight
    if structure.get_volume() - ref_structure.get_volume() > 5:
        raise Exception("Volumes of input structures differ too much. Structures likely incompatible")
    # optimize_pos(structure, ref_structure)

    matched_ref = [False]*len(ref_structure)
    matched_source = [False]*len(structure)
    deviation_sum = 0
    deviations = np.zeros((len(structure), 2, 3))
    center_of_deviation = np.zeros(3)
    rms = 0
    for j, search_atom in enumerate(structure):
        min_dev_pos = None
        min_dev_norm = np.inf
        min_dev = None
        min_ind = 0
        for i, ref_atom in enumerate(ref_structure):
            if search_atom.symbol != ref_atom.symbol or matched_ref[i]:
                continue
            test_deviation = geometry.get_distances(search_atom.position, ref_atom.position, structure.cell, pbc=True)[0][0][0]
            if np.linalg.norm(test_deviation) < min_dev_norm:
                min_dev_norm = np.linalg.norm(test_deviation)
                min_dev_pos = ref_atom.position
                min_dev = test_deviation
                min_ind = i

        if min_dev_norm < threshold:
            deviation_sum += min_dev
            center_of_deviation += min_dev*min_dev_pos
            deviations[j, 0, :] = min_dev
            deviations[j, 1, :] = min_dev_pos
            matched_ref[min_ind] = True
            matched_source[j] = True

    for i in range(len(matched_ref)):
        if not matched_ref[i]:
            deviations[i, 0, :] = vacancy_weight
            deviations[i, 1, :] = vacancy_weight*np.array([1, 1, 1])
            deviation_sum += vacancy_weight
            center_of_deviation = vacancy_weight*structure.get_positions()[i, :]
    for i in range(len(matched_source)):
        if not matched_source[i]:
            deviations[i, 0, :] = substitution_weight
            deviations[i, 1, :] = substitution_weight*np.array([1, 1, 1])
            deviation_sum += substitution_weight
            center_of_deviation = substitution_weight*structure.get_positions()[i, :]

    center_of_deviation /= deviation_sum
    return structure.cell.scaled_positions(center_of_deviation)
