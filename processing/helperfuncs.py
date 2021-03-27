import glob
import re, os
import bz2, gzip

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
        charge_set.add(0) # unpolarized ground state always present
        spin_set.add(0)

        # match with format cp1-s2 for example and extract 1 and 2
        charged_pattern = re.compile(".*\.(c(?:n|p)?\d)\..*")
        polarized_pattern = re.compile(".*\.(c(?:n|p)?\d)-s(\d)\..*")
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
    pattern = re.compile(".*\.{}-([klv])(\d)\..*".format(ground_ID))
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
	else:
		return results[0]

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
