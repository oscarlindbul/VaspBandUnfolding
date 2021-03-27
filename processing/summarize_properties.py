import numpy as np
import argparse
import subprocess, glob, os, re
from helperfuncs import open_file_with_check as open_file

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

parser = argparse.ArgumentParser("")
parser.add_argument("paths", metavar="paths", nargs="+", help="Paths to directories in order ground, (relax), final")
parser.add_argument("--gamma", dest="gamma", default=False, action="store_true")
parser.add_argument("-bands", dest="bands", default=None, nargs=4)
parser.add_argument("--force-recalc", dest="recalc", default=False, action="store_true", help="Forces recalculation of input data such as TDM")
parser.add_argument("-wav-prefix", dest="wav_prefix", default="WAVECAR", help="Prefix to use for searching wavecars")
args = parser.parse_args()
ground_path = args.paths[0]
if len(args.paths) < 3:
	relax_path = None
	final_path = args.paths[1]
else:
	relax_path, final_path = args.paths[1:]
gamma_option = ""
if args.gamma:
	gamma_option = "--gamma"

script_path="/afs/pdc.kth.se/home/o/oscarlb/Public/py_tools/VaspBandUnfolding/processing/"

## Find structural differences
ground_cont = search_for(ground_path + "/" + "*CONTCAR*")
ex_cont = search_for(final_path + "/" + "*CONTCAR*")
output = os.popen('python {}/structure_compare.py {} {}'.format(script_path, ground_cont, ex_cont)).read()

pr_search = re.search("Participation: ([\d.]+)", output)
structdiff_mean = float(pr_search.group(1))
max_search = re.search("Position difference \(abs max\): ([\d.]+) A", output)
structdiff_max = float(max_search.group(1))
sum_search = re.search("Position difference \(RMS\): ([\d.]+) A", output)
structdiff_sum = float(sum_search.group(1))

## Find ZPL
ground_osz = search_for(ground_path + "/" + "*OSZICAR*")
ex_osz = search_for(final_path + "/" + "*OSZICAR*")
with open_file(ground_osz, "r") as reader:
	tot_e_match = re.search("E0= ([\-\d.E+]+)", reader.readlines()[-1])
	ground_E = float(tot_e_match.group(1))
with open_file(ex_osz, "r") as reader:
	tot_e_match = re.search("E0= ([\-\d.E+]+)", reader.readlines()[-1])
	ex_E = float(tot_e_match.group(1))
ZPL = ex_E - ground_E

## Find polarization & lifetime
ground_wav = search_for(ground_path + "/" + args.wav_prefix + "*")
ground_eig = search_for(ground_path + "/" + "*EIGENVAL*")
ex_wav = search_for(final_path + "/*" + args.wav_prefix + "*")
ex_eig = search_for(final_path + "/" + "*EIGENVAL*")
tdm_file = glob.glob(final_path + "/" + "*tdm*")

if len(tdm_file) < 1 or args.recalc:
	recomp_ground = ".bz2" in ground_wav
	recomp_ex = ".bz2" in ex_wav
	### Calculate tdm
	if recomp_ground:
		print("decompressing ground wavecar...")
		os.system("bunzip2 {}".format(ground_wav))
		ground_wav = ground_wav.replace(".bz2", "") 
	if recomp_ex:
		print("decompressing ex wavecar...")
		os.system("bunzip2 {}".format(ex_wav))
		ex_wav = ex_wav.replace(".bz2", "")
	output = final_path + "/tdm.npz"
	os.system("python {}/calculate_tdm.py {} {} {} -out {}".format(script_path, ground_wav, ex_wav, gamma_option, output))
	tdm_file = glob.glob(final_path + "/" + "*tdm*")
	if recomp_ground:
		print("recompressing ground wavecar...")
		os.system("bzip2 {}".format(ground_wav))
	if recomp_ex:
		print("recompressing ex wavecar...")
		os.system("bzip2 {}".format(ex_wav))
else:
	tdm_file = [search_for(final_path + "/" + "*tdm*")]
tdm_file = tdm_file[0]

# find transition
if args.bands is None:
		ground_occup = []
		ex_occup = []
		for f, occup in zip((ground_eig, ex_eig), (ground_occup, ex_occup)):
			with open_file(f, "r") as reader:
				lines = reader.readlines()[8:]
				for l in lines:
					if len(l.split()) < 5:
						break
					occ = list(map(float,l.split()[3:5]))
					occup.append(occ)
		ground_occup = np.array(ground_occup)
		ex_occup = np.array(ex_occup)
		print(len(ground_occup), len(ex_occup))
		ids = np.where(ground_occup != ex_occup)
		from_band = ids[0][0]
		to_band = ids[0][1]
		from_spin = ids[1][0]
		to_spin = ids[1][1]
else:
	from_band = int(args.bands[0])
	from_spin = int(args.bands[1])
	to_band = int(args.bands[2])
	to_spin = int(args.bands[3])
# parse data
tdm_data = os.popen("python {}/calc_lifetime.py {} {} {} {} {} {} {} {}".format(script_path, tdm_file, from_spin, 0, from_band+1, to_spin, 0, to_band+1, ZPL)).read()
pol_search = re.search("Dipole Moment: \(([-+e\d.]+),([-+e\d.]+),([-+e\d.]+)\)", tdm_data)
pol = list(map(float, pol_search.group(1,2,3)))
lifetime_search = re.search("Lifetime: ([-+e\d.]+)", tdm_data)
lifetime = float(lifetime_search.group(1))*1e9 # nanoseconds
osc_search = re.search("Oscillator strength: ([-+e\d.]+)", tdm_data)
oscillator_strength = float(osc_search.group(1))
angle_search = re.search("Pol angle .+: ([-+e\d.]+)", tdm_data)
angle = float(angle_search.group(1))

## Find energy relaxation
relax_approx = None
relax_exact = None
if relax_path is not None:
	relax_osz = search_for(relax_path + "/*OSZICAR*")
	with open_file(relax_osz, "r") as reader:
		relax_search = re.findall("E0= ([\-\d.E+]+)", "\n".join(reader.readlines()))
	first, last = float(relax_search[0]), float(relax_search[-1])
	relax_approx = np.abs(last - first)
	relax_exact = np.abs(ex_E - first)

print("Transition: {} to {} in channel {}".format(from_band+1, to_band+1, from_spin))
print("ZPL: {:.4f}".format(ZPL))
print("Polarization: {:.4E} {:.4E} {:.4E}".format(pol[0], pol[1], pol[2]))
print("Pol angle: {:.2f}".format(angle))
print("Lifetime: {:.3f}".format(lifetime))
print("Oscillator Strength: {:.3e}".format(oscillator_strength))
print("Struct dev (PR): {:.4f}".format(structdiff_mean))
print("Struct dev (max): {:.4f}".format(structdiff_max))
print("Struct dev (RMS): {:.4f}".format(structdiff_sum))
if relax_path is not None:
	print("E relax approx: {:.4f}".format(relax_approx))
	print("E relax exact: {:.4f}".format(relax_exact))

