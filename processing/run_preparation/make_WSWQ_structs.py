import os, re
from pathlib import Path
from shutil import copyfile
import argparse, glob

parser = argparse.ArgumentParser("Prepare runs for nonradiative capture calculations")
parser.add_argument("ground_state", metavar="ground path", help="Path to ground state")
parser.add_argument("ground_runs", metavar="ground runs", help="Path to ground state run files")
parser.add_argument("out_path", metavar="output path", help="Path to directory to contain output files")
args = parser.parse_args()


# equilibrium structures from your first-principles calculation
real_ground_dir = Path(args.ground_state)
ground_dirs = glob.glob(args.ground_runs + "/*")
ground_dirs.sort(key=lambda x: int(x.split("/")[-1]))

#run_files = ['KPOINTS', 'POTCAR', 'INCAR', 'POSCAR', 'WAVECAR', 'CHGCAR']
cc_dir = Path(args.out_path)

run_files = ["INCAR", "KPOINTS", "POTCAR", "POSCAR", "WAVECAR"]

#real_ground_dir = Path(args.ground_state) / str(len(ground_dirs)//2)
for i, d in enumerate(ground_dirs):
	print("Copying files {}".format(i+1))
	working_dir = cc_dir / str(i)
	if not os.path.exists(working_dir):
		os.mkdir(str(working_dir))

	# write structure and copy necessary input files
	for f in run_files:
		f_out = f
		copyfile(str(real_ground_dir / f), str(working_dir / f_out))
	
	copyfile(str(Path(d) / "WAVECAR"), str(working_dir / "WAVECAR.qqq"))
	with open(str(working_dir / 'INCAR'), "r") as reader:
		lines = reader.readlines()
	for i in range(len(lines)):
		if "ALGO" in lines[i]:
			lines[i] = "ALGO = None\n"
	lines.append("LWSWQ = True")
	with open(str(working_dir / 'INCAR'), "w") as writer:
		writer.writelines(lines)
