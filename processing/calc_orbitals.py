from VaspBandUnfolding.vaspwfc import vaspwfc

import argparse
import glob
from helperfuncs import search_for

parser = argparse.ArgumentParser("Calculates the pseudo-wavefunctions of given bands and saves them in vesta format")
parser.add_argument("bands", metavar="bands", nargs="+", type=int, help="state bands to extract")
parser.add_argument("spin", metavar="spin", type=int, help="Spin channel to use, 1(up) or 2(down)")
parser.add_argument("-cell", dest="cell", default=None, type=str, help="Supercell to use as reference")
parser.add_argument("-wavecar", dest="wavecar", default=None, type=str, help="Waecar to take wavefunctions from")
parser.add_argument("--gamma", dest="gamma", default=False, action="store_true", help="Specify if wavefunction is gamma")
parser.add_argument("-gam_half", dest="gamma_half", default="x", type=str, help="Type of gamma cutoff, x or z")
parser.add_argument("-kp", dest="k_point", default=1, type=int, help="Which k-point to use")
#parser.add_argument("--density", dest="density", default=False, action="store_true", help="Calculate the charge density instead of wavefunction")
parser.add_argument("--imag", dest="imaginary", default=False, action="store_true", help="Write out imaginary part of wavefunction (default is to write real part only)")
parser.add_argument("--density", dest="calc_density", default=False, action="store_true", help="Write out the state density instead of wavefunction")

args = parser.parse_args()

if args.cell is None:
	cells = search_for("*CONTCAR*")
	if cells is None:
		cells = search_for("*POSCAR*")
	if cells is None:
		raise ValueError("No valid CONTCAR or POSCAR found! Specify with -cell")
	args.cell = cells

if args.wavecar is None:
	wavs = search_for("*WAVECAR*")
	if wavs is None:
		raise Exception("No valid WAVECAR found! Specify with -wavecar")
	args.wavecar = wavs

print("Calculating orbitals:")
print("Structure file: {}".format(args.cell))
print("Wavefunction file: {}".format(args.wavecar))

wav = vaspwfc(args.wavecar, lgamma=args.gamma, gamma_half=args.gamma_half)
for band in args.bands:
	wavefunc = wav.wfc_r(args.spin, args.k_point, band)
	argument = wavefunc
	if args.calc_density:
		argument = wavefunc.real**2 + wavefunc.imag**2
	wav.save2vesta(argument, lreal=not args.imaginary, poscar=args.cell, prefix="".join(["wav_s"+str(args.spin)+"_", str(band)]))
