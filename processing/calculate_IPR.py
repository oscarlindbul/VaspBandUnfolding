#!/usr/bin/env python

import numpy
import argparse
import glob

from VaspBandUnfolding.vaspwfc import vaspwfc

parser = argparse.ArgumentParser(description="Calculate the IPR of given VASP output")
parser.add_argument("-IBZ", dest="ibzkpt", default="", type=str)
parser.add_argument("-OUT", dest="outcar", default="", type=str)
parser.add_argument("-occup", dest="occupation", default="occupation.ground", type=str)
parser.add_argument("-spin1", dest="spin1", default="gamma.1.ground", type=str)
parser.add_argument("-spin2", dest="spin2", default="gamma.2.ground", type=str)
parser.add_argument("-WAV", dest="wav", default="", type=str)
band_group = parser.add_mutually_exclusive_group()
parser.add_argument("-range", dest="ipr_range", default=15, type=int)
parser.add_argument("--gamma", dest="gamma", default=False, action="store_true")
band_group.add_argument("-center", dest="spin_center", default=None, type=int)
band_group.add_argument("-bands", dest="bands", default=None, type=int, nargs="+", help="Manual specification of bands to calculate IPR")
parser.add_argument("-q", dest="quiet", default=False, action="store_true")

input = parser.parse_args()

#if input.ibzkpt == "":
#    input.ibzkpt = glob.glob("*IBZKPT*")[0]
#if input.outcar == "":
#    input.outcar = glob.glob("*OUTCAR*")[0]
if input.wav == "":
    input.wav = glob.glob("*WAVECAR*")[0]

# read files
wav = vaspwfc(input.wav, lgamma=input.gamma)

if input.spin_center is None and input.bands is None:
	occupation=open(input.occupation,'r')
	spin_1=open(input.spin1,'r')
	spin_2=open(input.spin2,'r')
	# find number of bands
	occupation_lines = occupation.readlines()
	spin_1_lines = spin_1.readlines()
	spin_2_lines = spin_2.readlines()

	bands = int(spin_2_lines[-1].split()[0])

	# find highest occupied band (HOB)
	HOB_1 = 0
	HOB_2 = 0

	for line in spin_1_lines:
		if line.split()[2]=='0.00000':
			HOB_1=int(line.split()[0])-1
			break

	for line in spin_2_lines:
		if line.split()[2]=='0.00000':
			HOB_2=int(line.split()[0])-1
			break

	# close files
	occupation.close()
	spin_1.close()
	spin_2.close()
else:
	if input.spin_center is not None:
		HOB_1 = input.spin_center
		HOB_2 = input.spin_center
		bands = [min(HOB_1,HOB_2) - b for b in range(-input.ipr_range,input.ipr_range+1)]
	elif input.bands is not None:
		bands = input.bands


# find and save IPR
ipr = wav.inverse_participation_ratio(bands=bands, quiet=input.quiet)
numpy.save("ipr.npy", ipr)
