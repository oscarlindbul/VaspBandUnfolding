#!/usr/bin/env python

from VaspBandUnfolding.vaspwfc import vaspwfc

import glob
import bz2
import argparse

import numpy as np

parser = argparse.ArgumentParser("Calculates the tdm data between given WAVECARS and corresponding input files")
parser.add_argument("from_wav", metavar="from_wav", help="initial state wavecar")
parser.add_argument("to_wav", metavar="to_wav", help="end state wavecar")
parser.add_argument("-bands", dest="bands", type=int, nargs="+", default=[1149, 1150, 1151, 1152, 1153, 1154, 1155], help="Bands to consider for transitions")
#parser.add_argument("-IBZ", dest="IBZKPT", default="IBZKPT", help="K-point structure file")
parser.add_argument("-out", dest="out_name", default="tdm", help="Name of output file")
parser.add_argument("--gamma", dest="gamma", default=False, action="store_true")
parser.add_argument("--same-wav", dest="same", default=False, action="store_true")
parser.add_argument("-states", dest="states", nargs=6, help="States to explicitly calculate tdm for, in the format (spin1, k1, band1, spin2, k2, band2)")
in_args = parser.parse_args()


#ibzkpt_path = in_args.IBZKPT

# find kpoints
#ibzkpt=open(ibzkpt_path,'r')
#ibzkpt.readline()
#kpoints=int(ibzkpt.readline())
#ibzkpt.close()

# get parchial charge bands
parchg = in_args.bands

# find wavecars
from_wav_path = in_args.from_wav
to_wav_path = in_args.to_wav

# from ground to tag step
from_wav = vaspwfc(from_wav_path, lgamma=in_args.gamma)
to_wav = vaspwfc(to_wav_path, lgamma=in_args.gamma)

# make output array to save dE, overlap and tdm
output=np.empty((2,1,len(parchg),len(parchg),5),dtype=np.complex_)

for spin in range(2):
	for k in [0]:
		for i in range(len(parchg)):
			for j in range(len(parchg)):
				tdm = [0, 0, (0, 0, 0)]
				if in_args.same and i != j:
					tdm = from_wav.TransitionDipoleMoment((spin+1, k+1, parchg[i]), (spin+1, k+1, parchg[j]))
				elif not in_args.same:
					tdm = from_wav.TransitionDipoleMomentBetweenDifferentWAVECAR(other=to_wav,ks_i= [spin+1, k+1, parchg[i]], ks_j= [spin+1, k+1, parchg[j]])
				output[spin,k,i,j,0]=tdm[0] # energy
				output[spin,k,i,j,1]=tdm[1] # wave func overlap
				output[spin,k,i,j,2]=tdm[2][0] # mu_x in debye
				output[spin,k,i,j,3]=tdm[2][1] # mu_y in debye
				output[spin,k,i,j,4]=tdm[2][2] # mu_z in debye

np.savez(in_args.out_name, bands=parchg, tdm=output)
