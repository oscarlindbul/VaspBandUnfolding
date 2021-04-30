import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
import sys
import argparse as arg

import glob
import bz2
import re
import os

pi = np.pi
eps_0 = constants.epsilon_0
e = constants.e
c = constants.speed_of_light
h = constants.h
hbar = constants.hbar
m_e = constants.m_e
debye_to_Cm = 3.33564e-30
n = 2.6473 # SiC infrared refractive index 
n = 1.65 # hBN c-perp index
n = 2.10 # hBN c-parallel index

parser = arg.ArgumentParser(description="Calculates the oscillation strength for excited state transitions")

parser.add_argument("tdm_file", metavar="tdm_file", help="Path to tdm data file")
parser.add_argument("transition_from", metavar="transition_from", nargs=3, help="Ground state identifiers (spin, k-point index, band)")
parser.add_argument("transition_to", metavar="transition_to", nargs=3, help="Excited state identifiers (spin, k-point index, band)")
parser.add_argument("zpl_energy", metavar="zpl_energy", help="ZPL energy of the given transition")
parser.add_argument("-refrac_ind", dest="refrac_ind", help="Refractive index of host material (affects the lifetime), default is n=2.6473 for SiC")

input = parser.parse_args()

n = float(input.refrac_ind)

## Load data
with open(input.tdm_file, "r") as data_reader:
	data = np.load(data_reader)
	bands = data["bands"]
	tdm_data = data["tdm"]

### Find transition band indices
from_band = int(input.transition_from[2])
to_band = int(input.transition_to[2])
from_ind = list(bands).index(from_band)
to_ind = list(bands).index(to_band)

# fail checks
assert from_ind >= 0
assert to_ind >= 0
assert input.transition_from[1] == input.transition_to[1]
assert input.transition_from[0] == input.transition_to[0]

zpl_value = float(input.zpl_energy)
spin_index = int(input.transition_from[0])
k_index = int(input.transition_from[1])

## calculate relevant tdm
tdm_vals = tdm_data[spin_index, k_index, from_ind, to_ind, :] # in Debye
tdm_vals = tdm_vals[2:]*debye_to_Cm # in Cm
tdm_tot = np.sum(np.abs(tdm_vals)**2) # (Cm)^2
tdm_angle_z = np.degrees(np.arccos(np.abs(np.dot(np.array((0,0,1)), tdm_vals))/(np.linalg.norm(tdm_vals))))
tdm_vals_deb = tdm_vals / debye_to_Cm

### Find ZPL value of state
energy_diff = zpl_value * e # J
freq = energy_diff / h

osc_strength = 8*(pi**2)*freq*m_e*tdm_tot/(3*e**2*h)
einstein_coeff = 2*(2*pi)**3*n*freq**3*tdm_tot / (3*eps_0*h*(c**3))


other_einstein_coeff = 2*pi*n*freq**2*e**2*osc_strength/(eps_0*m_e*c**3)
ratio_check = einstein_coeff/other_einstein_coeff
if np.abs(ratio_check - 1) > 1e-6:
    raise Exception("Error!!!!!! Einstein ratio={}".format(ratio_check))

print("TDM calcs for spin channel {}, k-point index {}, from band {} to {}".format(spin_index, k_index, from_band, to_band))
print("Lifetime: {} s".format(1/einstein_coeff))
print("Rate of emission: {} per s".format(einstein_coeff))
print("Oscillator strength: {}".format(osc_strength))
print("Dipole Moment: ({},{},{})".format(np.real(tdm_vals_deb[0]), np.real(tdm_vals_deb[1]), np.real(tdm_vals_deb[2])))
print("Pol angle (for z-axis): {}".format(tdm_angle_z))
