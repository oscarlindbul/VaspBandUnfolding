import numpy as np
import xml.etree.ElementTree as ET
from argparse import ArgumentParser as Parser

parser = Parser("Reads the dielectric tensor and oscillator strengths from a BSE calculation")
parser.add_argument("xml_file", metavar="vasprun.xml", help="xml file of output")
args = parser.parse_args()

try:
    root = ET.parse(args.xml_file).getroot()
except:
    # try to fix incomplete xml-file from BSE
    with open(args.xml_file, "a") as writer:
        writer.write(r"</modeling>" + "\n")
    root = ET.parse(args.xml_file).getroot()

eps_imag_root = root.find("./dielectricfunction/imag/array")
vals = eps_imag_root.findall("./set/r")
energies = np.zeros(len(vals))
eps_imag = np.zeros((len(vals), 3, 3))
eps_real = np.zeros((len(vals), 3, 3))
for i,r in enumerate(vals):
    v = r.text.split()
    energies[i] = float(v[0])
    eps_imag[i,0,1] = float(v[1+3+0])
    eps_imag[i,1,2] = float(v[1+3+1])
    eps_imag[i,2,1] = float(v[1+3+2])
    eps_imag[i,:,:] -= np.transpose(eps_imag[i,:,:])
    for k in range(3):
        eps_imag[i,k,k] = float(v[1+k])

eps_real_root = root.find("./dielectricfunction/real/array")
vals = eps_real_root.findall("./set/r")
for i,r in enumerate(vals):
    v = r.text.split()
    eps_real[i,0,1] = float(v[1+3+0])
    eps_real[i,1,2] = float(v[1+3+1])
    eps_real[i,2,1] = float(v[1+3+2])
    eps_real[i,:,:] += np.transpose(eps_real[i,:,:])
    for k in range(3):
        eps_real[i,k,k] = float(v[1+k])

osc_vals = root.findall("./varray[@name='opticaltransitions']/v")
osc_E = np.zeros(len(osc_vals))
osc_strengths = np.zeros(len(osc_vals))
for i,osc in enumerate(osc_vals):
    v = osc.text.split()
    osc_E[i] = float(v[0])
    osc_strengths[i] = float(v[1])

np.savez("spectral.npz", epsE=energies, oscE=osc_E, eps_imag=eps_imag, eps_real=eps_real, osc=osc_strengths)
