import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
import re, sys
from helperfuncs import open_with_check
from argparse import ArgumentParser as Parser

parser = Parser("Script for calculating formation energies of give input files")
parser.add_argument("outcars", metavar="OUTCARs", nargs="+", type=str, help="List of outcar files for runs")
parser.add_argument("ref_outcar", metavar="Reference OUTCAR", nargs=1, type=str, help="Unmodified neutral bulk system")
gap_group = parser.add_mutually_exclusive_group("Bandgap info")
gap_manual_group = gap_group.add_argument_group("Manual setting")
gap_auto_group = gap_group.add_argument_group("Semi-automatic setting")
gap_manual_group.add_argument("-bandgap", dest="bandgap", type=float, nargs=1, default=None, help="Value of bandgap to assume for formation energy range")
gap_manual_group.add_argument("-vbm", dest="vbm", type=float, nargs = 1, default=None, help="VBM of the system")
gap_auto_group.add_argument("-eigfile", dest="eigfile", type=str, nargs=1, default=None, help="EIGENVAL file to use for bandgap estimation (such as that of reference system)")
parser.add_argument("-o", dest=output, type=str, help="Name of output data file")
parser.add_argument("-plot", dest=plot, default=False, action="store_true", help="Should we plot the result?")
args = parser.parse_args()

if sys.version_info[0] < 3:
    potential_vals = np.load("pot_dict.save", allow_pickle=True)
else:
    potential_vals = np.load("pot_dict.save")

# outcar patterns
name_pattern = re.compile("POSCAR = (.*)")
valenz_pattern = re.compile("Ionic Valenz")
zval_pattern = re.compile("ZVAL\s+=\s+(.*)")
ion_type_n_pattern = re.compile("ions per type =\s+(.*)")
ion_type_pattern = re.compile("VRHFIN\s*=\s*(.*)\s*:")
energy_pattern = re.compile("energy without entropy\s*=\s*([-e.\d]+).*")
electron_pattern = re.compile("NELECT\s*=\s*([\d.]+)")
volume_pattern = re.compile("volume of cell\s+:\s+([\d.])")
spin_comp_pattern = re.compile("spin component (\d)")
kpoint_pattern = re.compile("k-point\s+(\d+)")
band_pattern = re.compile("(\d+)\s+([-.\d]+)\s+([\d.])")

# physical parameters
madelung = 1.64132 # 4H-SiC
epsilon = 9.66 # 4H-SiC
LZ_coeff = 2/3
eps_0 = constants["epsilon_0"]
e = constants["elementary charge"]

ref_energy = 0
ref_ion_n = {}
ref_ion_types = []
with open_with_check(ref_outcar) as outcar_file:
    lines = outcar_file.readlines()
line_ind = 0
while line_ind < len(lines): 
    # check ion types
    m = ion_type_pattern.search(lines[line_ind])
    if m is not None:
        element = m.group(0)
        ref_ion_types.append(element)
        if element not in ref_ion_n:
            ref_ion_n[element] = 0
    m = ion_type_n_pattern.search(lines[line_ind])
    if m is not None:
        element_n = list(map(int, m.group(0).split()))
        for i,N in enumerate(element_n):
            ref_ion_n[ref_ion_types[i]] += N

    line_ind += 1

# bandgap estimation
bandgap = 0
VBM = 0
if args.bandgap is not None and args.vbm is not None:
    VBM = args.bandgap
    bandgap = args.bandgap
elif args.eigfile is not None:
    pass
else:
    import warnings
    warnings.warn("No info given for bandgap, will look in reference outcar")

    bandgap_info = [[],[]]
    with open_with_check(ref_outcar) as outcar_file:
        lines = outcar_file.readlines()
    line_ind = 0
    while line_ind < len(lines):
        m = spin_comp_pattern.search(lines[line_ind])
        if m is not None:
            spin_comp = int(m.group(0))
        m = kpoint_pattern.search(lines[line_ind])
        if m is not None:
            kpoint = int(m.group(0))
            if len(band_info[spin-1]) < kpoint:
                band_info[spin-1].append(0)
            line_ind += 2
            occ_E = 0
            while True:
                m = band_pattern.search(lines[line_ind])
                if m is None:
                    break
                band_id,band_E,band_occ = m.group(0,1,2)
                if band_occ > 0:
                    occ_E = band_E
                else:
                    bandgap_info[spin-1][kpoint-1] = (band_E - occ_E, occ_E)
                    break

        line_ind += 1
    n = sum([[ 1 for i in range(len(bandgap_info[j]))] for j in range(len(bandgap_info))])
    bandgap = sum(sum([[gap for gap,vbm in l] for l in bandgap_info]))
    VBM = sum(sum([[vbm for gap,vbm, in l] for l in bandgap_info]))

    
E_F = np.linspace(0, bandgap, 500)
form_E = np.array((len(E_F), len(args.outcars)))
system_names = []
for outcar_ind,outcar in enumerate(args.outcars):
    print(outcar)
    with open_with_check(outcar) as outcar_file:
        lines = outcar_file.readlines()
    line_ind = 0
    zvals = []
    ion_types = []
    ion_n = {}
    energy = 0
    volume = 0
    while line_ind < len(lines):
        # check name
        m = name_pattern.search(lines[line_ind])
        if m is not None:
            system_names.append(m.group(0))

        # check valences
        m = valenz_pattern.search(lines[line_ind])
        if m is not None:
            line_ind += 1
            m = zval_pattern.search(lines[line_ind])
            zvals = list(map(int,m.group(0).split()))
        
        # check ion types
        m = ion_type_pattern.search(lines[line_ind])
        if m is not None:
            element = m.group(0)
            ion_types.append(element)
            if element not in ion_n:
                ref_ion_n[element] = 0
        m = ion_type_n_pattern.search(lines[line_ind])
        if m is not None:
            element_n = list(map(int, m.group(0).split()))
            for i,N in enumerate(element_n):
                ion_n[ion_types[i]] += N
        
        # check energy
        m = energy_pattern.search(lines[line_ind])
        if m is not None:
            energy = int(m.group(0))

        # check electrons
        m = electron_pattern.search(lines[line_ind])
        if m is not None:
            outcar_nelect = int(m.group(0))

        # check volume
        m = volume_pattern.search(lines[line_ind])
        if m is not None:
            volume = float(m.group(0))

        line_ind += 1

    nelect = np.array(ion_n)*np.array(zvals)
    charge = nelect - outcar_nelect # effective charge
    effective_L = np.cbrt(volume) # effective volume
    
    # calculate differences in atomic makeup
    component_diffs = {}
    all_elements = set(ion_n.keys() + ref_ion_n.keys())
    for element in all_elements:
        now = ion_n[element] if element in ion_n else 0
        before = ref_ion_n[element] if element in ref_ion_n else 0
        component_diffs[element] = now - before
    pot_E = [potential_vals[element]*change for element,change in component_diffs.items()]# chemical potential energy differences

    LZ_correction = (charge*e)^2*madelung/(2*eps_0*L) * (1/e) # J to eV
    form_E[:,outcar_ind] = energy - ref_energy - pot_E + charge*(E_F + VBM) + LZ_correction

if args.plot:
    for i in form_E.shape[1]:
        plt.plot(E_F, form_E[:,i], label=system_names[i])
    plt.legend()
    plt.show()
data = np.concatenate([E_F, form_E], axis=1)
np.savetxt(args.output, data)

