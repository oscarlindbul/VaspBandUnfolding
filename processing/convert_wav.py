import sys
sys.path.append("..")
#from VaspBandUnfolding.vaspwfc import vaspwfc as vaspwfc_normal
from VaspBandUnfolding.cythonize.vaspwfc import vaspwfc as vaspwfc_c
import argparse

parser = argparse.ArgumentParser("Converts gamma wavecar to complex wavecar")
parser.add_argument("gam_wav", metavar="gam_wav")
parser.add_argument("-out", dest="output", default="WAVECAR.conv", help="Output file name")
parser.add_argument("-axis", dest="axis", default="x", help="Axis of which to extrapolate complex components (binary specific)")
args = parser.parse_args()

wfc = vaspwfc_c(args.gam_wav, lgamma=True, gamma_half=args.axis) # efficient
#wfc = vaspwfc_normal("data/ground_gamma/WAVECAR.G_gamma.ground", lgamma=True, gamma_half="x") #non-compiled
wfc.write_std_wavecar(out=args.output)
