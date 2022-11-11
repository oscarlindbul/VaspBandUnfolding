import sys
sys.path.append("..")
import argparse

parser = argparse.ArgumentParser("Converts gamma wavecar to complex wavecar")
parser.add_argument("wav", metavar="wav")
parser.add_argument("mode", metavar="format", choices=["std", "gam"], help="The final output format (gam or std)")
parser.add_argument("-out", dest="output", default="WAVECAR.conv", help="Output file name")
parser.add_argument("-axis", dest="axis", default="x", help="Axis of which to extrapolate complex components (binary specific)")
parser.add_argument("-p", dest="parallel", type=int, default=1, help="Number of cores available for parallel tasks")
parser.add_argument("-c", dest="compiled", default=False, action="store_true", help="Run with compiled cython binaries")
args = parser.parse_args()

if args.mode == "std":
    from VaspBandUnfolding.cythonize.vaspwfc import vaspwfc as vaspwfc_c
    wfc = vaspwfc_c(args.wav, lgamma=True, gamma_half=args.axis) # efficient
    #wfc = vaspwfc_normal("data/ground_gamma/WAVECAR.G_gamma.ground", lgamma=True, gamma_half="x") #non-compiled
    wfc.write_std_wavecar(out=args.output)
else:
    from VaspBandUnfolding.vaspwfc import vaspwfc as vaspwfc_normal
    wfc = vaspwfc_normal(args.wav, lgamma=False, gamma_half=args.axis)
    wfc.write_gamma_wavecar(out=args.output)
