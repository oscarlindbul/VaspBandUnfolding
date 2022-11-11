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

if args.compiled: 
	from VaspBandUnfolding.cythonize.vaspwfc import vaspwfc
else:
	from VaspBandUnfolding.vaspwfc import vaspwfc
    

if args.mode == "std": 
	wfc = vaspwfc(args.wav, lgamma=True, gamma_half=args.axis, omp_num_threads=args.parallel) # efficient
	wfc.write_std_wavecar(out=args.output)
else:
	wfc = vaspwfc(args.wav, lgamma=False, gamma_half=args.axis, omp_num_threads=args.parallel)
	wfc.write_gamma_wavecar(out=args.output)
