from VaspBandUnfolding import vaspwfc
import numpy as np
import argparse
import warnings
warnings.simplefilter("ignore")

parser = argparse.ArgumentParser("Compare properties and calculate similarity of two Kohn-Sham wavefunctions")
parser.add_argument("wavs", metavar="wavs", type=str, nargs=2, help="WAVECARs")
parser.add_argument("wfc1", metavar="wfc1", nargs=3, type=int, help="Wavefunction 1 (spin, k-point, index)")
parser.add_argument("wfc2", metavar="wfc2", nargs=3, type=int, help="Wavefunction 2 (spin, k-point, index)")
parser.add_argument("-g1", dest="gamma_1", default=False, action="store_true", help="First WAVECAR is in gamma-format")
parser.add_argument("-g2", dest="gamma_2", default=False, action="store_true", help="Second WAVECAR is in gamma-format")
parser.add_argument("--phase", dest="show_phase", default=False, action="store_true", help="Show the difference in phase (degrees)")
args = parser.parse_args()

gammas = [args.gamma_1, args.gamma_2]
wavs = args.wavs

wfc_1 = vaspwfc.vaspwfc(wavs[0], lgamma=gammas[0])
wfc_2 = vaspwfc.vaspwfc(wavs[1], lgamma=gammas[1])

func1 = wfc_1.wfc_r(*args.wfc1)
func2 = wfc_2.wfc_r(*args.wfc2)

similarity = np.sum(np.conjugate(func1)*func2)
angle = np.angle(similarity, deg=True)
power = np.abs(similarity)

similarity = np.abs(np.sum(np.conjugate(func1)*func2))
if args.show_phase:
	print("Wavefunction similarity: {:.3f} {:2f}".format(similarity, angle))
else:
	print("Wavefunction similarity: {:.3f}".format(similarity))


