import sys
sys.path.append("../../..")
#from VaspBandUnfolding.vaspwfc import vaspwfc as vaspwfc_normal
from VaspBandUnfolding.cythonize.vaspwfc import vaspwfc as vaspwfc_c

wfc = vaspwfc_c("WAVECAR.gamma.ground", lgamma=True, gamma_half="x") # efficient
#wfc = vaspwfc_normal("data/ground_gamma/WAVECAR.G_gamma.ground", lgamma=True, gamma_half="x") #non-compiled
wfc.write_std_wavecar(out="WAVECAR.ground")
