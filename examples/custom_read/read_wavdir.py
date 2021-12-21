import numpy as np
from VaspBandUnfolding.vaspwfc_custom import vaspwfc

wav = vaspwfc("WAVDIR/WAVHEAD")
phi = wav.wfc_r(ion=1,ispin=2,ikpt=1,iband=5)
wav.save2vesta(phi, poscar="POSCAR")
