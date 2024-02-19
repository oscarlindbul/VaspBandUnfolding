import sys
import time
#sys.path.append("../..") # Add path to library if not in PYTHONPATH
import VaspBandUnfolding.vaspwfc

######## Normal Usage

from VaspBandUnfolding.vaspwfc import vaspwfc as vaspwfc_normal

print("Converting wavecar with standard method")
wfc = vaspwfc_normal("WAVECAR.ground", lgamma=False, gamma_half="x") #non-compiled

start = time.perf_counter_ns()
wfc.write_gamma_wavecar(out="WAVECAR.gamma.conv")
elapsed = time.perf_counter_ns() - start
print("Time elapsed: {} s".format(elapsed/1e9))

#exit()
######## Compiled Usage

from VaspBandUnfolding.cythonize.vaspwfc import vaspwfc as vaspwfc_c

print("Converting with compiled version")
wfc = vaspwfc_c("WAVECAR.ground", lgamma=False, gamma_half="x", omp_num_threads=32) # efficient

start = time.perf_counter_ns()
wfc.write_gamma_wavecar(out="WAVECAR.gamma.conv")
elapsed = time.perf_counter_ns() - start
print("Time elapsed: {} s".format(elapsed/1e9))
