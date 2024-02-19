import sys
import time
sys.path.append("../../..") # Add path to library if not in PYTHONPATH

######## Normal Usage

from VaspBandUnfolding.vaspwfc import vaspwfc as vaspwfc_normal

print("Converting wavecar with standard method")
wfc = vaspwfc_normal("WAVECAR.gamma.ground", lgamma=True, gamma_half="x") #non-compiled

start = time.perf_counter_ns()
wfc.write_std_wavecar(out="WAVECAR.conv")
elapsed = time.perf_counter_ns() - start
print("Time elapsed: {} s".format(elapsed/1e9))

######## Compiled Usage

from VaspBandUnfolding.cythonize.vaspwfc import vaspwfc as vaspwfc_c

print("Converting with compiled version")
wfc = vaspwfc_c("WAVECAR.gamma.ground", lgamma=True, gamma_half="x") # efficient

start = time.perf_counter_ns()
wfc.write_std_wavecar(out="WAVECAR.conv")
elapsed = time.perf_counter_ns() - start
print("Time elapsed: {} s".format(elapsed/1e9))
