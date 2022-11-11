import numpy as np
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser("Calculates the ipr of phonon modes in data file")
parser.add_argument("data_file", metavar="data_file", help="File containing phonon modes (yaml file)")
parser.add_argument("-LW", dest="locality_weight_limit", type=float, default=np.infty, help="Threshold for distance to origin in determining phonon locality weight")
args = parser.parse_args()

data = np.load(args.data_file)

#def meV2cm_inv(mev):
#	return mev/(4.13567*15.633302)
#def cm_inv2meV(cm_inv):
#	return cm_inv/(33.356*15.633302)
 
vasp=15.633302

nphonon = len(data["freqs"])
mode_iprs = np.zeros(nphonon)
mode_freqs = data["freqs"] #*4.13567*vasp # THz to meV
mode_eigs = data["eigs"]
lw_vals = np.zeros(nphonon)
lattice_points = data["atoms"]
lattice = data["lattice"]
for i in range(lattice_points.shape[0]):
	point = lattice_points[i,:]
	cartesian_point = np.zeros((1,3))
	for x, vector in zip(point, np.split(lattice, 3, axis=0)):
		if x > 0.5:
			x -= 1
		cartesian_point += x*np.reshape(vector, (1,3))
	lattice_points[i,:] = cartesian_point
defect_pos = np.array([0,0,0])
for i in range(nphonon):
	disp_4 = 0
	locality_disp = 0
	disp_total = 0
	for j, mode in enumerate(mode_eigs[i]):
		disp = np.dot(mode, mode)
		disp_4 += disp**2
		disp_total += disp
		if np.linalg.norm(lattice_points[j,:] - defect_pos) < args.locality_weight_limit:
			locality_disp += disp
	mode_iprs[i] = disp_4
	lw_vals[i] = locality_disp/disp_total

fig = plt.figure()
plt.scatter(mode_freqs, mode_iprs, s=10)
ax = plt.gca()
secax = ax.twiny()
secax.scatter(mode_freqs*521.47083116/64.65414, mode_iprs, s=10)
secax.set_xlabel("cm^{-1}")
ax.set_xlabel("meV")
ax.set_ylabel("IPR")
ax.grid(True)

#plt.figure()
#plt.scatter(mode_freqs, lw_vals, s=1)
#ax = plt.gca()
#secax = ax.twiny()
#secax.scatter(meV2cm_inv(mode_freqs), lw_vals, s=1)
#secax.set_xlabel("cm^{-1}")
#ax.set_xlabel("meV")
#ax.set_ylabel("loc_ratio")
#ax.grid(True)
#
plt.figure()
plt.scatter(mode_freqs, mode_iprs*lw_vals, s=15,c='C4')
ax = plt.gca()
secax = ax.twiny()
secax.scatter(mode_freqs*521.47083116/64.65414, mode_iprs*lw_vals, s=1)
secax.set_xlabel("cm^{-1}")
ax.set_xlabel("meV")
ax.set_ylabel("IPR*loc_ratio")
ax.grid(True)

plt.savefig('file.pdf', format='pdf')
plt.show()
