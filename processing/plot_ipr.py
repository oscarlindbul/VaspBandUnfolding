import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse as arg

from helperfuncs import open_file_with_check

parser = arg.ArgumentParser(description="Plots the IPR plot of given defect and state")

parser.add_argument("ipr_file", metavar="ipr_file", type=str, help="name of the ipr data file")
parser.add_argument("--all", dest="plot_all", default=False, action="store_true")
parser.add_argument("-vsE", dest="versus_E", default=False, action="store_true", help="Plots instead the ipr versus Kohn-Sham energy")

input = parser.parse_args()

ipr_filepath = input.ipr_file

with open_file_with_check(ipr_filepath) as ipr_file:
	ipr_data = np.load(ipr_file)

if input.plot_all:
	for k in range(ipr_data.shape[1]):
		fig = plt.figure()
		ax = fig.gca()
		ax.set_title("k-point: " + str(k))
		for spin_channel in range(ipr_data.shape[0]):
			ipr_channel = ipr_data[spin_channel, k, :, :]
			bands = np.arange(ipr_channel.shape[0])+1
			energies = ipr_channel[:,1]
			x_vals = bands if not input.versus_E else energies
			bands = np.arange(ipr_channel.shape[0])+1
			vals = (-1)**(spin_channel)*ipr_channel[:,2]
			ax.plot(x_vals, vals, label="spin ch " + str(spin_channel))
		max_val = np.max(np.abs(vals))
		test_vals = np.abs(vals) / max_val > 1e-5
		indices = np.argwhere(test_vals==True)
		offset = 10
		x_min = indices[0] - offset
		x_max = indices[-1] + offset
		ax.legend()
		ax.set_xlabel("Bands")
		ax.set_ylabel("IPR")
		ax.set_xlim([x_min, x_max])
	plt.show()
else:
	plt.figure()
	# spins, k-points, bands, (kpath?, band, ipr)
	for spin_channel in range(ipr_data.shape[0]):
		# take spin channel and average over k-points
		ipr_channel = ipr_data[spin_channel, 0, :, :]
		bands = np.arange(ipr_channel.shape[0])+1
		energies = ipr_channel[:,1]
		x_vals = bands if not input.versus_E else energies
		vals = (-1)**(spin_channel)*ipr_channel[:,2]
		plt.plot(x_vals, vals, label="spin ch " + str(spin_channel))
	max_val = np.max(np.abs(vals))
	test_vals = np.abs(vals) / max_val > 1e-5
	indices = np.argwhere(test_vals==True)
	offset = 10
	x_min = indices[0] - offset
	x_max = indices[-1] + offset
	plt.legend()
	plt.xlabel("Bands")
	plt.ylabel("IPR")
	#plt.xlim([x_min, x_max])
	plt.show()
