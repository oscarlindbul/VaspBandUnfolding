import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse as arg

from helperfuncs import open_file_with_check

parser = arg.ArgumentParser(description="Plots the IPR plot of given defect and state")

parser.add_argument("ipr_file", metavar="ipr_file", type=str, help="name of the ipr data file")
parser.add_argument("--all", dest="plot_all", default=False, action="store_true")
parser.add_argument("-vsE", dest="versus_E", default=False, action="store_true", help="Plots instead the ipr versus Kohn-Sham energy")
parser.add_argument("-title", dest="title", default=None, help="Title of plot")

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
            ax.plot(x_vals, vals, "-o", label="spin ch " + str(spin_channel))
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
        if input.title is not None:
            ax.set_title(input.title)
    plt.show()
else:
    plt.figure()
    # spins, k-points, bands, (kpath?, band, ipr)
    for spin_channel in range(ipr_data.shape[0]):
        # take spin channel and average over k-points
        ipr_channel = np.mean(ipr_data[spin_channel, :, :, :],axis=0)
        bands = np.arange(ipr_channel.shape[0])+1
        energies = ipr_channel[:,1]
        selection = energies > 0
        band_selection = bands[selection]
        x_vals = bands[selection] if not input.versus_E else energies[selection]
        vals = (-1)**(spin_channel)*ipr_channel[:,2][selection]
        if input.versus_E:
            plt.scatter(vals, x_vals, label="spin ch " + str(spin_channel))
            plt.ylabel("Energy")
            plt.xlabel("IPR")
            max_val = np.max(np.abs(vals))*1.5
            plt.xlim([-max_val, max_val])
            # write band labels
            for i,p in enumerate(zip(vals, x_vals)):
                plt.annotate(str(band_selection[i]), p)
        else:
            plt.plot(x_vals, vals, "-o", label="spin ch " + str(spin_channel))
            plt.xlabel("Bands")
            plt.ylabel("IPR")
            max_val = np.max(np.abs(vals))*1.2
            plt.ylim([-max_val, max_val])
    plt.legend()
    if input.title is not None:
        plt.title(input.title)
    plt.show()
