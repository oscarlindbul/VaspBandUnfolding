from argparse import ArgumentParser as Parser
import matplotlib.pyplot as plt
import numpy as np

parser = Parser("Plots a histogram of contributions to a specific BSE eigenvalue")
parser.add_argument("data_file", metavar="data file", help="File with extracted FATBAND data")
parser.add_argument("transition_ind", metavar="Eigenvalue index", type=int, default = 0, help="Index of transition (from 1)")
parser.add_argument("-bands", dest="plot_bands", default=False, action="store_true", help="Flag to plot vs bands instead of energy")
args = parser.parse_args()

print("Reading data...")
data = np.load(args.data_file, allow_pickle=True)
eh_pairs = data["eh_pairs"][args.transition_ind-1]
print("Data read!")

E_h = np.array([vals[1] for vals in eh_pairs])
E_e = np.array([vals[2] for vals in eh_pairs])
band_from = np.array([vals[4] for vals in eh_pairs])
band_to = np.array([vals[5] for vals in eh_pairs])
coupling = np.array([vals[3] for vals in eh_pairs])

if args.plot_bands:
    min_x, max_x = (np.min(band_from) - 1, np.max(band_to)+1)

    x_vals = np.arange(min_x, max_x)
    sigma = 0.005
    gauss = lambda x, c, sigma: np.exp(-(x - c)**2/(2*sigma**2))
    from_spectrum = np.zeros(x_vals.shape)
    to_spectrum = np.zeros(x_vals.shape)
    for i in range(len(eh_pairs)):
        from_spectrum += coupling[i] * gauss(x_vals, band_from[i], sigma)
        to_spectrum += coupling[i] * gauss(x_vals, band_to[i], sigma)

    plt.plot(x_vals, from_spectrum, label="Hole")
    plt.plot(x_vals, to_spectrum, label="Electron")
    plt.show()
else:
    min_x, max_x = (np.min(E_h)*0.9, np.max(E_e)*1.1)

    x_vals = np.linspace(min_x*0.9, max_x*1.1, 1000)
    sigma = 0.005
    gauss = lambda x, c, sigma: np.exp(-(x - c)**2/(2*sigma**2))
    from_spectrum = np.zeros(x_vals.shape)
    to_spectrum = np.zeros(x_vals.shape)
    for i in range(len(eh_pairs)):
        from_spectrum += coupling[i] * gauss(x_vals, E_h[i], sigma)
        to_spectrum += coupling[i] * gauss(x_vals, E_e[i], sigma)

    plt.plot(x_vals, from_spectrum, label="Hole")
    plt.plot(x_vals, to_spectrum, label="Electron")
    plt.show()


