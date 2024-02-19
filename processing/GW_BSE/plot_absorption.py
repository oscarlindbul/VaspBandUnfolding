import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as constants
from argparse import ArgumentParser as Parser

parser = Parser("Plots the imaginary part of the dielectric function for optical spectrum")
parser.add_argument("data_file", metavar="data file", help="npz file containing extracted spectral data")
parser.add_argument("-lims", dest="interval", nargs=2, default=[0,3.2], type=float, help="Limits of the plot")
args = parser.parse_args()


optic_data = np.load(args.data_file)

E_abs = optic_data["epsE"]
E_osc = optic_data["oscE"]
eps = optic_data["eps_imag"]
osc = optic_data["osc"]


eps_abs = np.sqrt(np.sum(eps**2, axis=(1,2)))
eps_perp = np.sqrt(sum([ eps[:, i, i]**2 for i in range(0,2)])) # x and y
eps_par = np.sqrt(eps[:, 2,2]**2)

with open("abs_data.txt", "w") as writer:
    for i in range(eps_abs.shape[0]):
        writer.write("{} {} {} {}\n".format(E_abs[i], eps_abs[i], eps_perp[i], eps_par[i]))


plt.plot(E_abs, eps_abs, label="total")
plt.plot(E_abs, eps_perp, label="E \perp c")
plt.plot(E_abs, eps_par, label="E || c")
plt.xlabel("Energy (eV)")
plt.ylabel("Im\epsilon")
plt.xlim(args.lims)
plt.legend()
plt.show()
