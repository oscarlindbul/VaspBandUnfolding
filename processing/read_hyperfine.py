import numpy as np
import argparse

parser = argparse.ArgumentParser("Reads hyperfine parameters and vectors of desired ions")
parser.add_argument("hyp_file", metavar="hyp_file", help="Hyperfine data file")
parser.add_argument("ions", metavar="ions", type=int, nargs="+", help="Which ions to read (start from 1)")
parser.add_argument("--xyz", dest="xyz", default=False, action="store_true", help="Print the hyperfine matrix in xyz basis")
args = parser.parse_args()

# data in format (ions, (0 if diagonalized matrix, 1 if xyz matrix, 2 if eigenvectors columnwise), 3x3-matrix) 
data = np.load(args.hyp_file)

if args.xyz:
    mat = data[args.ions[0]-1,1,:,:]
    print("{} {} {} {} {} {}".format(mat[0,0], mat[1,1], mat[2,2], mat[0,1], mat[0,2], mat[1,2]))
else:
    for ion in args.ions:
        print("For ion {}:".format(ion))
        eigs = np.diag(data[ion-1,0,:,:])
        for label, eig in zip(["A1","A2","A3"], eigs):
            print("{}: {:.3g} MHz".format(label, eig))
        print("Vectors:")
        print(data[ion-1,2,:,:])
        print("")


