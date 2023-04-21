import os, sys
if sys.version_info[0] < 3:
    print("Requires python 3, current version: {}".format(".".join(map(str,sys.version_info))))
    exit()
import numpy as np
from pathlib import Path
from shutil import copyfile
from pymatgen.core import Structure
import argparse, glob

# borrowed from nonrad
def get_cc_structures(ground, excited, displacements, remove_zero = True):
    """Generate the structures for a CC diagram.
       Parameters
       ----------
       ground : pymatgen.core.structure.Structure pymatgen structure corresponding to the ground (final) state
       excited : pymatgen.core.structure.Structure
            pymatgen structure corresponding to the excited (initial) state
       displacements : list(float)
            list of displacements to compute the perturbed structures. Note: the
            displacements are for only one potential energy surface and will be
            applied to both (e.g. displacements=np.linspace(-0.1, 0.1, 5)) will
            return 10 structures 5 of the ground state displaced at +-10%, +-5%,
            and 0% and 5 of the excited state displaced similarly)
            remove_zero : bool
            remove 0% displacement from list (default is True)
      Returns
      -------
      ground_structs = list(pymatgen.core.structure.Struture)
         a list of structures corresponding to the displaced ground state
      excited_structs = list(pymatgen.core.structure.Structure)
         a list of structures corresponding to the displaced excited state
    """
    displacements = np.array(displacements)
    
    if remove_zero:
        displacements = displacements[displacements != 0.]
    ground_structs = ground.interpolate(excited, nimages=displacements)
    excited_structs = ground.interpolate(excited, nimages=(displacements + 1.))
    return ground_structs, excited_structs

parser = argparse.ArgumentParser("Prepare runs for nonradiative capture calculations")
parser.add_argument("ground_state", metavar="ground path", help="Path to ground state run-files")
parser.add_argument("ex_state", metavar="excited path", help="Path to excited state run-files")
parser.add_argument("-N", dest="n_disloc", default=5, type=int, help="Number of perturbations to calculate")
parser.add_argument("-disloc", dest="disloc", default=0.5, type=float, help="Amount of perturbation between the two structures [0,1]")
parser.add_argument("-o", dest="output", default="interpolation", help="Parent directory of the new runs")
args = parser.parse_args()


# equilibrium structures from your first-principles calculation
ground_files = Path(args.ground_state)
CONTCAR = glob.glob(args.ground_state + "/CONTCAR*")[0]
ground_struct = Structure.from_file(CONTCAR)
excited_files = Path(args.ex_state)
CONTCAR = glob.glob(args.ex_state + "/CONTCAR*")[0]
excited_struct = Structure.from_file(CONTCAR)

# output directory that will contain the input files for the CC diagram
cc_dir = Path('.')
if not os.path.exists(str(cc_dir / args.output)):
    os.mkdir(str(cc_dir / args.output))

# displacements as a percentage, this will generate the displacements
# -50%, -37.5%, -25%, -12.5%, 0%, 12.5%, 25%, 37.5%, 50%
displacements = np.linspace(-args.disloc, args.disloc, args.n_disloc)

# note: the returned structures won't include the 0% displacement, this is intended
# it can be included by specifying remove_zero=False
ground, excited = get_cc_structures(ground_struct, excited_struct, displacements)

run_files = ['KPOINTS', 'INCAR', 'WAVECAR', "run_vasp"]
for i, struct in enumerate(ground):
    print("Copying ground files {}".format(i+1))
    working_dir = cc_dir / args.output / str(i)
    if not os.path.exists(working_dir):
        os.mkdir(str(working_dir))

    # write structure and copy necessary input files
    struct.to(filename=str(working_dir / 'POSCAR'), fmt='poscar')
    for f in run_files:
        copyfile(str(ground_files / f), str(working_dir / f))
        if f == "INCAR":
            with open(str(working_dir / f), "r") as reader:
                lines = reader.readlines()
            for i in range(len(lines)):
                if "NSW" in lines[i]:
                    lines[i] = "NSW = 0\n"
            with open(str(working_dir / f), "w") as writer:
                writer.writelines(lines)

