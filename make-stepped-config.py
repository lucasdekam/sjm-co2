"""
Make config and geomopt
"""

from ase.build import fcc211
from ase.optimize import BFGS
from ase.constraints import FixAtoms
from ase import Atoms
from ase.io import read
import numpy as np
from gpaw import GPAW, PW

# Parameters for the surface
surface = fcc211("Au", size=(3, 4, 3), a=4.20, vacuum=12)

# Read the CO2 molecule and add it to the surface
co2 = read("CO2.traj")
co2 = Atoms([atom for atom in co2 if atom.symbol in ["C", "O"]])
co2.positions[:, 2] += 2.5

# Fix the bottom layer of the surface
indices = np.argsort(surface.positions[:, 2])[:12]
surface.set_constraint(FixAtoms(indices=indices))

surface = surface + co2

# Define the GPAW calculator with a plane-wave basis and 4x4x1 k-points
calc = GPAW(
    mode=PW(400),  # Plane wave cutoff energy in eV
    xc="RPBE",  # Exchange-correlation functional
    kpts=(4, 3, 1),  # k-point grid
    txt="optimization.log",
    poissonsolver={"dipolelayer": "xy"},
)  # Output log file

# Attach the calculator to the surface
surface.calc = calc

# Perform geometry optimization
opt = BFGS(surface, logfile="optimization.log", trajectory="optimization.traj")
opt.run(fmax=0.02)  # Converge forces to 0.02 eV/Å

# Save the optimized geometry
surface.write("optimized_surface.traj")
