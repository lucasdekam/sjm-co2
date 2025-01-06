"""
Make config and geomopt
"""

import numpy as np

from ase.build import fcc211, molecule
from ase.optimize import BFGS
from ase.constraints import FixAtoms

from ase.units import Pascal, m

from gpaw.utilities import h2gpts
from gpaw.solvation.sjm import SJM, SJMPower12Potential
from gpaw.solvation import (
    EffectivePotentialCavity,
    LinearDielectric,
    GradientSurface,
    SurfaceInteraction,
)

# Parameters for the surface
surface = fcc211("Au", size=(3, 4, 4), a=4.20, vacuum=10)

# Read the CO2 molecule and add it to the surface
co = molecule("CO")
top_gold_atom = surface[np.argmax(surface.positions[:, 2])]
co.positions += top_gold_atom.position + np.array([0, 0, 3])

# Fix the bottom layer of the surface
indices = np.argsort(surface.positions[:, 2])[:24]
surface.set_constraint(FixAtoms(indices=indices))

surface = surface + co

# Implicit solvent parameters (to SolvationGPAW).
epsinf = 78.36  # dielectric constant of water at 298 K
gamma = 18.4 * 1e-3 * Pascal * m
cavity = EffectivePotentialCavity(
    effective_potential=SJMPower12Potential(H2O_layer=False),
    temperature=298.15,  # K
    surface_calculator=GradientSurface(),
)
dielectric = LinearDielectric(epsinf=epsinf)
interactions = [SurfaceInteraction(surface_tension=gamma)]

# The calculator
calc = SJM(
    # General GPAW parameters.
    mode="fd",
    gpts=h2gpts(0.17, surface.get_cell(), idiv=8),
    kpts=(4, 3, 1),
    xc="RPBE",
    maxiter=1000,
    # Implicit solvent parameters.
    cavity=cavity,
    dielectric=dielectric,
    interactions=interactions,
)

# Attach the calculator to the surface
surface.calc = calc
surface.calc.set(sj={"target_potential": 3.00}, txt=f"optimization_3.00V.txt")

# Perform geometry optimization
opt = BFGS(surface, logfile="optimization.log", trajectory="optimization.traj")
opt.run(fmax=0.02)  # Converge forces to 0.02 eV/Å

# Save the optimized geometry
surface.write("optimized_surface.traj")
