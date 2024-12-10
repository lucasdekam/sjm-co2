from ase.io import read
from ase.units import Pascal, m
import numpy as np
from gpaw.utilities import h2gpts
from gpaw.solvation.sjm import SJM, SJMPower12Potential
from gpaw.solvation import (
    EffectivePotentialCavity,
    LinearDielectric,
    GradientSurface,
    SurfaceInteraction,
)

atoms = read("../optimized_surface.traj")

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
    gpts=h2gpts(0.17, atoms.get_cell(), idiv=8),
    kpts=(4, 3, 1),
    xc="RPBE",
    maxiter=1000,
    # Implicit solvent parameters.
    cavity=cavity,
    dielectric=dielectric,
    interactions=interactions,
)
atoms.calc = calc

potential_range = np.linspace(3, 5, 5)
energies = np.zeros(potential_range.shape)

for i, pot in enumerate(potential_range):
    atoms.calc.set(sj={"target_potential": pot}, txt=f"slab_co2_{pot:.2f}V.txt")
    energies[i] = atoms.get_potential_energy()
    #    atoms.calc.write(f"slab_co2_{pot:.2f}V.binary", mode="all")

np.savetxt(
    "slab_co2_energies.txt",
    np.vstack([potential_range, energies]).T,
    header="potential-vacuum / V, energy [eV]",
)
