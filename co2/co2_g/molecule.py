from ase.io import read

from gpaw import GPAW

atoms = read("CO2_g.traj")
atoms.pbc = (True, True, False)
calc = GPAW(mode="fd", xc="RPBE", txt="molecule.txt")
atoms.calc = calc

# Run the calculation.
e = atoms.get_potential_energy()
print(e, " eV")
