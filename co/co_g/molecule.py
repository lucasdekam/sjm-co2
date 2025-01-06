from ase.build import molecule

from gpaw import GPAW

atoms = molecule("CO", vacuum=10)
atoms.pbc = (True, True, False)
calc = GPAW(mode="fd", xc="RPBE", txt="molecule.txt")
atoms.calc = calc

# Run the calculation.
e = atoms.get_potential_energy()
print(e, " eV")
