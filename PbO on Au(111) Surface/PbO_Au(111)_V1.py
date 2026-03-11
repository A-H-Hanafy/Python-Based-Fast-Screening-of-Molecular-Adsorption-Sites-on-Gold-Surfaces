import numpy as np
from ase import Atoms, Atom
from ase.build import bulk, fcc111
from ase.optimize import BFGS
from ase.constraints import FixAtoms
from ase.filters import UnitCellFilter  
from ase.io import write
from chgnet.model.model import CHGNet
from chgnet.model.dynamics import CHGNetCalculator
import warnings
warnings.filterwarnings("ignore", message="Skipping unhashable information")

def get_relaxed_energy(atoms, label):
    atoms.calc = calc
    if len(atoms) > 2:
        mask = [atom.tag == 3 for atom in atoms]
        atoms.set_constraint(FixAtoms(mask=mask))

    dyn = BFGS(atoms, logfile=None)
    dyn.run(fmax=0.05)
    write(f"{label}_relaxed.xyz", atoms)
    return atoms.get_potential_energy()

# Load Model
model = CHGNet.load()
calc = CHGNetCalculator(model=model)

atoms = bulk('Au', 'fcc', a=4.0, cubic=True)
atoms.calc = calc

# Allow the cell to relax its volume/shape
ucf = UnitCellFilter(atoms)
dyn = BFGS(ucf, logfile=None)
dyn.run(fmax=0.01)

# The lattice constant a
a_opt = atoms.get_cell()[0, 0]
print('The lattice constant a = ',a_opt)

#Setup Au(111) Slab
slab = fcc111('Au', size=(4, 4, 3), vacuum=10.0, a=a_opt)
z_surf = np.max(slab.positions[:, 2])
e_slab = get_relaxed_energy(slab.copy(), "pure_slab")
print('The Energy of Gold Surface Au(111) = ',e_slab)

#Setup PbO Molecular (Horizontal)
bond_length = 2.3
pbo_mol = Atoms('PbO', positions=[[0, 0, 0], [bond_length, 0, 0]])
pbo_mol.center(vacuum=10.0)
e_pbo_free = get_relaxed_energy(pbo_mol, "isolated_pbo")

print('The Energy of PbO (Isolated) = ',e_pbo_free)

# d is the distance between nearest neighbor Au atoms on the (111) surface
d = a_opt / np.sqrt(2)

manual_sites = {
    'ontop':  [0.0, 0.0],
    'bridge': [d / 2.0, 0.0],
    'fcc':    [d / 2.0, d / (2.0 * np.sqrt(3))],
    'hcp':    [d / 2.0, -d / (2.0 * np.sqrt(3))]
}

print(f"\n{'Site':<10} | {'Orientation':<12} | {'E_total (eV)':<10} | {'E_ads (eV)':<10}")
print("-" * 50)

results = []
initial_height = 2.3

for site_name, (sx, sy) in manual_sites.items():
    system = slab.copy()

    # Place Pb
    pb_pos = np.array([sx, sy, z_surf + initial_height])
    system.append(Atom('Pb', position=pb_pos))

    # Place O 
    o_pos = pb_pos + np.array([bond_length, 0.0, 0.0])
    system.append(Atom('O', position=o_pos))

    system.set_pbc(True)
    system.wrap()
    label = f"flat_{site_name}"

    try:
        e_total = get_relaxed_energy(system, label)
        e_ads = e_total - (e_slab + e_pbo_free)
        print(f"{site_name:<10} | {'Flat':<12} | {e_total:>10.4f} | {e_ads:>10.4f}")
        results.append((e_ads, site_name))
    except Exception as e:
        print(f"{site_name:<10} | {'Flat':<12} | FAILED: {e}")

if results:
    best_e, best_site = min(results)
    print("-" * 40)
    print(f"Most Stable: {best_site} (Flat) at {best_e:.4f} eV")
