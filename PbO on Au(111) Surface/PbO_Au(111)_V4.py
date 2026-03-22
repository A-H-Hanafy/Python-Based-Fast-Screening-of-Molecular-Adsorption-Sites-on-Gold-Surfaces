import numpy as np
from ase import Atoms, Atom
from ase.build import bulk, fcc111
from ase.optimize import BFGS
from ase.constraints import FixAtoms
from ase.filters import UnitCellFilter
from ase.io import write
from ase.io.trajectory import Trajectory
from chgnet.model.model import CHGNet
from chgnet.model.dynamics import CHGNetCalculator
import warnings

warnings.filterwarnings("ignore", message="Skipping unhashable information")

def get_relaxed_energy(atoms, label):
    atoms.calc = calc
    #write(f"{label}_initial.vasp", atoms, format='vasp')

    traj = Trajectory(f"/content/res/{label}.traj", 'w', atoms)

    if len(atoms) > 2:
        z_min = np.min(atoms.positions[:, 2])
        mask = [pos[2] < (z_min + 1.0) for pos in atoms.positions]
        atoms.set_constraint(FixAtoms(mask=mask))

    dyn = BFGS(atoms, logfile=None)
    dyn.attach(traj.write, interval=1)
    dyn.run(fmax=0.05)

    #write(f"{label}_final.vasp", atoms, format='vasp')
    return atoms.get_potential_energy()

model = CHGNet.load()
calc = CHGNetCalculator(model=model)

#gold_bulk = bulk('Au', 'fcc', a=4.0, cubic=True)
#gold_bulk.calc = calc
#ucf = UnitCellFilter(gold_bulk)
#dyn = BFGS(ucf, logfile=None)
#dyn.run(fmax=0.01)

#a_opt = gold_bulk.get_cell()[0, 0]
#print(f'Optimized Lattice Constant a = {a_opt:.4f}')

slab = Atoms(
    symbols=['Au']*48,
    positions=[
        [2.162792418, 0.416273963, 5.000027010], [-0.720892382, 2.081170159, 7.354685128],
        [0.720950007, 1.248722043, 9.709145886], [5.046592448, 0.416273963, 5.000027010],
        [2.162907605, 2.081170159, 7.354685128], [3.604750037, 1.248722043, 9.709145886],
        [7.930392134, 0.416273963, 5.000027010], [5.046707978, 2.081170159, 7.354685128],
        [6.488550067, 1.248722043, 9.709145886], [10.814192164, 0.416273963, 5.000027010],
        [7.930508008, 2.081170159, 7.354685128], [9.372350097, 1.248722043, 9.709145886],
        [0.720892425, 2.913718011, 5.000027010], [-2.162792397, 4.578614244, 7.354685128],
        [-0.720950007, 3.746166128, 9.709145886], [3.604692455, 2.913718011, 5.000027010],
        [0.721007590, 4.578614244, 7.354685128], [2.162850022, 3.746166128, 9.709145886],
        [6.488492141, 2.913718011, 5.000027010], [3.604807963, 4.578614244, 7.354685128],
        [5.046650052, 3.746166128, 9.709145886], [9.372292170, 2.913718011, 5.000027010],
        [6.488607993, 4.578614244, 7.354685128], [7.930450082, 3.746166128, 9.709145886],
        [-0.721007762, 5.411162394, 5.000027010], [-3.604692240, 7.076058032, 7.354685128],
        [-2.162850022, 6.243610213, 9.709145886], [2.162792268, 5.411162394, 5.000027010],
        [-0.720892253, 7.076058032, 7.354685128], [0.720950007, 6.243610213, 9.709145886],
        [5.046591954, 5.411162394, 5.000027010], [2.162908120, 7.076058032, 7.354685128],
        [3.604750037, 6.243610213, 9.709145886], [7.930391984, 5.411162394, 5.000027010],
        [5.046708150, 7.076058032, 7.354685128], [6.488550067, 6.243610213, 9.709145886],
        [-2.162907777, 7.908606479, 5.000027010], [-5.046592255, 9.573502117, 7.354685128],
        [-3.604750037, 8.741054298, 9.709145886], [0.720892253, 7.908606479, 5.000027010],
        [-2.162792268, 9.573502117, 7.354685128], [-0.720950007, 8.741054298, 9.709145886],
        [3.604691939, 7.908606479, 5.000027010], [0.721008106, 9.573502117, 7.354685128],
        [2.162850022, 8.741054298, 9.709145886], [6.488491969, 7.908606479, 5.000027010],
        [3.604808135, 9.573502117, 7.354685128], [5.046650052, 8.741054298, 9.709145886]
    ],
    cell=[[11.5352, 0, 0], [-5.7676, 9.98977, 0], [0, 0, 19.7092]],
    pbc=[True, True, True]
)

e_slab = get_relaxed_energy(slab.copy(), "pure_slab")

bond_length = 2.3
pbo_mol = Atoms('PbO', positions=[[0, 0, 0], [bond_length, 0, 0]])
pbo_mol.center(vacuum=10.0)
e_pbo_free = get_relaxed_energy(pbo_mol, "isolated_pbo")

z_levels = sorted(list(set(np.round(slab.positions[:, 2], 2))))
z_surf = z_levels[-1]
z_sub = z_levels[-2]

surf_atoms = slab[slab.positions[:, 2] > (z_surf - 0.5)]
#center_xy = np.diag(slab.get_cell())[:2] / 2
#dist_to_center = np.linalg.norm(surf_atoms.positions[:, :2] - center_xy, axis=1)
#ref_idx = np.argmin(dist_to_center)


all_surf_indices = [i for i, atom in enumerate(slab) if atom.position[2] > (z_surf - 0.5)]

ref_idx = all_surf_indices[4]

ref_pos = slab.positions[ref_idx]

dists = []
for idx in all_surf_indices:
    d = np.linalg.norm(slab.positions[idx][:2] - ref_pos[:2])
    dists.append((idx, d))

sorted_neighbors = sorted(dists, key=lambda x: x[1])
n1_idx = sorted_neighbors[0][0] 
n2_idx = sorted_neighbors[1][0] 
n3_idx = sorted_neighbors[2][0] 

#ref_idx = surf_atoms[4].index # For Example
#ref_pos = surf_atoms.positions[ref_idx]

#dists = np.linalg.norm(surf_atoms.positions - ref_pos, axis=1)
neighbor_idx = sorted_neighbors[1][0]#np.argsort(dists)[1]
neighbor_pos = slab.positions[neighbor_idx]
bridge_site = (ref_pos[:2] + neighbor_pos[:2]) / 2

third_idx = sorted_neighbors[2][0]#np.argsort(dists)[2]
third_pos = slab.positions[third_idx]

n4_idx = sorted_neighbors[3][0]
fourth_pos = slab.positions[n4_idx]


hollow_pos = (ref_pos[:2] + neighbor_pos[:2] + third_pos[:2]) / 3
other_hollow_pos = (ref_pos[:2] + neighbor_pos[:2] + fourth_pos[:2]) / 3

is_hcp = False
for s_pos in sub_atoms.positions:
    if np.linalg.norm(s_pos[:2] - hollow_pos) < 0.6:
        is_hcp = True
        break

final_sites = {
    'ontop': ref_pos[:2],
    'bridge': bridge_site,
    'hcp' if is_hcp else 'fcc': hollow_pos,
    'fcc' if is_hcp else 'hcp': other_hollow_pos
}

print("\n" + "="*40)
print("="*40)
print(f"Reference Atom (On-top) Index: {ref_idx}")
print(f"Reference Atom (On-top) Coordinates (XY):  {ref_pos}")

print(f"Neighbor Atom 2 Index:         {neighbor_idx}")
print(f"Neighbor Atom 2 Coordinates (XY):  {neighbor_pos}")

print(f"Neighbor Atom 3 Index:         {third_idx}")
print(f"Neighbor Atom 3 Coordinates (XY):  {third_pos}")

print(f"Bridge Site Coordinates (XY):  {bridge_site}")

print(f"Hollow Site Coordinates (XY):  {hollow_pos}")
print(f"Hollow Site Type Identified:   {'HCP' if is_hcp else 'FCC'}")
print("="*40 + "\n")


angles = np.arange(0, 360, 180)
print(f"\n{'Site':<10} | {'Angle':<8} | {'E_total (eV)':<12} | {'E_ads (eV)':<10}")
print("-" * 55)

results = []
initial_height = 2.3

for site_name, (sx, sy) in final_sites.items():
    for angle in angles:
        system = slab.copy()
        mol = Atoms('PbO', positions=[[0, 0, 0], [bond_length, 0, 0]])
        mol.rotate(angle, 'z', center='COM')
        mol.positions += [sx, sy, z_surf + initial_height]
        system.extend(mol)
        system.wrap()
        system.set_pbc(True)

        label = f"PbO_{site_name}_deg{angle}"
        try:
            e_total = get_relaxed_energy(system, label)
            e_ads = e_total - (e_slab + e_pbo_free)
            print(f"{site_name:<10} | {angle:<8} | {e_total:>12.4f} | {e_ads:>10.4f}")
            results.append((e_ads, site_name, angle))
        except Exception as e:
            print(f"{site_name:<10} | {angle:<8} | FAILED")

if results:
    best = min(results)
    print("-" * 55)
    print(f"Most Stable Site: {best[1]} at {best[2]}° (E_ads: {best[0]:.4f} eV)")
