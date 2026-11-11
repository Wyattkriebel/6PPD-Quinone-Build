"""Analyze already-unwrapped, atom-mapped NPZ trajectories; coordinates in angstrom.

This is analysis, not an MD engine. It never repairs periodic imaging silently.
NPZ requires protein [frames,atoms,3], ligand [frames,atoms,3], time_ps [frames],
and unwrapped=True. Protein atoms should be the same ordered backbone selection.
"""
import argparse
import json
from pathlib import Path
import numpy as np


def align(mobile, reference):
    """Least-squares proper rotation for row-vector coordinates, no reflection."""
    mobile = np.asarray(mobile, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if mobile.shape != reference.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError('Alignment arrays must have equal [atoms,3] shapes')
    if len(mobile) < 3 or np.linalg.matrix_rank(mobile-mobile.mean(0)) < 2:
        raise ValueError('Alignment requires at least three non-collinear atoms')
    center, target = mobile.mean(0), reference.mean(0)
    u, _, vt = np.linalg.svd((mobile-center).T @ (reference-target))
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    return rotation, center, target


def analyze(protein, ligand, time_ps, unwrapped, contact_cutoff=4.5):
    """Frame-zero backbone fit; ligand RMSD is atom-order-specific, not symmetry corrected."""
    protein, ligand, times = map(lambda a: np.asarray(a, dtype=float), (protein, ligand, time_ps))
    if not unwrapped:
        raise ValueError('Unwrap trajectory and image ligand with protein before analysis')
    for label, xyz in [('protein', protein), ('ligand', ligand)]:
        if xyz.ndim != 3 or xyz.shape[-1] != 3 or xyz.shape[1] == 0 or not np.isfinite(xyz).all():
            raise ValueError(f'{label} must be finite [frames,atoms,3] coordinates')
    if len(protein) < 2 or len(protein) != len(ligand) or times.shape != (len(protein),):
        raise ValueError('Need >=2 frames with matched coordinate/time dimensions')
    if not np.isfinite(times).all() or np.any(np.diff(times) <= 0) or contact_cutoff <= 0:
        raise ValueError('Times must increase and cutoff must be positive')
    aligned_p, aligned_l = [], []
    for p, l in zip(protein, ligand):
        rotation, center, target = align(p, protein[0])
        aligned_p.append((p-center) @ rotation+target)
        aligned_l.append((l-center) @ rotation+target)
    p, l = np.array(aligned_p), np.array(aligned_l)
    rmsd = np.sqrt(np.mean(np.sum((p-p[0])**2, axis=-1), axis=1))
    ligand_rmsd = np.sqrt(np.mean(np.sum((l-l[0])**2, axis=-1), axis=1))
    displacement = np.linalg.norm(l.mean(1)-l[0].mean(0), axis=1)
    rmsf = np.sqrt(np.mean(np.sum((p-p.mean(0))**2, axis=-1), axis=0))
    contacts = np.array([np.linalg.norm(a[:,None,:]-b[None,:,:], axis=-1) < contact_cutoff for a,b in zip(p,l)])
    # Means use frames, not time weights; irregular sampling is disclosed.
    return {'coordinate_unit':'angstrom', 'time_unit':'ps', 'reference_frame':0,
            'n_frames':len(p), 'contact_cutoff_angstrom':contact_cutoff,
            'uniform_time_sampling':bool(np.allclose(np.diff(times), np.diff(times)[0])),
            'time_ps':times.tolist(), 'backbone_rmsd':rmsd.tolist(),
            'backbone_rmsf':rmsf.tolist(), 'ligand_pose_rmsd':ligand_rmsd.tolist(),
            'ligand_centroid_displacement':displacement.tolist(),
            'contact_occupancy_by_atom_pair':contacts.mean(0).tolist(),
            'first_half_mean_ligand_rmsd':float(ligand_rmsd[:len(p)//2].mean()),
            'second_half_mean_ligand_rmsd':float(ligand_rmsd[len(p)//2:].mean()),
            'limitations':['First/second-half comparison is descriptive, not proof of convergence.',
                          'Contacts are distances, not hydrogen bonds or binding affinity.',
                          'Protein selection is caller supplied; use fixed backbone or pocket atoms.',
                          'RMSD uses fixed atom mapping; symmetric ligand permutations are not optimized.']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    with np.load(args.input, allow_pickle=False) as data:
        result = analyze(data['protein'], data['ligand'], data['time_ps'], bool(data['unwrapped']))
    result['source_file'] = str(args.input.resolve())
    import hashlib
    result['source_sha256'] = hashlib.sha256(args.input.read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)

if __name__ == '__main__':
    main()
