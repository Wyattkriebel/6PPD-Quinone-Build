"""Transparent, geometry-only analysis. Confidence and affinity are not inferred."""
import argparse
import json
from pathlib import Path
import numpy as np

AA = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(), 'ARNDCQEGHILKMFPSTWYV'))
RADII = {'C':1.70, 'N':1.55, 'O':1.52, 'S':1.80, 'P':1.80, 'F':1.47, 'CL':1.75, 'BR':1.85, 'I':1.98}

def read_pdb(path):
    """Read one PDB model. Blank/A altloc allowed; B-or-later locations rejected."""
    atoms, keys, models = [], set(), 0
    for line in Path(path).read_text().splitlines():
        if line.startswith('MODEL'):
            models += 1
            if models > 1:
                raise ValueError('Multiple PDB models: split explicitly before analysis')
        if line[:6].strip() not in {'ATOM','HETATM'}:
            continue
        if len(line) < 54:
            raise ValueError('Truncated PDB coordinate record')
        if line[16:17] not in {' ', '', 'A'}:
            raise ValueError('Alternate locations must be resolved explicitly')
        name = line[12:16].strip()
        element = line[76:78].strip().upper() if len(line) >= 78 else ''
        if not element:
            # Protein CA is carbon, not calcium. Unknown HETATM elements must be supplied.
            if line[:6].strip() != 'ATOM':
                raise ValueError('HETATM requires explicit element field')
            element = name.lstrip('0123456789')[0]
        key = (line[21:22], int(line[22:26]), line[26:27], name)
        if key in keys:
            raise ValueError(f'Duplicate atom identity {key}')
        keys.add(key)
        xyz = np.array([float(line[i:i+8]) for i in (30,38,46)])
        if not np.isfinite(xyz).all():
            raise ValueError('Nonfinite coordinates')
        atoms.append(dict(key=key, residue=key[:3], name=name, resname=line[17:20].strip(), element=element, xyz=xyz))
    if not atoms:
        raise ValueError('No atoms found')
    return atoms

def aligned_rmsd(mobile, reference):
    """Kabsch RMSD in angstrom; rows correspond by identity, never nearest neighbors."""
    p, q = np.asarray(mobile,float), np.asarray(reference,float)
    if p.shape != q.shape or p.ndim != 2 or p.shape[1] != 3 or len(p) < 3:
        raise ValueError('Need >=3 corresponding points, shape (N,3)')
    if not np.isfinite(p).all() or not np.isfinite(q).all():
        raise ValueError('Nonfinite coordinates')
    p, q = p-p.mean(0), q-q.mean(0)
    u, _, vt = np.linalg.svd(p.T @ q)
    correction = np.eye(3)
    correction[-1,-1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    return float(np.sqrt(np.mean(np.sum((p @ rotation-q)**2,axis=1))))

def backbone_rmsd(mobile_path, reference_path):
    a, b = ({x['key']:x['xyz'] for x in read_pdb(path) if x['resname'] in AA and x['name'] in {'N','CA','C','O'}} for path in (mobile_path, reference_path))
    if a.keys() != b.keys():
        raise ValueError('Backbone atom identities differ; provide explicit residue mapping first')
    keys = sorted(a)
    return aligned_rmsd([a[k] for k in keys], [b[k] for k in keys])

def analyze_pdb(path, ligand_resname='LIG', contact_cutoff=4.5, overlap_cutoff=0.8):
    atoms = read_pdb(path)
    protein = [x for x in atoms if x['resname'] in AA and x['element'] not in {'H','D'}]
    ligand = [x for x in atoms if x['resname']==ligand_resname and x['element'] not in {'H','D'}]
    if not protein or not ligand:
        raise ValueError('Protein and requested ligand heavy atoms are required')
    if len({x['residue'] for x in ligand}) != 1:
        raise ValueError('Exactly one ligand residue is required; select the intended copy')
    unknown = {x['element'] for x in protein+ligand} - RADII.keys()
    if unknown:
        raise ValueError(f'No supported van der Waals radius for {unknown}')
    p, l = np.array([x['xyz'] for x in protein]), np.array([x['xyz'] for x in ligand])
    distance = np.linalg.norm(p[:,None,:]-l[None,:,:],axis=2)
    overlap = np.array([RADII[x['element']] for x in protein])[:,None]+np.array([RADII[x['element']] for x in ligand])[None,:]-distance
    contacts=[]
    for i,j in np.argwhere(distance <= contact_cutoff):
        contacts.append({'protein_atom':list(protein[i]['key']), 'ligand_atom':list(ligand[j]['key']), 'distance_A':float(distance[i,j]), 'overlap_A':float(overlap[i,j]), 'clash':bool(overlap[i,j] > overlap_cutoff)})
    residues = {}
    for x in protein:
        residues.setdefault(x['residue'],x['resname'])
    sequences = {}
    for key, res in residues.items():
        sequences[key[0]] = sequences.get(key[0],'')+AA[res]
    sequence = ''.join(sequences.values())
    ca = np.array([x['xyz'] for x in protein if x['name']=='CA'])
    return {'source_structure':str(Path(path).resolve()), 'ligand_resname':ligand_resname,
            'sequences_by_chain':sequences, 'residue_count':len(residues),
            'amino_acid_fractions':{a:sequence.count(a)/len(sequence) for a in sorted(set(AA.values()))},
            'ca_radius_of_gyration_A':float(np.sqrt(np.mean(np.sum((ca-ca.mean(0))**2,axis=1)))) if len(ca) else None,
            'protein_heavy_atoms':len(protein), 'ligand_heavy_atoms':len(ligand),
            'minimum_contact_A':float(distance.min()), 'contact_pairs':len(contacts),
            'clash_pairs':int((overlap > overlap_cutoff).sum()), 'contacts':contacts,
            'parameters':{'contact_cutoff_A':contact_cutoff,'vdw_overlap_cutoff_A':overlap_cutoff},
            'plddt':None,'pocket_volume_A3':None,'ligand_burial_fraction':None,'hydrogen_bonds':None,
            'interpretation':'Geometry only; no binding or foldability claim. B factors are not interpreted as confidence.'}

def pareto_front(rows, objectives):
    """Return nondominated row indices; exclude rows missing any objective. direction=min/max."""
    if not objectives or any(v not in {'min','max'} for v in objectives.values()):
        raise ValueError('Explicit min/max objectives required')
    values=[]
    for i,row in enumerate(rows):
        if any(row.get(k) is None for k in objectives):
            continue
        vector=np.array([float(row[k])*(1 if d=='min' else -1) for k,d in objectives.items()])
        if np.isfinite(vector).all():
            values.append((i,vector))
    return [i for i,v in values if not any(np.all(w<=v) and np.any(w<v) for j,w in values if j!=i)]

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pdb'); parser.add_argument('--ligand',default='LIG'); parser.add_argument('--output',required=True)
    args=parser.parse_args()
    output=Path(args.output); output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x') as handle:
        json.dump(analyze_pdb(args.pdb,args.ligand),handle,indent=2,allow_nan=False)
