"""Audit reference PDB atom integrity only; does not establish chemical correctness."""
import argparse
import hashlib
import json
from pathlib import Path
from metrics import read_pdb

def audit(path,resname):
    path=Path(path); atoms=read_pdb(path)
    ligand=[a for a in atoms if a['resname']==resname]
    if not ligand:raise ValueError('Requested ligand absent')
    if len({a['residue'] for a in ligand})!=1:raise ValueError('Multiple ligand copies')
    names=[a['name'] for a in ligand]
    if len(names)!=len(set(names)):raise ValueError('Nonunique atom names')
    return {'evidence_class':'reference_input_integrity_only','source_path':str(path.resolve()),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'resname':resname,'atom_count':len(ligand),'heavy_atom_count':sum(a['element'] not in {'H','D'} for a in ligand),'atom_names_unique':True,'coordinates_finite':True,'atoms':[{'name':a['name'],'element':a['element'],'xyz_A':a['xyz'].tolist(),'identity':list(a['key'])} for a in ligand],'not_validated':['bond orders','stereochemistry','protonation','generation fidelity','binding','reproduction performance']}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('pdb');p.add_argument('--ligand',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    with Path(a.output).open('x') as f:json.dump(audit(a.pdb,a.ligand),f,indent=2)
