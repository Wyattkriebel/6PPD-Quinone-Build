#!/usr/bin/env python3
"""Prepare stereospecific molecules; SDF is chemical authority, PDB is coordinates only."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys


def prepare_record(record, destination, seed=17, count=8):
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem, rdMolDescriptors
    ligand_id = record['ligand_id']
    if not re.fullmatch(r'[A-Za-z0-9_-]+', ligand_id):
        raise ValueError('Unsafe ligand_id')
    molecule = Chem.MolFromSmiles(record['smiles'])
    if molecule is None:
        raise ValueError('Invalid SMILES')
    if len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError('Only single connected molecules supported')
    if any(label == '?' for _, label in Chem.FindMolChiralCenters(molecule, includeUnassigned=True)):
        raise ValueError('Specify every tetrahedral stereocenter; prepare isomers separately')
    formula = rdMolDescriptors.CalcMolFormula(molecule)
    inchikey = Chem.MolToInchiKey(molecule)
    for key, observed in [('expected_formula', formula), ('expected_inchikey', inchikey)]:
        if record.get(key) and record[key] != observed:
            raise ValueError(f'{key} mismatch: {observed}')
    canonical = Chem.MolToSmiles(molecule, isomericSmiles=True)
    molecule = Chem.AddHs(molecule)
    atom_table = []
    for index, atom in enumerate(molecule.GetAtoms()):
        name = f'{atom.GetSymbol()}{index + 1}'
        atom.SetAtomMapNum(index + 1)
        info = Chem.AtomPDBResidueInfo()
        info.SetName(name.rjust(4)); info.SetResidueName('LIG')
        info.SetResidueNumber(1); info.SetChainId('L'); info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
        atom_table.append({'index': index, 'map_number': index + 1, 'name': name,
                           'element': atom.GetSymbol(), 'formal_charge': atom.GetFormalCharge()})
    params = AllChem.ETKDGv3(); params.randomSeed = seed; params.numThreads = 1
    params.pruneRmsThresh = 0.5; params.enforceChirality = True
    ids = list(AllChem.EmbedMultipleConfs(molecule, numConfs=count, params=params))
    if not ids or not AllChem.MMFFHasAllMoleculeParams(molecule):
        raise ValueError('Embedding failed or MMFF parameters incomplete')
    optimizations = AllChem.MMFFOptimizeMoleculeConfs(molecule, numThreads=1, maxIters=2000)
    accepted = [(conf_id, energy) for conf_id, (status, energy) in zip(ids, optimizations) if status == 0]
    if not accepted:
        raise ValueError('No MMFF-converged conformer')
    out = Path(destination) / ligand_id
    out.mkdir(parents=True, exist_ok=False)
    energies = []
    with Chem.SDWriter(str(out / 'conformers.sdf')) as writer:
        for conf_id, energy in accepted:
            molecule.SetProp('_Name', f'{ligand_id}_c{conf_id:03d}')
            molecule.SetProp('atom_names_json', json.dumps(atom_table))
            molecule.SetProp('mmff_energy_kcal_mol', str(energy))
            writer.write(molecule, confId=conf_id)
            Chem.MolToPDBFile(molecule, str(out / f'c{conf_id:03d}.pdb'), confId=conf_id)
            pdb_check = Chem.MolFromPDBFile(str(out / f'c{conf_id:03d}.pdb'), sanitize=False, removeHs=False)
            if pdb_check is None or [a.GetPDBResidueInfo().GetName().strip() for a in pdb_check.GetAtoms()] != [a['name'] for a in atom_table]:
                raise ValueError('PDB atom name/order roundtrip failed')
            energies.append({'conformer_id':conf_id, 'mmff_energy_kcal_mol':energy})
    reread = list(Chem.SDMolSupplier(str(out / 'conformers.sdf'), removeHs=False))
    if len(reread) != len(accepted) or any(m is None for m in reread):
        raise ValueError('SDF roundtrip failed')
    for mol in reread:
        if [a.GetAtomMapNum() for a in mol.GetAtoms()] != list(range(1, mol.GetNumAtoms()+1)):
            raise ValueError('Atom map roundtrip failed')
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        if Chem.MolToInchiKey(Chem.RemoveHs(mol)) != inchikey:
            raise ValueError('Chemical/stereochemical identity changed during roundtrip')
    metadata = dict(record, formula=formula, inchikey=inchikey, canonical_smiles=canonical,
        atom_map=atom_table, seed=seed, rdkit_version=rdBase.rdkitVersion,
        requested_conformers=count, embedded_conformers=len(ids), converged_conformers=len(accepted),
        optimization_results=[{'conformer_id':int(cid), 'status':int(status), 'energy_kcal_mol':float(energy)} for cid, (status, energy) in zip(ids, optimizations)],
        conformers=energies, parameterization_status='NOT_MD_PARAMETERIZED',
        energy_interpretation='MMFF gas-phase energies; not aqueous populations or binding energies',
        sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir())})
    (out / 'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--gate', type=Path)
    parser.add_argument('--seed', type=int, default=17)
    parser.add_argument('--conformers', type=int, default=8)
    args = parser.parse_args()
    if not 1 <= args.conformers <= 100 or not 0 <= args.seed < 2**31:
        parser.error('Conformers must be 1–100; seed must fit signed 32-bit integer')
    manifest = json.loads(args.input.read_text())
    # All nonfixture molecules require the Phase I gate, regardless of caller labels.
    fixture_allowlist = {'CCO', 'C[C@H](O)C(=O)O', 'C[C@@H](O)C(=O)O'}
    fixture = manifest.get('purpose') == 'test_fixture' and all(r['smiles'] in fixture_allowlist for r in manifest['records'])
    if not fixture:
        if args.gate is None:
            parser.error('Phase I gate required for research ligand preparation')
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
        from pipeline import validate_gate
        validate_gate(args.gate)
    for record in manifest['records']:
        prepare_record(record, args.output_dir, args.seed, args.conformers)

if __name__ == '__main__':
    main()
