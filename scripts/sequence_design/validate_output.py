#!/usr/bin/env python3
"""Independent LigandMPNN backbone/FASTA interface checks, not fold validation.

PDB checks establish names/elements/coordinates only, never a chemical bond graph.
Water and H/D are excluded explicitly. Chain ordering for FASTA is sorted chain ID,
matching the pinned upstream parse_PDB implementation.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

AA = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(), 'ARNDCQEGHILKMFPSTWYV'))
BACKBONE = {'N', 'CA', 'C', 'O'}
WATER = {'HOH', 'WAT', 'DOD'}


def parse(path):
    atoms, residues, models = {}, {}, 0
    for line in Path(path).read_text().splitlines():
        if line.startswith('MODEL'):
            models += 1
            if models > 1:
                raise ValueError('Multiple models require explicit selection')
        if line[:6].strip() not in {'ATOM', 'HETATM'}:
            continue
        if line[16:17] not in {' ', 'A'}:
            raise ValueError('Resolve alternate locations before interface validation')
        key = (line[21:22], int(line[22:26]), line[26:27], line[12:16].strip())
        element = line[76:78].strip().upper()
        resname = line[17:20].strip()
        if not element:
            raise ValueError('Explicit PDB element fields required')
        if element in {'H','D'} or resname in WATER:
            continue
        if key in atoms:
            raise ValueError(f'Duplicate atom identity: {key}')
        xyz = tuple(float(line[i:i+8]) for i in (30,38,46))
        if not all(math.isfinite(v) for v in xyz):
            raise ValueError('Nonfinite coordinates')
        atoms[key] = {'xyz':xyz, 'element':element, 'resname':resname}
        if resname in AA:
            residue = key[:3]
            if residue in residues and residues[residue] != resname:
                raise ValueError('Conflicting residue names within one residue')
            residues[residue] = resname
    if not atoms or not residues:
        raise ValueError('No protein atoms found')
    # Preserve within-chain residue encounter order; upstream input contains this order.
    sequences = {chain: ''.join(AA[name] for key,name in residues.items() if key[0] == chain)
                 for chain in sorted({key[0] for key in residues})}
    backbone = {k:v for k,v in atoms.items() if v['resname'] in AA and k[3] in BACKBONE}
    for residue in residues:
        if {k[3] for k in backbone if k[:3] == residue} != BACKBONE:
            raise ValueError(f'Incomplete N/CA/C/O backbone at {residue}')
    context = {k:v for k,v in atoms.items() if v['resname'] not in AA}
    return backbone, context, sequences, residues


def read_fasta(path):
    records = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line: continue
        if line.startswith('>'):
            records.append([line[1:], ''])
        elif not records:
            raise ValueError('FASTA sequence before header')
        else:
            records[-1][1] += line
    if not records: raise ValueError('Empty FASTA')
    return records


def validate(reference, designed, fasta, design_id=1, tolerance=0.002):
    if tolerance < 0 or not math.isfinite(tolerance): raise ValueError('Invalid tolerance')
    ref_bb, ref_lig, ref_seq, ref_res = parse(reference)
    out_bb, out_lig, out_seq, out_res = parse(designed)
    if ref_res.keys() != out_res.keys(): raise ValueError('Protein chain/residue/insertion-code identities changed')
    if ref_bb.keys() != out_bb.keys(): raise ValueError('Backbone atom identities changed')
    if ref_lig.keys() != out_lig.keys(): raise ValueError('Nonwater context atom identities changed')
    maximum = {}
    for name, before, after in [('backbone',ref_bb,out_bb), ('context',ref_lig,out_lig)]:
        displacements = []
        for key, atom in before.items():
            target = after[key]
            if atom['element'] != target['element']:
                raise ValueError(f'{name} element changed at {key}')
            if name == 'context' and atom['resname'] != target['resname']:
                raise ValueError(f'Context residue identity changed at {key}')
            displacement = math.dist(atom['xyz'],target['xyz'])
            if displacement > tolerance:
                raise ValueError(f'{name} coordinate changed at {key}: {displacement:.6g} A')
            displacements.append(displacement)
        maximum[name] = max(displacements, default=0.)
    records = read_fasta(fasta)
    native = ':'.join(ref_seq.values())
    if records[0][1] != native: raise ValueError('Native FASTA differs from reference protein sequence/chain order')
    matches = [(header,sequence) for header,sequence in records[1:]
               if re.search(r'(?:^|,\s*)id='+str(design_id)+r'(?:,|$)', header)]
    if len(matches) != 1: raise ValueError('Expected exactly one matching FASTA design ID')
    if matches[0][1] != ':'.join(out_seq.values()):
        raise ValueError('Designed FASTA differs from output PDB residue identities/chain order')
    return {'status':'passed', 'scope':'LigandMPNN sequence/backbone interface; not structure prediction',
            'design_id':design_id, 'designed_sequence':matches[0][1], 'chain_sequences':out_seq,
            'protein_residues':len(out_res), 'backbone_atoms':len(out_bb), 'context_heavy_atoms':len(out_lig),
            'context_residue_names':sorted({v['resname'] for v in out_lig.values()}),
            'coordinate_tolerance_A':tolerance, 'maximum_coordinate_displacement_A':maximum,
            'residue_mapping':[{'chain':k[0], 'residue_number':k[1], 'insertion_code':k[2],
                                'input_resname':ref_res[k], 'output_resname':out_res[k]} for k in ref_res],
            'context_atom_mapping':[{'chain':k[0], 'residue_number':k[1], 'insertion_code':k[2],
                                    'atom_name':k[3], 'element':out_lig[k]['element']} for k in out_lig],
            'limitations':['PDB does not verify bond orders, chemical graph, protonation or stereochemistry.',
                           'Backbone coordinates are copied, not independently predicted.',
                           'Water and hydrogen/deuterium atoms excluded; no solvent preservation claim.',
                           'Zero context heavy atoms would test sequence-only interface, not ligand-conditioned function.'],
            'sources':{str(Path(p).resolve()):hashlib.sha256(Path(p).read_bytes()).hexdigest()
                       for p in (reference,designed,fasta)}}


def validate_stats(path, expected_residues, design_id):
    """Validate locally generated tensors using restricted deserialization."""
    import torch
    data = torch.load(path, map_location='cpu', weights_only=True)
    probabilities = data['sampling_probs']
    logs = data['log_probs']
    sequences = data['generated_sequences']
    if probabilities.ndim != 3 or probabilities.shape[1] != expected_residues:
        raise ValueError('Probability tensor dimensions inconsistent with reference')
    if (sequences.shape != probabilities.shape[:2] or logs.shape[:2] != sequences.shape
            or probabilities.shape[2] != 20 or logs.ndim != 3 or logs.shape[2] != 21):
        raise ValueError('Sequence/probability/log-probability shapes inconsistent')
    if not 1 <= design_id <= len(sequences):
        raise ValueError('Requested design absent from statistics')
    if not torch.isfinite(probabilities).all() or not torch.isfinite(logs).all():
        raise ValueError('Nonfinite model probabilities')
    if (probabilities < 0).any():
        raise ValueError('Negative sampling probabilities')
    active = data['mask'].bool() & data['chain_mask'].bool()
    sums = probabilities[:,active,:].sum(-1)
    if not torch.allclose(sums,torch.ones_like(sums),atol=1e-4):
        raise ValueError('Sampling distributions do not normalize at designed residues')
    return {'status':'passed', 'sampling_probability_shape':list(probabilities.shape),
            'log_probability_shape':list(logs.shape),
            'alphabet_note':'Sampling excludes X (20 amino acids); model log-probabilities include X (21).',
            'seed':data['seed'], 'temperature':data['temperature'],
            'source_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            'interpretation':'Sequence-model statistics; not structure confidence or binding affinity.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--designed', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--design-id',type=int,default=1)
    parser.add_argument('--stats',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=validate(args.reference,args.designed,args.fasta,args.design_id)
    if args.stats:
        result['model_statistics']=validate_stats(args.stats,result['protein_residues'],args.design_id)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as handle: json.dump(result,handle,indent=2)

if __name__ == '__main__': main()
