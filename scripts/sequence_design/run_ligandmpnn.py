#!/usr/bin/env python3
"""Run pinned upstream LigandMPNN without changing its files; record compatibility shim."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import runpy
import subprocess
import sys
import time


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',required=True)
    parser.add_argument('--checkpoint',required=True)
    parser.add_argument('--input',required=True)
    parser.add_argument('--output-dir',required=True)
    parser.add_argument('--seed',type=int,default=111)
    parser.add_argument('--count',type=int,default=2)
    parser.add_argument('--temperature',type=float,default=0.1)
    args=parser.parse_args()
    if not 1 <= args.count <= 8 or args.seed <= 0 or args.temperature <= 0:
        parser.error('Pilot requires 1–8 sequences, a positive nonzero seed and positive temperature')
    repo=Path(args.repo).resolve(); output=Path(args.output_dir).resolve()
    output.mkdir(parents=True,exist_ok=True)
    provenance=output/'adapter_provenance.json'
    if provenance.exists():
        raise FileExistsError('Refusing to overwrite prior inference provenance')
    input_path=Path(args.input).resolve(); checkpoint=Path(args.checkpoint).resolve()
    argv=[str(repo/'run.py'),'--model_type','ligand_mpnn','--pdb_path',str(input_path),
          '--checkpoint_ligand_mpnn',str(checkpoint),'--out_folder',str(output),
          '--seed',str(args.seed),'--temperature',str(args.temperature),'--batch_size','1',
          '--number_of_batches',str(args.count),'--ligand_mpnn_use_atom_context','1',
          '--pack_side_chains','0','--save_stats','1']
    record={'purpose':'official_ligandmpnn_cpu_functionality_test','python':sys.version,
            'platform':platform.platform(),'argv':argv,'seed':args.seed,
            'packages':{d.metadata['Name']:d.version for d in importlib.metadata.distributions()},
            'input_sha256':sha256(input_path),'checkpoint_sha256':sha256(checkpoint),
            'upstream_python_sources':{str(p.relative_to(repo)):sha256(p) for p in sorted(repo.rglob('*.py'))},
            'git_commit':subprocess.run(['git','-C',str(repo),'rev-parse','HEAD'],capture_output=True,text=True,check=True).stdout.strip(),
            'adaptations':[],'device_requested':'cpu','status':'failed'}
    os.environ['CUDA_VISIBLE_DEVICES']=''
    started=time.time()
    try:
        import numpy as np
        # Bundled OpenFold uses the historical alias; int is its exact replacement.
        if 'int' not in np.__dict__:
            np.int=int
            record['adaptations'].append('In-process numpy.int = int for bundled OpenFold; upstream files unchanged')
        sys.path.insert(0,str(repo))
        sys.argv=argv
        runpy.run_path(str(repo/'run.py'),run_name='__main__')
        record['status']='success'
    except BaseException as error:
        record['error']=repr(error)
        raise
    finally:
        record['runtime_seconds']=time.time()-started
        with provenance.open('x') as handle:
            json.dump(record,handle,indent=2)

if __name__=='__main__': main()
