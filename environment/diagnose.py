#!/usr/bin/env python3
"""Read-only host diagnostics. Run inside WSL; never installs or changes drivers."""
import argparse
import datetime
import glob
import importlib.metadata
import json
import os
import platform
from pathlib import Path
import shutil
import subprocess
import sys

def command(argv):
    try:
        p = subprocess.run(argv, capture_output=True, text=True, timeout=20)
        return {"argv": argv, "returncode": p.returncode, "stdout": p.stdout, "stderr": p.stderr}
    except (OSError, subprocess.TimeoutExpired) as e:
        return {"argv": argv, "error": str(e)}

def audit():
    packages = {}
    for name in ['numpy','PyYAML','rdkit','torch','streamlit','openmm','snakemake','rc-foundry']:
        try: packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError: packages[name] = None
    result = {"timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
      "scope": "Only the host executing this script; not a remote user's machine",
      "os": platform.platform(), "python": sys.version, "python_executable": sys.executable,
      "wsl_detected": 'microsoft' in platform.release().lower(),
      "cpu_count": os.cpu_count(), "disk_bytes": dict(zip(['total','used','free'], shutil.disk_usage('.'))),
      "packages": packages, "gpu_devices": glob.glob('/dev/nvidia*'),
      "commands": [command(x) for x in [['lscpu'], ['free','-b'], ['nvidia-smi'], ['nvcc','--version'], ['conda','env','list'], ['docker','--version']]]}
    if packages['torch']:
        result['torch_probe'] = command([sys.executable, '-c',
          "import torch,json; d={'version':torch.__version__,'cuda_runtime':torch.version.cuda,'cuda_available':torch.cuda.is_available(),'compiled_arches':torch.cuda.get_arch_list()}; d['devices']=[{'name':torch.cuda.get_device_name(i),'capability':torch.cuda.get_device_capability(i),'total_memory':torch.cuda.get_device_properties(i).total_memory} for i in range(torch.cuda.device_count())]; print(json.dumps(d))"])
    result['warning'] = 'CUDA visibility is not validation of model kernels or memory requirements. Test each pinned model separately.'
    return result

if __name__ == '__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--out', required=True); args=parser.parse_args()
    path=Path(args.out); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as f: json.dump(audit(),f,indent=2)
    print(path)
