#!/usr/bin/env python3
"""Fetch one official MIT-licensed LigandMPNN checkpoint and verify observed hash."""
import argparse
import hashlib
from pathlib import Path
import os
import tempfile
import urllib.request

URL = 'https://files.ipd.uw.edu/pub/ligandmpnn/ligandmpnn_v_32_010_25.pt'
EXPECTED = '161cd264061fda9680cbb940255522ae42f2966c552d045d87913d9452a80970'

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='external/LigandMPNN/model_params/ligandmpnn_v_32_010_25.pt')
    args = parser.parse_args(); output = Path(args.output)
    if output.exists():
        if digest(output) != EXPECTED:
            raise ValueError('Existing checkpoint hash differs; preserve it and investigate')
        print('Verified existing checkpoint:', output); return
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix='checkpoint-', dir=output.parent)
    os.close(fd); temporary = Path(name)
    try:
        with urllib.request.urlopen(URL, timeout=60) as response, temporary.open('wb') as handle:
            count = 0
            while chunk := response.read(1024 * 1024):
                count += len(chunk)
                if count > 20_000_000:
                    raise ValueError('Download exceeds expected size budget')
                handle.write(chunk)
        if digest(temporary) != EXPECTED:
            raise ValueError('Downloaded checkpoint hash differs from recorded version')
        # Hard link provides exclusive creation; never overwrite a concurrent file.
        os.link(temporary, output)
        print('Downloaded and verified:', output)
    finally:
        temporary.unlink(missing_ok=True)

if __name__ == '__main__': main()
