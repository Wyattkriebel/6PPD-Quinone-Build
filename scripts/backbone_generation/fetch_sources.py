#!/usr/bin/env python3
"""Retrieve exact upstream checkouts without installing packages or weights.

Existing mismatched checkouts are rejected, never reset. --check-only is offline.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def command(args, cwd=None):
    return subprocess.check_output(args, cwd=cwd, text=True, stderr=subprocess.STDOUT).strip()


def verify_references(manifest):
    for record in manifest['files']:
        path = ROOT / record['path']
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != record['sha256']:
            raise ValueError(f'Reference hash mismatch: {path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only', action='store_true')
    parser.add_argument('--only', choices=['foundry', 'RFdiffusion', 'LigandMPNN'])
    args = parser.parse_args()
    manifest = json.loads((ROOT / 'configs/reference_examples/source_manifest.json').read_text())
    verify_references(manifest)
    print('Reference checksums verified')
    records = manifest['repositories']
    for name, record in records.items():
        if args.only and name != args.only:
            continue
        destination = ROOT / 'external' / name
        if not destination.exists():
            if args.check_only:
                print(f'NOT INSTALLED: {name} @ {record["commit"]}')
                continue
            destination.mkdir(parents=True)
            command(['git', 'init', str(destination)])
            command(['git', 'remote', 'add', 'origin', record['url']], destination)
            command(['git', 'fetch', '--depth', '1', 'origin', record['commit']], destination)
            command(['git', 'checkout', '--detach', 'FETCH_HEAD'], destination)
        actual = command(['git', 'rev-parse', 'HEAD'], destination)
        if actual != record['commit']:
            raise RuntimeError(f'{destination}: existing HEAD differs; use a fresh directory, never silently reset')
        dirty = command(['git', 'status', '--porcelain', '--untracked-files=no'], destination)
        if dirty:
            raise RuntimeError(f'{destination}: tracked files modified; preserve changes and investigate')
        print(f'VERIFIED: {name} @ {actual}')


if __name__ == '__main__':
    main()
