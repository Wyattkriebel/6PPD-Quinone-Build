#!/usr/bin/env python3
"""Small, fail-closed workflow runner; commands are argv lists, never shell strings."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import time
import yaml

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('x') as handle:
        json.dump(data, handle, indent=2)


def audit():
    """Report this execution host, not an inaccessible user's computer."""
    info = {'scope': 'execution_host_only', 'platform': platform.platform(),
            'python': sys.version, 'executable': sys.executable, 'cpu_count': os.cpu_count(),
            'disk': dict(zip(['total','used','free'], shutil.disk_usage(ROOT))),
            'tools': {name: shutil.which(name) for name in ['nvidia-smi','conda','docker','snakemake','git']}}
    if Path('/proc/meminfo').exists():
        info['meminfo'] = Path('/proc/meminfo').read_text()
    if shutil.which('nvidia-smi'):
        p = subprocess.run(['nvidia-smi', '--query-gpu=name,memory.total,driver_version,compute_cap', '--format=csv'], capture_output=True, text=True, timeout=15)
        info['gpu'] = {'returncode': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr}
    else:
        info['gpu'] = {'status': 'nvidia-smi_unavailable', 'cuda_compatibility': 'unverified'}
    return info


def validate_gate(path):
    path = Path(path)
    gate = json.loads(path.read_text())
    review = gate.get('independent_review', {})
    if gate.get('phase1_complete') is not True or review.get('approved') is not True:
        raise ValueError('Phase I complete and approved independent review required')
    if not review.get('reviewer') or not review.get('producer') or review['reviewer'] == review['producer']:
        raise ValueError('Reviewer must be named and different from producer')
    needed = {'original_inference', 'ligand_reproduction', 'sequence_complex_validation'}
    validations = {item['name']: item for item in gate.get('validations', [])}
    for name in needed:
        item = validations.get(name, {})
        if item.get('status') != 'passed' or not item.get('evidence') or item.get('fixture', False):
            raise ValueError(f'Missing real passed evidence: {name}')
        for evidence in item['evidence']:
            source = (path.parent / evidence['path']).resolve()
            if not source.is_file() or digest(source) != evidence['sha256']:
                raise ValueError(f'Evidence missing or checksum changed: {source}')
    return gate


def load_config(path):
    config = yaml.safe_load(Path(path).read_text())
    if not isinstance(config, dict) or config.get('schema_version') != 1:
        raise ValueError('schema_version: 1 required')
    if config.get('purpose') not in ('test_fixture', 'reproduction', 'phase2'):
        raise ValueError('purpose must be test_fixture, reproduction, or phase2')
    if not isinstance(config.get('seed'), int):
        raise ValueError('Integer seed required')
    budget = config.get('resources', {})
    if budget.get('cloud_spend_cap_usd', 0) != 0:
        raise ValueError('This runner does not provision paid compute; cloud cap must remain zero')
    if not 0 < budget.get('timeout_seconds', 600) <= 86400:
        raise ValueError('timeout_seconds must be in (0,86400]')
    seen = set()
    for stage in config.get('stages', []):
        name = stage['id']
        if not name.replace('_','').isalnum() or name in seen:
            raise ValueError('Stage IDs must be unique alphanumeric/underscore strings')
        if not set(stage.get('depends_on', [])).issubset(seen):
            raise ValueError('Dependencies must precede the dependent stage')
        if not isinstance(stage.get('argv'), list) or not stage['argv'] or not all(isinstance(x,str) for x in stage['argv']):
            raise ValueError('argv must be a nonempty list of strings')
        for output in stage.get('outputs', []):
            if Path(output).is_absolute() or '..' in Path(output).parts:
                raise ValueError('Outputs must stay inside the stage directory')
        if not stage.get('outputs'):
            raise ValueError('Every stage must declare outputs')
        seen.add(name)
    if not seen:
        raise ValueError('No stages configured')
    return config


def plan(config, run_dir):
    values = {'root': str(ROOT), 'python': sys.executable, 'seed': str(config['seed']), 'run_dir': str(run_dir)}
    result = []
    for stage in config['stages']:
        substitutions = dict(values, output_dir=str(run_dir / stage['id']))
        item = dict(stage)
        item['argv'] = [arg.format_map(substitutions) for arg in stage['argv']]
        item['inputs'] = [str((ROOT / arg.format_map(substitutions)).resolve()) for arg in stage.get('inputs', [])]
        result.append(item)
    return result


def run(config_path, run_dir, retry_failed=False):
    config_path = Path(config_path).resolve()
    config = load_config(config_path)
    run_dir = Path(run_dir).resolve()
    if config['purpose'] == 'phase2':
        validate_gate(ROOT / config.get('phase1_gate', 'reports/phase1_gate.json'))
    stages = plan(config, run_dir)
    config_hash = digest(config_path)
    run_dir.mkdir(parents=True, exist_ok=True)
    header_path = run_dir / 'run.json'
    if header_path.exists():
        if json.loads(header_path.read_text())['config_sha256'] != config_hash or json.loads(header_path.read_text()).get('runner_sha256') != digest(__file__):
            raise ValueError('Existing run uses a different configuration or runner; choose a new run directory')
    else:
        write_json(header_path, {'config_sha256': config_hash, 'config': config, 'hardware': audit(), 'purpose': config['purpose'], 'runner_sha256': digest(__file__), 'packages': subprocess.run([sys.executable, '-m', 'pip', 'list', '--format=json'], capture_output=True, text=True, timeout=30).stdout})
    for stage in stages:
        output_dir = run_dir / stage['id']
        manifest = output_dir / '_execution.json'
        if output_dir.exists() and retry_failed:
            previous_status = json.loads(manifest.read_text()).get('status') if manifest.exists() else 'interrupted'
            if previous_status != 'success':
                archive = run_dir / 'attempts'
                archive.mkdir(exist_ok=True)
                output_dir.rename(archive / f'{stage["id"]}_{time.time_ns()}')
        inputs = {}
        for source in stage['inputs']:
            if not Path(source).is_file():
                raise FileNotFoundError(f'Blocked {stage["id"]}: missing input {source}')
            inputs[source] = digest(source)
        if manifest.exists():
            previous = json.loads(manifest.read_text())
            if previous['status'] != 'success' or previous['inputs'] != inputs or previous['argv'] != stage['argv']:
                raise ValueError('Stage changed or failed; use --retry-failed for failures or new run for changed inputs')
            if not all(Path(p).is_file() and digest(p) == sha for p,sha in previous['outputs'].items()):
                raise ValueError('Saved output missing or modified')
            print(f'Resume: verified {stage["id"]}')
            continue
        if output_dir.exists():
            raise ValueError('Incomplete stage directory retained; use --retry-failed to archive and retry')
        output_dir.mkdir()
        started = time.time()
        record = {'stage': stage['id'], 'purpose': config['purpose'], 'argv': stage['argv'], 'inputs': inputs,
                  'seed': config['seed'], 'started_unix': started, 'status': 'failed', 'outputs': {}}
        try:
            if shutil.disk_usage(run_dir).free < config.get('resources', {}).get('minimum_free_disk_gb', 1) * 1e9:
                raise RuntimeError('Insufficient free disk')
            env = dict(os.environ, OMP_NUM_THREADS=str(config.get('resources', {}).get('cpu_threads', 2)),
                       OPENBLAS_NUM_THREADS=str(config.get('resources', {}).get('cpu_threads', 2)), PYTHONHASHSEED=str(config['seed']))
            with (output_dir / 'stdout.log').open('x') as stdout, (output_dir / 'stderr.log').open('x') as stderr:
                p = subprocess.Popen(stage['argv'], cwd=ROOT, env=env, stdout=stdout, stderr=stderr, start_new_session=True)
                try:
                    record['returncode'] = p.wait(timeout=config.get('resources', {}).get('timeout_seconds', 600))
                except subprocess.TimeoutExpired:
                    os.killpg(p.pid, signal.SIGKILL)
                    p.wait()
                    raise RuntimeError('Stage exceeded configured walltime; process group terminated')
            if record['returncode'] != 0:
                raise RuntimeError('Command failed; see stderr.log')
            for output in stage['outputs']:
                matches = sorted(output_dir.glob(output))
                if not matches:
                    raise RuntimeError(f'No files match declared output: {output}')
                for source in matches:
                    if not source.is_file() or not source.stat().st_size:
                        raise RuntimeError(f'Missing or empty declared output: {source}')
                    record['outputs'][str(source)] = digest(source)
            record['status'] = 'success'
        except Exception as error:
            record['error'] = str(error)
            raise
        finally:
            record['runtime_seconds'] = time.time() - started
            write_json(manifest, record)
    return run_dir


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['audit','validate','dry-run','run','validate-gate'])
    parser.add_argument('--config', default=str(ROOT / 'configs/smoke.yaml'))
    parser.add_argument('--run-dir', default=str(ROOT / 'results/smoke'))
    parser.add_argument('--gate')
    parser.add_argument('--retry-failed', action='store_true')
    args = parser.parse_args()
    try:
        if args.action == 'audit':
            print(json.dumps(audit(), indent=2))
        elif args.action == 'validate-gate':
            validate_gate(args.gate)
            print('Evidence files and recorded review verified; scientific claims require human assessment.')
        elif args.action == 'run':
            print(run(args.config, args.run_dir, args.retry_failed))
        else:
            config = load_config(args.config)
            print(json.dumps(plan(config, Path(args.run_dir).resolve()), indent=2))
    except (ValueError, KeyError, FileNotFoundError, RuntimeError) as error:
        print(f'BLOCKED: {error}', file=sys.stderr)
        return 2
    return 0

if __name__ == '__main__':
    sys.exit(main())
