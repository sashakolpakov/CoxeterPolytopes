#!/usr/bin/env python3
"""Check the 57/17/15 companion's bytes, links and declared replay inputs.

Metadata only: no scientific-module imports, enumeration, Lean, TeX, network
access, certificate regeneration or claim of mathematical verification.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'repro/dimension_bounds/manifest.json'
MAX_FILE = 32 * 1024 * 1024
MAX_TOTAL = 256 * 1024 * 1024


def require(condition, message):
    if not condition:
        raise ValueError(message)


def local_path(value):
    require(isinstance(value, str) and value and not Path(value).is_absolute(),
            'Expected repository-relative path')
    require('..' not in Path(value).parts, 'Parent traversal refused: ' + value)
    path = (ROOT / value).resolve()
    require(ROOT in path.parents, 'Path escapes repository: ' + value)
    return path


def digest(path):
    require(path.is_file() and path.stat().st_size <= MAX_FILE,
            'Missing or oversized file: ' + str(path))
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def literal_pins(source):
    tree = ast.parse(source.read_text())
    candidates = [node for node in tree.body if isinstance(node, ast.Assign)
                  and any(isinstance(target, ast.Name) and target.id == 'PINS'
                          for target in node.targets)]
    require(len(candidates) == 1, 'Expected one literal PINS mapping')
    pins = ast.literal_eval(candidates[0].value)
    require(isinstance(pins, dict), 'PINS must be a dictionary')
    return pins


def validate():
    data = json.loads(MANIFEST.read_text())
    require(data['schema_version'] == 1, 'Unknown companion schema')
    expected = {'FV1', 'FV2', 'FV3', 'FV4', 'C1', 'C2',
                'A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6'}
    items = data['proof_items']
    require(len(items) == len(expected) and {x['id'] for x in items} == expected,
            'Incomplete or duplicate proof catalogue')
    artifacts = data['artifacts']
    require(len(artifacts) <= 3000, 'Artifact inventory limit')
    total = 0
    for name, binding in artifacts.items():
        path = local_path(name)
        require(path.stat().st_size == binding['bytes'], 'Size mismatch: ' + name)
        total += binding['bytes']
        require(total <= MAX_TOTAL, 'Total hashing budget exceeded')
        require(digest(path) == binding['sha256'], 'Hash mismatch: ' + name)
    # This is an index-membership check, not a scan of unpublished local files.
    tracked = subprocess.run(['git', 'ls-files', '-z', '--', *artifacts],
                             cwd=ROOT, check=True, capture_output=True,
                             timeout=10).stdout.decode().split('\0')
    require(set(artifacts) <= set(tracked), 'A selected evidence file is untracked')
    supplement = local_path(data['written_proof']['source']).read_text()
    labels = set(re.findall(r'\\label\s*\{([^}]+)\}', supplement))
    article = local_path('arxiv/hyperbolic_coxeter_dimension_bounds/main.tex').read_text()
    article_labels = set(re.findall(r'\\label\s*\{([^}]+)\}', article))
    for entry in items:
        for label in entry['article_labels']:
            require(label in article_labels,
                    entry['id'] + ': missing article label ' + label)
        for label in entry['full_proof_labels']:
            require(label in labels, entry['id'] + ': missing proof label ' + label)
        for key in entry['checks']:
            require(key in data['checks'], 'Unknown checker group: ' + key)
        for module in entry['lean_modules']:
            require(module in data['lean']['modules'], 'Unknown Lean module: ' + module)
    for link in data['navigation_paths']:
        require(local_path(link).is_file(), 'Missing navigation link: ' + link)
    for key, check in data['checks'].items():
        links = check['inputs'] + check['receipts']
        if check['worker'] is not None:
            links.append(check['worker'])
        for link in links:
            require(local_path(link).is_file(), key + ': missing source/evidence link ' + link)
    replays = data['arithmetic_replays']
    require(set(replays) == {'A1', 'A2', 'A3', 'A4'},
            'Expected exactly four arithmetic replay groups')
    snapshot_files = 0
    for key, check in replays.items():
        source = local_path(check['worker'])
        require(digest(source) == check['worker_sha256'], key + ': worker changed')
        inputs = check['snapshot_inputs']
        pins = literal_pins(source)
        require(set(pins) <= set(inputs), key + ': incomplete literal pin closure')
        for logical, value in inputs.items():
            local_path(logical)
            require(value['source'] in artifacts, key + ': input missing from inventory')
            require(artifacts[value['source']]['sha256'] == value['sha256'],
                    key + ': snapshot hash disagrees with inventory')
            if logical in pins:
                require(pins[logical] == value['sha256'], key + ': historical pin mismatch')
        command = check['historical_command']
        part_count = 2 if key == 'A2' else 3
        require(len(command) == 6 + 4 * part_count and
                command[:4] == [check['worker'], '--expected-source-sha256',
                                check['worker_sha256'], 'combine'] and
                command[-2] == '--output', key + ': unexpected combine command')
        parts = []
        for i in range(4, len(command) - 2, 4):
            require(command[i] == '--part' and command[i + 2] == '--part-sha256',
                    key + ': unexpected certificate arguments')
            part = command[i + 1]
            require(part in inputs, key + ': missing certificate part')
            require(inputs[part]['sha256'] == command[i + 3],
                    key + ': certificate argument hash differs from snapshot binding')
            parts.append(part)
        require(len(set(parts)) == part_count, key + ': duplicate certificate part')
        snapshot_files += len(inputs) + 1
    for module, entry in data['lean']['modules'].items():
        require(entry['source'] in artifacts, module + ': source is unbound')
        for name in entry.get('separate_execution_records', []):
            require(name in artifacts, module + ': execution record is unbound')
    return {'status': 'PASS_COMPANION_METADATA_ONLY',
            'manifest_sha256': digest(MANIFEST),
            'validator_sha256': digest(Path(__file__)),
            'proof_items': len(items), 'bound_public_files': len(artifacts),
            'bound_public_bytes': total, 'arithmetic_snapshot_groups': len(replays),
            'arithmetic_snapshot_files_with_repetition': snapshot_files,
            'lean_module_entries': len(data['lean']['modules']),
            'scientific_workers_executed': 0, 'lean_compilations': 0,
            'new_axiom_audits': 0, 'mathematical_verification': False,
            'full_transitive_regeneration': False,
            'scope': 'Selected tracked bytes, proof labels, links, and literal arithmetic input closure only.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', help='Optional NEW metadata receipt under results/')
    args = parser.parse_args()
    report = validate()
    text = json.dumps(report, indent=2) + '\n'
    if args.output:
        path = local_path(args.output)
        require(path.parent == ROOT / 'results', 'Receipt must be directly under results/')
        with path.open('x') as stream:
            stream.write(text)
    print(text, end='')


if __name__ == '__main__':
    main()
