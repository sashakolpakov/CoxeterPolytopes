#!/usr/bin/env python3
"""Audit a matching capped Lean object; never compile a missing object.

This child requires the unchanged outer bounded research runner.  It obtains
the shared lock itself and refuses if any dependency needs rebuilding.  A
separate successful compilation is therefore a prerequisite, not a fallback.
"""
import os
import sys

if not __debug__ or sys.flags.optimize:
    raise SystemExit('ASSERTIONS_REQUIRED: Lean audit cannot run with Python -O')
if os.environ.get('PYTHONOPTIMIZE', '') not in ('', '0'):
    raise SystemExit('ASSERTIONS_REQUIRED: remove PYTHONOPTIMIZE before verification')

import argparse
import fcntl
import json
import re
import signal

import run_lean_batches_signal_safe_20261018 as safe


def validate_axiom_output(output, commands):
    """Require one actual result per requested name, with no unlisted axiom."""
    expected = []
    for command in commands:
        match = re.fullmatch(r'#print axioms ([A-Za-z_][A-Za-z0-9_.]*)', command)
        if match is None:
            raise RuntimeError('Unsupported axiom command: ' + command)
        expected.append(match.group(1))
    if len(set(expected)) != len(expected):
        raise RuntimeError('Duplicate requested axiom declaration')
    text = output.decode('utf-8', errors='strict')
    if 'sorryAx' in text:
        raise RuntimeError('Actual axiom output contains sorryAx')
    pattern = re.compile(r"^'([^'\n]+)'\s+(?:depends on axioms:\s*\[([^\]]*)\]"
                         r'|does not depend on any axioms)[ \t]*(?=\n|$)', re.MULTILINE)
    headers = re.compile(r"^'[^'\n]+'\s+(?:depends on axioms:"
                         r'|does not depend on any axioms)', re.MULTILINE)
    matches = list(pattern.finditer(text))
    if {m.start() for m in matches} != {m.start() for m in headers.finditer(text)}:
        raise RuntimeError('Malformed actual axiom result')
    allowed = {'propext', 'Classical.choice', 'Quot.sound'}
    results = {}
    for match in matches:
        name = match.group(1)
        if name in results:
            raise RuntimeError('Duplicate actual axiom result: ' + name)
        body = (match.group(2) or '').strip()
        axioms = [value.strip() for value in body.split(',')] if body else []
        if any(not re.fullmatch(r'[A-Za-z_][A-Za-z0-9_.]*', value)
               for value in axioms) or len(set(axioms)) != len(axioms):
            raise RuntimeError('Malformed actual axiom list for ' + name)
        if set(axioms) - allowed:
            raise RuntimeError('Unapproved axioms for ' + name + ': ' + repr(axioms))
        results[name] = axioms
    resolved = []
    for query in expected:
        # Lean prints a fully qualified name even when the source requests a
        # short name inside a namespace.  Accept only a unique suffix match;
        # a fully qualified request must match exactly.
        candidates = [name for name in results if name == query or
                      (not query.startswith('CoxeterProofs.')
                       and name.endswith('.' + query))]
        if len(candidates) != 1:
            raise RuntimeError('Missing or ambiguous actual axiom result: ' + query)
        resolved.append(candidates[0])
    if len(set(resolved)) != len(resolved) or set(resolved) != set(results):
        raise RuntimeError('Actual axiom results differ from requested declarations')
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('module')
    args = parser.parse_args()
    if not re.fullmatch(r'CoxeterProofs(?:\.[A-Za-z0-9_]+)+', args.module):
        parser.error('Choose one named module')
    owner = os.environ.get('COXETER_CHECK_OWNER', '')
    if not owner.isdigit() or int(owner) <= 1 or os.getppid() != int(owner):
        raise RuntimeError('OUTER_RSS_MONITOR_REQUIRED: use run_proof_batch lean-audit')
    base = safe.base
    if (base.ROOT / 'repro/RESEARCH_PAUSED').exists():
        raise RuntimeError('Research is paused; no audit was started')
    base.STATE.mkdir(exist_ok=True)
    with (base.STATE / 'exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        base.stop_group = safe.stop_group
        for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
            signal.signal(signum, base.interrupt)
        safe.resources()
        states = {name: base.cached(name) for name in base.dependency_order(args.module)}
        if states.get(args.module) != 'KERNEL_COMPILED':
            raise RuntimeError('SEPARATE_COMPILATION_REQUIRED: no matching capped target object')
        if any(state is None for state in states.values()):
            raise RuntimeError('SEPARATE_COMPILATION_REQUIRED: dependency needs rebuilding')
        source = base.source(args.module)
        with source.open() as stream:
            printed = [line.strip() for line in stream
                       if line.lstrip().startswith('#print axioms ')]
        if not printed:
            raise RuntimeError('NO_ACTUAL_AXIOM_COMMANDS: audit would print no declarations')
        original_guarded_run = base.guarded_run
        parsed = {}

        def validated_run(*arguments, **options):
            result = original_guarded_run(*arguments, **options)
            if result['output_bytes'] != len(result['output']):
                raise RuntimeError('Actual axiom output is truncated')
            parsed.update(validate_axiom_output(result['output'], printed))
            return result

        # Validation happens before the unchanged batch publishes its audit
        # ledger.  Calling safe.guarded_run here would recurse through base.
        base.guarded_run = validated_run
        try:
            base.batch(args.module, audit=True)
            print(json.dumps(dict(status='PASS_SEPARATE_ACTUAL_AXIOMS_PARSED',
                module=args.module, compilation_performed=False,
                printed_commands=printed, actual_axioms=parsed,
                ordinary_theorem_parameters_remain_external=True,
                cache_provenance=states)), flush=True)
        finally:
            base.guarded_run = original_guarded_run
            base.clean_unused()


if __name__ == '__main__':
    try:
        main()
    except (RuntimeError, ValueError, OSError, KeyboardInterrupt) as error:
        print(str(error), file=sys.stderr)
        sys.exit(1)
