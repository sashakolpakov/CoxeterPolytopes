#!/usr/bin/env python3
"""Replay the mixed-corner package serially; optional single-paper TeX build.

Uses the common resource gate and exclusive computation lock. No Lean runs.
"""
import argparse
import fcntl
from hashlib import file_digest
import json
from pathlib import Path
import re
import sys
from run_lean_batches_20261018 import guarded_run, resources, STATE

ROOT = Path(__file__).resolve().parents[1]
STEM = 'right_angled_mixed_corner_constraints'
CHECKS = (
    'check_right_angled_mixed_weighted_triangles_research_20261018',
    'check_right_angled_mixed_low_facet_rigidity_research_20261018',
    'check_right_angled_mixed_incidence_feasibility_research_20261018',
    'check_right_angled_mixed_complete_five_corner_research_20261018',
    'check_right_angled_mixed_opposition_moments_research_20261018',
)


def digest(path):
    with path.open('rb') as stream:
        return file_digest(stream, 'sha256').hexdigest()


def clean_tex():
    for ext in ('aux', 'fls', 'fdb_latexmk', 'log', 'out', 'toc', 'synctex.gz'):
        (ROOT/f'manuscripts/build/{STEM}.{ext}').unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', action='store_true')
    args = parser.parse_args()
    STATE.mkdir(exist_ok=True)
    ledger = []
    with (STATE/'exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        before = resources()
        for name in CHECKS:
            source = ROOT/f'scripts/{name}.py'
            old_hash = digest(source)
            result = guarded_run([sys.executable, '-B', str(source)],
                                 cwd=ROOT, timeout_seconds=59)
            assert digest(source) == old_hash, source
            output = result.pop('output').decode(errors='replace').strip()
            assert result['output_bytes'] < 65536, 'Unexpectedly large checker output'
            ledger.append(dict(checker=name, source_sha256=old_hash, **result))
            print(output, flush=True)
        pdf = None
        if args.build:
            source = ROOT/f'manuscripts/{STEM}.tex'
            old_hash = digest(source)
            try:
                result = guarded_run(['latexmk', '-pdf', '-interaction=nonstopmode',
                    '-halt-on-error', '-outdir=manuscripts/build', str(source)],
                    cwd=ROOT, timeout_seconds=59)
                assert digest(source) == old_hash, source
                result.pop('output')
                log = (ROOT/f'manuscripts/build/{STEM}.log').read_text(errors='replace')
                warnings = [line for line in log.splitlines()
                            if re.search(r'Warning|Overfull|Underfull', line)]
                assert not warnings, warnings[:10]
                pdf = ROOT/f'manuscripts/build/{STEM}.pdf'
                assert pdf.is_file() and pdf.stat().st_size > 0
                ledger.append(dict(manuscript=STEM, source_sha256=old_hash,
                    pdf_sha256=digest(pdf), final_pass_warnings=0, **result))
            finally:
                clean_tex()
        record = dict(status='PASS_SERIAL_MIXED_CORNER_PACKAGE', checks=ledger,
                      resources=before, lean_compiled=False,
                      universal_geometry_lean_verified=False)
        (ROOT/'docs/mixed_corner_package_audit_20261018.json').write_text(
            json.dumps(record, indent=2)+'\n')
        print(json.dumps(dict(status=record['status'], exact_checkers=len(CHECKS),
              pdf_built=pdf is not None, lean_compiled=False,
              resources=resources(start=False))), flush=True)


if __name__ == '__main__':
    main()
