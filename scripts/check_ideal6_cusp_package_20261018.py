#!/usr/bin/env python3
"""Serial ideal-six exact replays and optional one-paper TeX batch; no Lean."""
import argparse
import fcntl
from hashlib import file_digest
import json
from pathlib import Path
import re
import sys
from run_lean_batches_20261018 import guarded_run, resources, STATE

ROOT = Path(__file__).resolve().parents[1]
STEM = 'ideal_right_angled_six_face_constraints'
CHECKS = (
    'fresh_right_angled_ideal6_carrier_propagation_20261018',
    'fresh_ideal_ra6_cusp_transport_20261018',
    'fresh_ideal_ra6_sharp_reciprocal_20261018',
)


def digest(path):
    with path.open('rb') as stream:
        return file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', action='store_true')
    args = parser.parse_args()
    before = resources()
    ledger = []
    STATE.mkdir(exist_ok=True)
    # New producers lock themselves. The inherited carrier checker is
    # protected in the parent; never double-lock a self-locking child.
    for name in CHECKS:
        source = ROOT/f'scripts/{name}.py'
        old_hash = digest(source)
        with (STATE/'exclusive.lock').open('a') as lock:
            if name == CHECKS[0]:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            resources()
            result = guarded_run([sys.executable, '-B', str(source), '--verify'],
                                 cwd=ROOT, timeout_seconds=59)
        assert digest(source) == old_hash
        output = result.pop('output').decode(errors='replace').strip()
        assert result['output_bytes'] < 65536
        print(output, flush=True)
        ledger.append(dict(checker=name, source_sha256=old_hash, **result))
    pdf_built = False
    if args.build:
        STATE.mkdir(exist_ok=True)
        with (STATE/'exclusive.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            resources()
            source = ROOT/f'manuscripts/{STEM}.tex'
            old_hash = digest(source)
            try:
                result = guarded_run(['latexmk', '-pdf', '-interaction=nonstopmode',
                    '-halt-on-error', '-outdir=manuscripts/build', str(source)],
                    cwd=ROOT, timeout_seconds=59)
                assert digest(source) == old_hash
                result.pop('output')
                log = (ROOT/f'manuscripts/build/{STEM}.log').read_text(errors='replace')
                warnings = [line for line in log.splitlines()
                            if re.search(r'Warning|Overfull|Underfull', line)]
                assert not warnings, warnings[:10]
                pdf = ROOT/f'manuscripts/build/{STEM}.pdf'
                assert pdf.is_file() and pdf.stat().st_size > 0
                ledger.append(dict(manuscript=STEM, source_sha256=old_hash,
                    pdf_sha256=digest(pdf), final_pass_warnings=0, **result))
                pdf_built = True
            finally:
                for ext in ('aux', 'fls', 'fdb_latexmk', 'log', 'out', 'toc', 'synctex.gz'):
                    (ROOT/f'manuscripts/build/{STEM}.{ext}').unlink(missing_ok=True)
    record = dict(status='PASS_SERIAL_IDEAL6_CUSP_PACKAGE', checks=ledger,
                  resources=before, lean_compiled=False,
                  universal_geometry_lean_verified=False)
    (ROOT/'docs/ideal6_cusp_package_audit_20261018.json').write_text(
        json.dumps(record, indent=2)+'\n')
    print(json.dumps(dict(status=record['status'], exact_checkers=len(CHECKS),
          pdf_built=pdf_built, lean_compiled=False,
          resources=resources(start=False))), flush=True)


if __name__ == '__main__':
    main()
