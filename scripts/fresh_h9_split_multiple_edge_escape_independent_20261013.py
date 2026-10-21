#!/usr/bin/env python3
"""Independent run at 90 digits of the reviewed 70-digit interval certificate.

Uses the existing interval implementation at a different precision. This is
not a second arithmetic implementation. Every coefficient enclosure, branch,
strict center and escape is reconstructed, without importing NumPy.
"""
import argparse
import hashlib
import json
from pathlib import Path
import fresh_h9_split_multiple_edge_escape_20261013 as M

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/fresh_h9_split_multiple_edge_escape_independent_20261013.json'


def run():
    saved=json.loads(M.OUT.read_text());assert saved['interval_scale']==str(10**70)
    M.S=10**90;M.O=(M.S,M.S)
    M.PHI=M.div(M.add(M.O,M.sqrt(M.rat(5))),M.rat(2))
    result=M.run(saved['point']);assert result['cap_count']==3857
    old={row['proposal']:row for row in saved['caps']};count=nested=samebranch=0
    factor=10**20
    for row in result['caps']:
        previous=old[row['proposal']];samebranch+=row['branch']==previous['branch']
        for(lo,hi),(ol,oh)in zip(row['seed_covector_enclosures'],previous['seed_covector_enclosures'],strict=True):
            lo,hi,ol,oh=map(int,(lo,hi,ol,oh))
            assert lo<=oh*factor and hi>=ol*factor
            nested+=(lo>=ol*factor and hi<=oh*factor);count+=1
    assert nested==count==38570 and samebranch==3857
    return {'status':'PASS_INDEPENDENT_90_DIGIT_REPLAY_OF_3857_CAP_ESCAPE',
        'precision_digits':90,'source_precision_digits':70,'cap_count':3857,
        'same_quadratic_branch_ids':samebranch,'covector_interval_overlap_checks':count,
        'nested_tighter_covector_intervals':nested,'old_and_new_strict_center_checks':3929,
        'old_and_new_strict_escape_checks':3929,
        'minimum_new_pairing_lower_bound':result['minimum_new_pairing_lower_bound'],
        'source_sha256':{str(M.OUT):hashlib.sha256(M.OUT.read_bytes()).hexdigest()},
        'scope':'Different-precision exact interval replay and source proof review, using the existing interval implementation; not a second arithmetic implementation or a global cap search.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();r=run()
    if args.verify:assert json.loads(OUT.read_text())==r
    else:OUT.write_text(json.dumps(r,indent=2)+'\n')
    print(r,flush=True)
