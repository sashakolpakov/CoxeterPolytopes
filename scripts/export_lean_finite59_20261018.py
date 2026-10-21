#!/usr/bin/env python3
"""Serialize every saved residual table at its actually used rank cap.

The 5148 tables are truncated, never optimized. Their 285643 required
last-item comparisons are recomputed by Lean's kernel. Contextual domain
coverage and the geometric interpretation remain separately documented.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/certificates/finite_volume58_arithmetic_20261018.json'
EXPECTED='4478785c83ea7da75752806b0cfc3f78e321b42b2bb6a267fe94442ac57dfe73'
DEST=ROOT/'formal/CoxeterProofs/Generated/Finite59.lean'

def ns(xs):
    return '['+','.join(map(str,xs))+']'

def render():
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==EXPECTED
    data=json.loads(SOURCE.read_text())
    for p,h in data['sources'].items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    caps={(tuple(r['key']),r['fixed_rank5_count']):r['residual_rank'] for r in data['expanded_rows']}
    lines=['import CoxeterProofs.FiniteTableArithmetic','','set_option maxRecDepth 100000',
        'set_option maxHeartbeats 0','namespace CoxeterProofs.Generated.Finite59',
        'open CoxeterProofs CoxeterProofs.FiniteTable']
    names=[]
    for k in range(0,len(data['tables']),64):
        name=f'chunk{k//64:02}'
        names.append(name)
        entries=[]
        for t in data['tables'][k:k+64]:
            cap=caps[tuple(t['key']),t['fixed_rank5_count']]
            items=[f'⟨{r},{w}⟩' for r,w in zip(t['item_ranks'],t['item_weights']) if r<=cap]
            entries.append(f'  ⟨{cap},['+','.join(items)+f'],{ns(t["values"][:cap+1])}⟩')
        lines.append(f'def {name} : List Table := [\n'+',\n'.join(entries)+']')
        lines.append(f'theorem {name}_checked : {name}.all check = true := by decide')
    lines.extend(['def tables : List Table := '+' ++ '.join(names),
        'theorem all_tables_checked : tables.all check = true := by',
        '  simp only [tables, List.all_append, '+', '.join(n+'_checked' for n in names)+', Bool.and_self]',
        'theorem table_count : tables.length = 5148 := by decide',
        'theorem every_table_upper (t : Table) (member : t ∈ tables)',
        '    (xs : List ChargedItem) (membership : ∀ a ∈ xs, a ∈ t.allowed)',
        '    (capacity : totalRank xs ≤ t.cap) :',
        '    totalCharge xs ≤ tableValue t.values (totalRank xs) :=',
        '  upper_of_check t ((List.all_eq_true.mp all_tables_checked) t member)',
        '    xs membership capacity',
        'end CoxeterProofs.Generated.Finite59',''])
    return '\n'.join(lines)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    rendered=render()
    if args.write:DEST.write_text(rendered)
    else:assert DEST.read_text()==rendered
    print('PASS export replay: 5148 used-cap residual tables; contextual/geometric bridges pending')
