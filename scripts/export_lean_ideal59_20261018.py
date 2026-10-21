#!/usr/bin/env python3
"""Export all 1890 terminal profiles as kernel-checkable Lean data.

Python only serializes saved integers. Lean recomputes every polynomial,
vertex count, index sum, rational envelope cost and both margin inequalities.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/certificates/finite_volume58_arithmetic_20261018.json'
EXPECTED = '4478785c83ea7da75752806b0cfc3f78e321b42b2bb6a267fe94442ac57dfe73'
DEST = ROOT / 'formal/CoxeterProofs/Generated/Ideal59.lean'

def ns(xs):
    return '[' + ','.join(map(str, xs)) + ']'

def high_parts(left, least=9):
    if left == 0:
        yield ()
    for r in range(least, left+1):
        for tail in high_parts(left-r, r):
            yield (r,)+tail

def render():
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == EXPECTED
    data = json.loads(SOURCE.read_text())
    for path, expected in data['sources'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == expected, path
    by_rank = {tuple(p['ranks']): p for p in data['ideal_profiles']}
    domain = [(low,)+rs if low else rs for low in range(9) for rs in high_parts(58-low)]
    assert len(domain) == len(by_rank) == 1890 and set(domain) == set(by_rank)
    profiles = [by_rank[rs] for rs in domain]
    lines = ['import CoxeterProofs.IdealProfileArithmetic', '',
        '/- Generated from the frozen dimension58 arithmetic export. -/',
        'set_option maxRecDepth 100000', 'set_option maxHeartbeats 0',
        'namespace CoxeterProofs.Generated.Ideal59',
        'open CoxeterProofs.IdealProfile']
    envelopes = [[0,1]] + [pair for _,pair in data['normalized_ideal_envelopes']]
    lines.append('def envelopes : List FractionData := [' +
        ','.join(f'⟨{a},{b}⟩' for a,b in envelopes) + ']')
    names = []
    for k in range(0,len(profiles),32):
        name = f'chunk{k//32:02}'
        names.append(name)
        lines.append(f'def {name} : List Profile := [')
        chunk = []
        for p in profiles[k:k+32]:
            a,b = p['normalized_cost']
            chunk.append(f'  ⟨{ns(p["ranks"])},{ns(p["coefficients"])},'
                f'{p["V"]},{p["G"]},⟨{a},{b}⟩⟩')
        lines.append(',\n'.join(chunk) + ']')
        lines.append(f'theorem {name}_checked : {name}.all (check envelopes) = true := by decide')
    lines.append('def profiles : List Profile := ' + ' ++ '.join(names))
    lines.extend(['theorem all_profiles_checked : profiles.all (check envelopes) = true := by',
        '  simp only [profiles, List.all_append, ' + ', '.join(n+'_checked' for n in names) +
        ', Bool.and_self]',
        'theorem profile_count : profiles.length = 1890 := by decide',
        'theorem domain_complete : profiles.map Profile.ranks = terminalDomain := by decide',
        'theorem every_profile_checked (p : Profile) (membership : p ∈ profiles) :',
        '    check envelopes p = true :=',
        '  (List.all_eq_true.mp all_profiles_checked) p membership',
        'theorem every_profile_gap (p : Profile) (membership : p ∈ profiles) :',
        '    91 * p.cost.den + 58 * p.vertices * p.cost.den +',
        '      2 * p.vertices * p.cost.num ≤ 5 * p.indexTotal * p.cost.den :=',
        '  absolute_gap_of_check envelopes p (every_profile_checked p membership)',
        'end CoxeterProofs.Generated.Ideal59', ''])
    return '\n'.join(lines)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    rendered = render()
    if args.write:
        DEST.parent.mkdir(parents=True,exist_ok=True)
        DEST.write_text(rendered)
    else:
        assert DEST.read_text() == rendered
    print('PASS export replay: 1890 terminal ideal profiles; geometric bridge pending')
