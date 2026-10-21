#!/usr/bin/env python3
"""Rank62 finite upper434 by fixed-vector shifts of frozen exact tables.

Standard-library-only cost/domain recount. No Bellman optimization is run.
--prepare is explicitly conditional and writes no final certificate.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_finite_rank15_rank16_rank62_budget_20261018'
OLD = 'fresh_finite_d17_cost120_rank62_budget_20261018'
R15 = 'fresh_finite_rank15_existing_half_transfer_20261018'
R16 = 'fresh_rank16_existing_payment_exchange_20261018'
SPARSE = 'fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018'
# Filled only after both authors announce their final source/data/report freeze.
FINAL_LOCAL_HASHES = {
    R15: ('5747b6dd3864e65db88e13ca58150c147dbd0a13fea6cc54b5adaf0b429359eb',
          'b9e0dc22af8c5af0d4008fdaf4e3a2d45b5763e7bfa2a7b08d0daf22b2ae50e2',
          '4fd9e9fe75061bed5c437da345322ebdf3f1dccb5417eb35e52a075f18e4a868'),
    R16: ('6ee75a84db15c5b10fd4b93ff7b23a1d99866c9712bc3496c792be01240e244e',
          '344ed205cb5c84a2bcf21178e906a3bc964bd9c5ed0b8646c44667b4401b5603',
          '2d5d5dc422017c6602d5b0fdc0e986e78f3c37e4ef6b3a53c3b47ba9881ce7f5'),
}
N = 62


def before_base(t):
    """Doubled actual same-component cost in the frozen438 rule."""
    k, r = t
    if r == 16:
        return 208 if k == 'B' else 228
    if r == 17:
        return {'A': 228, 'B': 210, 'D': 240}[k]
    if r <= 15:
        return r * (r - 1)
    large = 16*r - (66 if k == 'D' and r >= 20 else 64 if k == 'D' else 68)
    return large - 2*max(0, r - (22 + {'A': 8, 'B': 9, 'D': 10}[k]))


def actual_base(t):
    if t[1] == 15:
        return 208
    if t == ('A', 16):
        return 226
    if t == ('D', 16):
        return 225
    if t == ('A', 17):
        return 230
    return before_base(t)


def upper_base(t):
    # A16 needs226, so the D16 representative is deliberately overallowed1.
    return 226 if t == ('D', 16) else actual_base(t)


def representative(t):
    k, r = t
    if t == ('A', 17):
        return t
    if r == 6:
        return ('E', 6)
    if r == 9:
        return ('B', 9)
    return ('A' if r <= 8 else 'D', r)


# Doubled total component-pair increments, including the final D13 upgrades.
SPARSE_CROSS = {
    (('A', 4), ('A', 12)): 8,
    (('A', 4), ('A', 13)): 8,
    (('D', 5), ('A', 12)): 4,
    (('A', 2), ('A', 13)): 4,
    (('A', 2), ('A', 14)): 4,
    (('E', 6), ('A', 11)): 8,
    (('A', 2), ('B', 9)): 2,
    (('A', 2), ('D', 10)): 4,
    (('A', 3), ('D', 10)): 2,
}


def cross(a, b):
    """Unchanged current actual coupon total, doubled; E6/D11 is already4."""
    value = 0
    for small, large in ((a, b), (b, a)):
        if small[1] <= 5:
            if large == ('D', 10):
                value += 2*small[1]
            elif large[0] == 'D' and 11 <= large[1] <= 16:
                value += 3*small[1]
            elif large in (('A', 13), ('A', 15)):
                value += small[1]
            elif large == ('A', 17):
                value += 2*small[1]
        if (small, large) in (
                (('A', 1), ('D', 11)), (('A', 1), ('D', 12)),
                (('A', 2), ('B', 10)), (('A', 4), ('D', 13)),
                (('A', 1), ('A', 11)), (('A', 2), ('A', 11))):
            value += 2
        if (small, large) == (('A', 2), ('D', 11)):
            value += 4
        if (small, large) == (('E', 6), ('D', 11)):
            value += 8
        if (small, large) == (('A', 1), ('B', 10)):
            value += 1
        if small in (('A', 1), ('A', 2)) and large == ('D', 14):
            value += small[1]
        value += SPARSE_CROSS.get((small, large), 0)
    return value


def cost(ts, base):
    return sum(base(t) for t in ts) + sum(cross(a, b)
        for i, a in enumerate(ts) for b in ts[i+1:])


def type_list():
    return ([('A', r) for r in range(1, N+1)] +
            [('B', r) for r in range(2, N+1)] +
            [('D', r) for r in range(4, N+1)] +
            [('E', 6), ('E', 7), ('E', 8), ('F', 4),
             ('H', 3), ('H', 4), ('I2*', 2)])


def domains():
    rows = []
    def visit(rank, left, counts):
        if rank == 18:
            c = dict(zip(range(9, 18), counts))
            key = (2*c[10]+3*sum(c[r] for r in range(11, 17))+2*c[17],
                   c[11]+c[12], c[13], c[11], c[14], c[9], c[10])
            rows.append((tuple(counts), N-left, key))
            return
        for count in range(left//rank+1):
            visit(rank+1, left-rank*count, counts+[count])
    visit(9, N, [])
    return sorted(rows)


def arithmetic(old):
    assert old['finite_maximum'] == '438'
    assert old['E6_D11_actual_total'] == '4'
    assert old['D13_zero_increments'] == {'A2_D10': '1', 'A3middle_D10': '1/2'}
    types = type_list()
    for a in types:
        assert upper_base(representative(a)) >= actual_base(a)
        for b in types:
            assert cross(a, b) == cross(b, a)
            assert cross(representative(a), representative(b)) >= cross(a, b)
    dom = domains()
    assert len(dom) == 1227
    tables = {tuple(t['key']): t for t in old['certificate']['tables']}
    assert len(tables) == 921 and set(tables) == {r[2] for r in dom}
    prior = {tuple(r['counts']): r for r in old['current_actual_comparisons']}
    assert len(prior) == len(dom) and set(prior) == {r[0] for r in dom}
    records = []; residual_checks = 0
    for counts, used, key in dom:
        large = [('B' if r == 9 else 'A' if r == 17 else 'D', r)
                 for r, c in zip(range(9, 18), counts) for _ in range(c)]
        assert sum(r for k, r in large) == used
        assert all(cross(a, b) == 0 for i, a in enumerate(large) for b in large[i+1:])
        for r in tuple(range(1, 9)) + tuple(range(17, N+1)):
            t = representative(('A' if r <= 8 else 'D', r))
            old_item = before_base(t) + sum(cross(t, b) for b in large)
            new_item = upper_base(t) + sum(cross(t, b) for b in large)
            assert old_item == new_item
            broad, extra, a4, d11, d14, b9, d10 = key
            encoded = before_base(t) + (r*broad if r <= 5 else 0)
            encoded += (2*extra+d14 if r == 1 else
                        4*d11+2*d14+2*b9+4*d10 if r == 2 else
                        2*d10 if r == 3 else 2*a4 if r == 4 else
                        8*d11 if r == 6 else 0)
            assert old_item == encoded
            residual_checks += 1
        table = tables[key]
        assert len(table['values']) == len(table['last']) == N+1
        rem = N-used; residual = []
        while rem:
            r = table['last'][rem]
            assert r in tuple(range(1, 9))+tuple(range(17, N+1)) and r <= rem
            residual.append(representative(('A' if r <= 8 else 'D', r)))
            rem -= r
        ts = large + residual; names = sorted(k+str(r) for k, r in ts)
        prev = prior[counts]
        assert names == prev['forest']
        old_value = cost(ts, before_base)
        assert old_value == prev['actual_doubled'] == prev['upper_doubled']
        assert old_value == sum(map(before_base, large)) + table['values'][N-used]
        c15, c16, a17 = counts[6], counts[7], counts[8]
        shift = -2*c15 - 2*c16 + 2*a17
        assert sum(map(upper_base, large))-sum(map(before_base, large)) == shift
        upper = old_value + shift
        actual = cost(ts, actual_base)
        assert actual <= upper == cost(ts, upper_base)
        assert upper == sum(map(upper_base, large)) + table['values'][N-used]
        records.append(dict(counts=list(counts),key=list(key),forest=names,
            previous_doubled=old_value,doubled_shift=shift,
            new_doubled_upper=upper,new_doubled_actual=actual))
    maximum = max(r['new_doubled_upper'] for r in records)
    assert maximum == 868
    witness = [('D', 17), ('D', 45)]
    assert sum(r for k, r in witness) == N and cost(witness, actual_base) == maximum
    attaining = [r for r in records if r['new_doubled_upper'] == maximum]
    return dict(ambient_dimension=N,finite_maximum='434',finite_curvature='434',finite_slack='0',
        actual_attaining_forest=['D17','D45'],actual_attaining_component_costs=['120','314'],
        actual_costs=dict(A15='104',B15='104',D15='104',A16='113',D16='225/2',A17='115',D17='120'),
        representative_D16_upper='113',new_cross_component_coupons=[],
        complete_large_vectors=len(dom),inherited_residual_tables=len(tables),
        simultaneous_type_pairs=len(types)**2,unchanged_residual_item_checks=residual_checks,
        inherited_table_sha256=sha256(json.dumps(old['certificate']['tables'],sort_keys=True,separators=(',',':')).encode()).hexdigest(),
        fixed_vector_shift_formula='-2*c15-2*c16+2*a17 (doubled)',
        alternatives=records,upper_attaining_alternatives=attaining,
        Bellman_optimization_rerun=False,Bellman_upper_certificate_inherited=True,
        new_dimension_bound=None)


def load(stem):
    d = json.loads((ROOT/f'results/{stem}.json').read_text())
    assert d['status'].startswith('PASS'), stem
    for p, h in d['sources'].items():
        assert sha256((ROOT/p).read_bytes()).hexdigest() == h, (stem,p)
    return d


def run(prepare=False):
    old = load(OLD)
    out = arithmetic(old)
    if prepare:
        out['status'] = 'PREPARED_CONDITIONAL_RANK15_RANK16_FINITE434'
        out['geometric_inputs_final'] = False
        return out
    assert set(FINAL_LOCAL_HASHES) == {R15,R16}, 'Wait for both local final freezes.'
    for stem, hashes in FINAL_LOCAL_HASHES.items():
        assert len(hashes) == 3
        for (folder, ext), h in zip((('scripts','py'),('results','json'),('results','md')), hashes):
            assert sha256((ROOT/f'{folder}/{stem}.{ext}').read_bytes()).hexdigest() == h
    a, b, sp = load(R15), load(R16), load(SPARSE)
    assert a['new_costs'] == {'A15':'104','B15':'104','D15':'104'}
    assert b['new_costs'] == {'A16':'113','D16':'225/2','A17':'115'}
    assert a['new_cross_component_coupons'] == b['new_cross_component_coupons'] == []
    assert a['simultaneous_rank16_exchange'] == b['new_costs']
    assert a['preserved_D16_alternative_pairs'] == [[1,14],[1,15]]
    assert a['preserved_D16_alternative_weight'] == '1/2'
    assert a['triangle_minimum'] == b['triangle_minimum'] == 2
    assert a['quadrilateral_minimum'] == b['quadrilateral_minimum'] == 1
    assert b['D16_half_pairs'] == [[2,13],[2,14],[2,15]]
    assert b['A17_restored_half_pairs'] == [[1,12],[4,15]]
    retained = {(1,14),(1,15)}
    d16 = next(r for r in sp['rules'] if r['kind']=='D' and r['rank']==16)
    weights = {tuple(r['pair']):F(r['weight']) for r in d16['weighted_pairs']}
    assert all(weights[p] == F(1,2) for p in retained)
    assert not retained & set(map(tuple,b['D16_half_pairs']))
    stems = (OLD,R15,R16,SPARSE)
    paths = [Path(__file__)] + [ROOT/f'{folder}/{s}.{ext}' for s in stems
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    out.update(status='PASS_FINITE_RANK15_RANK16_COMPLETE_RANK62_BUDGET',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        geometric_inputs_final=True,rank15_required_D16_halves_retained=sorted(retained),
        scope='Proved local rank15/rank16 exchange, including actual A17 restoration, is applied by exact fixed-vector shifts to all frozen438 alternatives. Current coupons are unchanged; conservative D16upper113 dominates actual225/2. Complete inherited Bellman certificates remain valid because every residual item is identical. No Bellman or global optimizer, geometric census, or ambient theorem is rerun or claimed.')
    return out


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--prepare',action='store_true');ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    assert not (args.prepare and args.verify)
    out=json.loads(json.dumps(run(args.prepare)))
    p=ROOT/f'results/{TAG}.json'
    if args.prepare:
        print(out['status'],'maximum',out['finite_maximum'],'shifted',out['complete_large_vectors'],flush=True)
    elif args.verify:
        assert out==json.loads(p.read_text());print(out['status'],out['finite_maximum'],flush=True)
    else:
        assert not p.exists();p.write_text(json.dumps(out,separators=(',',':'))+'\n')
        print(out['status'],out['finite_maximum'],flush=True)
