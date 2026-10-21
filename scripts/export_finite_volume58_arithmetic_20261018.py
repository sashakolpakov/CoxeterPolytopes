#!/usr/bin/env python3
"""Export frozen dimension58 arithmetic as plain integer/rational data.

No optimization, local geometry or convolution is run. The export is not a
proof: a proof assistant must verify the stated domains and recurrences.
"""
import argparse,json
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import fresh_a15_d15_ambient59_forest_closure_20261018 as F

ROOT=Path(__file__).resolve().parents[1]
TAG='finite_volume58_arithmetic_20261018'
BINDER='fresh_a15_d15_dimension58_20261018'
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def rat(x):
    q=Fraction(x);return [q.numerator,q.denominator]
def run():
    p=f'results/{F.TAG}.json';saved=json.loads((ROOT/p).read_text())
    for path,h in saved['sources'].items():assert digest(path)==h,path
    binder=json.loads((ROOT/f'results/{BINDER}.json').read_text())
    assert binder['dimension_bound']==58 and binder['sources'][p]==digest(p)
    rep=F.C.C.P.Z.B.B.representative
    types=[t for t in F.C.C.P.Z.B.B.types()if t[1]<=59]
    reps=sorted(set(map(rep,types)))
    ti={t:i for i,t in enumerate(types)};ri={t:i for i,t in enumerate(reps)}
    assert len(types)==180
    tabs=[]
    for t in saved['finite']['tables']:
        key=tuple(t['key']);k5=t['fixed_A5_count'];rs=[r for r in F.ITEMS if not(t['exclude2']and r==2)]
        tabs.append(dict(key=t['key'],fixed_rank5_count=k5,exclude_rank2=t['exclude2'],
            item_ranks=rs,item_weights=[F.item(r,key,k5)for r in rs],values=t['values'],last_ranks=t['last']))
    assert len(tabs)==5148
    rows=[]
    for r in saved['finite']['expanded_alternatives']:
        cs=r['counts'];fixed=[('B'if n==9 else'A'if n==17 else'D',n)for n,k in zip(range(9,18),cs)for _ in range(k)]+[('A',5)]*r['fixed_A5_count']
        fixedcost=sum(F.upper_base(t)for t in fixed)+sum(F.upper_cross(a,b)for i,a in enumerate(fixed)for b in fixed[i+1:])
        rows.append(dict(counts=cs,key=r['key'],fixed_rank5_count=r['fixed_A5_count'],
            fixed_representative_indices=[ri[t]for t in fixed],fixed_cost=fixedcost,
            residual_rank=59-sum(t[1]for t in fixed),context_options=r['context_options'],
            winning_option=r['winning_option'],upper=r['sixteenths_upper'],actual_literal=r['sixteenths_actual']))
    ideal=saved['ideal'];profiles=[]
    for r in ideal['rows']:
        profiles.append(dict(ranks=r['ranks'],coefficients=r['h_polynomial'],V=r['link_vertices'],G=r['G'],
            curvature=rat(r['curvature']),normalized_cost=rat(r['normalized_cost']),
            normalized_margin=rat(r['normalized_margin']),absolute_margin=rat(r['absolute_margin'])))
    source_paths=[f'scripts/{F.TAG}.py',p,f'results/{F.TAG}.md',f'scripts/{BINDER}.py',f'results/{BINDER}.json',f'results/{BINDER}.md',
        'scripts/export_finite_volume58_arithmetic_20261018.py']
    return dict(schema='finite-volume58-arithmetic-v1',status='EXPORTED_FROZEN_DATA_NOT_A_FORMAL_PROOF',
        sources={p:digest(p)for p in source_paths},ambient_dimension=59,link_rank=58,finite_scale=16,
        actual_types=types,representative_types=reps,representative_map=[ri[rep(t)]for t in types],
        actual_bases=[F.actual_base(t)for t in types],upper_bases=[F.upper_base(t)for t in reps],
        actual_cross=[[F.actual_cross(a,b)for b in types]for a in types],
        upper_cross=[[F.upper_cross(a,b)for b in reps]for a in reps],
        residual_representatives=[[r,ri[rep(('A'if r<=8 else'D',r))]]for r in F.ITEMS],
        fixed_large_ranks=list(range(9,18)),fixed_rank17_type='A17',residual_rank17_type='D17',
        old_context=dict(guard_rank10=True,minimum_A1=2,scaled_coefficient=4,formula='j*(j-1)'),
        new_context=dict(guard_rank9=True,minimum_A1=1,scaled_coefficient=16,formula='j*k5'),
        fixed_vector_count=961,expanded_vector_count=2574,tables=tabs,expanded_rows=rows,
        finite_upper_scaled=6264,finite_upper=rat('783/2'),actual_attainment_claimed=False,
        normalized_ideal_envelopes=[[int(r),rat(c)]for r,c in ideal['current_envelopes'].items()],
        terminal_minimum_large_rank=9,terminal_profile_count=1890,ideal_profiles=profiles,
        normalized_margin_minimum=rat('1/340'),absolute_margin_minimum=rat('91/2'),
        finite_curvature=rat('783/2'),extremal_slack_sum=870,global_correction=4350,
        exclusions=['universal local geometry','supporting-root census completeness','all-rank concentration theorem',
            'cusp-factor geometric interpretation','shape-curvature geometric interpretation','higher-dimension tail'],
        intended_verification='Recompute complete domains, type domination, item identities, contextual coverage, Bellman upper inequalities, polynomial products and exact margins in the kernel; prove the excluded internal lemmas separately.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'data/certificates/{TAG}.json';out=json.loads(json.dumps(run()))
    if args.verify:assert out==json.loads(p.read_text())
    else:
        assert not p.exists();p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],len(out['tables']),len(out['ideal_profiles']),flush=True)
