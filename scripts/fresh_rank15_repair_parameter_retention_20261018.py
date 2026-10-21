#!/usr/bin/env python3
"""Move the two quarter repairs along their nonnegative half-sum segment.

The four required branches use half+x+y. Old routes already had both
rank15 orbits zero. No forest optimum or new dimension bound is computed.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
import fresh_a1_d5_witness_context_retention_20261018 as W

ROOT=W.ROOT;Q=W.C.Q;C=W.C;R=W.R
TAG='fresh_rank15_repair_parameter_retention_20261018'
OLD=F(1,4);TARGET=F(1,2);ENDPOINTS=(F(),F(1,2))
STEMS=(W.TAG,Q.TAG,R.TAG,Q.N.TAG,Q.N.S,Q.N.A,Q.N.D)

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def actual_base(t,x=TARGET):
    assert 0<=x<=F(1,2)
    return W.actual_base(t)+(128*(x-OLD)if t==('D',15)else-128*(x-OLD)if t==('A',15)else 0)
actual_cross=W.actual_cross
context=W.context
def value(fs,x=TARGET):
    return sum(actual_base(t,x)for t in fs)+sum(actual_cross(a,b)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs)

def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    q=docs[Q.TAG];assert q['restored_weight']=='1/4'
    assert q['A15_restored_pairs']==[[0,13],[1,14]]and q['D15_restored_pairs']==[[1,13],[1,14]]
    variants=[];required=[]
    for p in q['whole_variants']:
        corners=[]
        for c in p['corners']:
            assert c['quarter_role']==Q.N.role(c)
            row=dict(marks=c['marks'],role=c['quarter_role'],current_lower=c['new'])
            weights={}
            for x in(*ENDPOINTS,OLD):
                y=F(1,2)-x
                w=x if c['quarter_role']=='D15'else y if c['quarter_role']=='A15'else F(c['new'])
                assert w>=0;weights[str(x)]=str(w)
            row['parameter_lower']=weights;corners.append(row)
        totals={str(x):sum(F(c['parameter_lower'][str(x)])for c in corners)for x in(*ENDPOINTS,OLD)}
        assert all(v>=1 for v in totals.values())
        if any(c['quarter_role']=='A16'for c in p['corners']):
            assert all(v==1 for v in totals.values());required.append([p['index'],p['extra_contact']])
        variants.append(dict(index=p['index'],extra_contact=p['extra_contact'],corners=corners,totals=totals))
    assert len(variants)==24 and required==[[673,None],[1558,None],[1646,13],[1933,3]]
    tables={}
    for typ,pairs in(('A15',q['A15_restored_pairs']),('D15',q['D15_restored_pairs'])):
        old={tuple(r['pair']):F(r['weight'])for r in q['weighted_tables'][typ]}
        assert all(old[tuple(p)]==OLD for p in pairs)
        target=dict(old)
        for p in pairs:target[tuple(p)]=TARGET if typ=='D15'else F(1,2)-TARGET
        expected=F(101)if typ=='D15'else F(100)
        assert sum(target.values())==expected==F(actual_base((typ[0],15)),64)
        tables[typ]=[dict(pair=p,weight=w)for p,w in sorted(target.items())]
    for x in(*ENDPOINTS,OLD):
        for t in C.TYPES:assert actual_base(t,x)>=0
        assert F(actual_base(('A',15),x),64)==101-2*x
        assert F(actual_base(('D',15),x),64)==100+2*x
    # The four rank13 seals never spend a restored A15/D15 long orbit.
    for seal in docs[R.TAG]['local_seals']:
        for c in seal['corners']:
            if F(c['parameter_values']['9/64'])==0:continue
            blocks=c['components'];ms=c['marks'];aa,bb=[next(b for b in blocks if m in b['nodes'])for m in ms]
            if aa is bb and aa['type']=='A15':
                assert not set(R.T.mark_pairs(R.A.matrix(seal['gram']),aa['nodes'],ms,15,'A'))&{(0,13),(1,14)}
    witness=[('A',16)]*2+[('D',26)]
    assert value(witness)==W.value(witness)==64*F(761,2)
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md'))}
    return W.enc(dict(status='PASS_RANK15_REPAIR_PARAMETER_FACE_RETENTION',sources={p:digest(p)for p in sorted(paths)},
        baseline=W.TAG,scale=64,parameter_interval=list(ENDPOINTS),chosen_D15_weight=TARGET,
        chosen_A15_weight='0',constraint='x+y=1/2; x,y>=0',whole_variants=variants,
        required_half_sum_variants=required,weighted_tables=tables,
        actual_costs={'A15':'100','D15':'101','B15':'97','A16':'441/4','D16':'439/4'},
        family_actual_costs={'A15':'101-2*x','D15':'100+2*x'},representative_rank15_cost='max(101-2*x,100+2*x)',
        chosen_rank15_representative_cost='101',finite_actual_types=len(C.TYPES),
        ideal_weights_changed=False,cross_or_context_weights_changed=False,
        triangle_minimum='2',quadrilateral_minimum='1',
        completed_dimension_bound=58,new_dimension_bound=None,forest_count_performed=False,
        global_maximum_claimed=False,
        scope='Affine endpoint retention of all24 sealed D14 variants and complete original-zero chosen-output argument, with current rank13 and guarded A1/D5 compatibility. Chosen endpoint restores only D15 halves and returns A15 orbits to their earlier zeros. Complete fresh counting is still required.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['actual_costs'],flush=True)
