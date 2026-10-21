#!/usr/bin/env python3
"""Exempt exactly the E6 trivalent-to-degree-two edge orbit from discount."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_h3_center_complete_triangle_repair_20261016 as A

M=A.M;S=A.S;ROOT=M.ROOT
TAG='fresh_sqrt5_h3_center_e6_orbit_exempt_20261016'
INITIAL=f'data/enumerations/{TAG}_initial.json'
CHECKPOINT=Path('/private/tmp/'+TAG+'_checkpoint.json')
PRIOR_INITIAL='data/enumerations/fresh_sqrt5_h3_center_complete_triangle_repair_20261016_initial.json'
ORIGINAL_ELIGIBLE=M.P.D.eligible
e6=M.MODELS['E6'];branch=next(u for u in e6 if len(e6[u])==3)
neighbor=next(u for u in e6[branch]if len(e6[u])==2)
EXEMPT=M.R.tree_key(e6,[branch,neighbor])
EXEMPT_PAIRS=[(a,b)for a,b in combinations(e6,2)if M.R.tree_key(e6,[a,b])==EXEMPT]
assert len(EXEMPT_PAIRS)==2

def eligible(policy,name,g,a,b):
    if policy=='long_adjacent3'and name=='E6'and M.R.tree_key(g,[a,b])==EXEMPT:return False
    return ORIGINAL_ELIGIBLE(policy,name,g,a,b)

def install():
    M.P.D.eligible=eligible
    # Imports may have formed compensation scouts; never reuse a cached
    # static/marked payment evaluated under the old eligibility predicate.
    for module in(M.P,M.B,M.C,M.D,M.E,M):
        for name in('static_component','static','components','numerical','floating_components'):
            f=getattr(module,name,None)
            if hasattr(f,'cache_clear'):f.cache_clear()
    changes=[]
    for name,g in M.MODELS.items():
        for a,b in combinations(g,2):
            old=ORIGINAL_ELIGIBLE('long_adjacent3',name,g,a,b);new=eligible('long_adjacent3',name,g,a,b)
            if old!=new:changes.append(dict(component=name,pair=[a,b],old=old,new=new))
    assert len(changes)==2 and all(r['component']=='E6'and r['old']and not r['new']for r in changes)
    return changes

def prepare_initial():
    path=ROOT/INITIAL
    if path.exists():return
    prior=json.loads((ROOT/PRIOR_INITIAL).read_text());vertices=[];changed=0
    for k,rec in enumerate(prior['vertex_cuts']):
        meta=rec['metadata'];profile=M.literal_profile(meta)
        old={int(j):F(z)for j,z in rec['row']if int(j)!=M.N}
        if F(rec['rhs']):old[-1]=-F(rec['rhs'])
        assert all(profile.get(j,0)==old.get(j,0)for j in set(profile)|set(old)if j!=1)
        difference=profile.get(1,0)-old.get(1,0);assert difference>=0 and difference.denominator==1
        changed+=bool(difference)
        row={j:z for j,z in profile.items()if j>=0};row[M.N]=F(-1)
        vertices.append(dict(row=[[j,str(z)]for j,z in sorted(row.items())],rhs=str(-profile.get(-1,F(0))),metadata=meta))
        if k%3000==0:print('REBUILT_VERTEX_CUTS',k,'changed',changed,flush=True)
    paths=[f'scripts/{TAG}.py',PRIOR_INITIAL,'scripts/fresh_sqrt5_marked_h3_center_cross_20261016.py']
    doc=dict(status='LITERAL_REBUILT_VERTEX_CUT_HANDOFF_NOT_A_CERTIFICATE',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},vertex_cuts=vertices,
        steps=[],latest_numeric_weights=None,changed_discount_profiles=changed,
        scope='All inherited forest metadata are rebuilt by literal pair enumeration under the new eligibility predicate. Only adjacent-discount coefficients change; old numerical weights/solve steps are not inherited.')
    path.write_text(json.dumps(doc,separators=(',',':'))+'\n')

def run(saved=None):
    changes=install()
    A.TAG=TAG;A.INITIAL=INITIAL;A.CHECKPOINT=CHECKPOINT
    result=A.run(saved)
    paths=['scripts/fresh_sqrt5_h3_center_complete_triangle_repair_20261016.py',
           'scripts/fresh_sqrt5_uniform_adjacent_discount_barrier_20261016.py',
           'results/fresh_sqrt5_uniform_adjacent_discount_barrier_20261016.json',PRIOR_INITIAL]
    result['sources'].update({p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths})
    result['status']='PASS_EXACT_E6_ORBIT_EXEMPT_PARTIAL_PAID_QUAD_PROPOSAL'
    result['eligibility_changes']=changes;result['exempt_pair_orbit']=EXEMPT
    result['variable_order'][1]=dict(kind='adjacent3_discount',minimum_component_rank=6,
        excluded_component='E6',excluded_pair_orbit=EXEMPT,excluded_pair_description='The two edges from the trivalent root to degree-two neighbors.')
    result['scope']='Same nonnegative bonus variables, distance-two discount0 and uniform adjacent discount at most2/3, except exactly one two-edge E6 orbit is exempt. All literal component totals, facet-root losses,996triangle states,37paid-quad stopping rows and inherited baseline-zero quad profiles are recomputed. Full marked finite-forest budgets retained; inherited cuts rebuilt by literal root pairs. Further paid-quadrilateral extension coverage remains open. No dimension-bound promotion.'
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    install()
    if not args.verify:
        prepare_initial()
        if not CHECKPOINT.exists():CHECKPOINT.write_bytes((ROOT/INITIAL).read_bytes())
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text())if args.verify else None
    result=json.loads(json.dumps(run(saved),default=str))
    if args.verify:assert result==saved
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    proof=result['verification'];number=proof.get('maximum',proof.get('lower_bound'))
    print(result['status'],result['certificate']['kind'],float(F(number)),flush=True)
