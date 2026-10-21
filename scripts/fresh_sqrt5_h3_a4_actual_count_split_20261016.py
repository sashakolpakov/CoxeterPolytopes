#!/usr/bin/env python3
"""Targeted actual-H3-count split of the H3 x A4-end payment only."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_h3_center_e6_orbit_exempt_20261016 as E

M=E.M;A=E.A;S=E.S;ROOT=M.ROOT
TAG='fresh_sqrt5_h3_a4_actual_count_split_20261016'
OLD_N=M.N;ORIGINAL_BONUS=M.bonus
BASE='fresh_sqrt5_h3_center_e6_orbit_exempt_20261016'
INITIAL=f'data/enumerations/{TAG}_initial.json'
PRIOR_INITIAL=f'data/enumerations/{BASE}_initial.json'
CHECKPOINT=Path('/private/tmp/'+TAG+'_checkpoint.json')
KEYS=sorted((a,b,c)for a in range(22)for b in range(3)for c in range(2,6)
            if a+2*b+3*c+4<=21 and ORIGINAL_BONUS('h3cross',M.A4_ENDPOINT,(a,b,c,0,1,0,0,0)))
SPLIT={key:OLD_N+i for i,key in enumerate(KEYS)};N=OLD_N+len(KEYS)

def bonus(kind,key,ctx):
    old=ORIGINAL_BONUS(kind,key,ctx)
    if kind!='h3cross'or key!=M.A4_ENDPOINT or ctx[2]<=1 or not old:return old
    k=(ctx[0],min(ctx[1],2),ctx[2])
    if k not in SPLIT:
        # Component-cost tables also form unused components that cannot fit
        # beside these reserved counts; no actual A4 corner can enter here.
        assert k[0]+2*k[1]+3*k[2]+4>21
        return {}
    return {SPLIT[k]:F(1)}

def install():
    E.install()
    for module in(M.B,M.C,M.D,M.E,M):module.bonus=bonus
    M.N=N

def rebuilt_initial():
    path=ROOT/INITIAL
    if path.exists():return
    prior=json.loads((ROOT/PRIOR_INITIAL).read_text());vertices=[];changed=0
    for k,rec in enumerate(prior['vertex_cuts']):
        meta=rec['metadata'];ctx=tuple(meta['context']);literal=M.literal_profile(meta)
        old={int(j):F(z)for j,z in rec['row']if int(j)!=OLD_N}
        if F(rec['rhs']):old[-1]=-F(rec['rhs'])
        original=ORIGINAL_BONUS('h3cross',M.A4_ENDPOINT,ctx);replacement=bonus('h3cross',M.A4_ENDPOINT,ctx)
        expected=old.copy()
        if original and replacement and original!=replacement:
            before=next(iter(original));after=next(iter(replacement));amount=expected.pop(before,F(0))
            if amount:expected[after]=expected.get(after,F(0))+amount;changed+=1
        assert literal==expected,(k,meta)
        # Reassemble the DP component sum independently of literal pair loops.
        vals,marks=M.components(ctx,meta['marked_hub']);assembled=marks[meta['marked_component'],meta['facet_root']]
        counts=M.Counter(dict(zip(M.RESERVED,ctx[:7])))
        if meta['marked_component']in M.RESERVED:counts[meta['marked_component']]-=1
        counts.update(meta['residual_components'])
        for name,count in counts.items():assembled=M.plus(assembled,vals[name],count)
        assert assembled==literal,(k,meta)
        row={j:z for j,z in literal.items()if j>=0};row[N]=F(-1)
        vertices.append(dict(row=[[j,str(z)]for j,z in sorted(row.items())],rhs=str(-literal.get(-1,F(0))),metadata=meta))
        if k%3000==0:print('COUNT_SPLIT_REBUILT',k,'changed',changed,flush=True)
    paths=[f'scripts/{TAG}.py',PRIOR_INITIAL,f'scripts/{BASE}.py']
    path.write_text(json.dumps(dict(status='LITERAL_COUNT_SPLIT_HANDOFF_NOT_A_CERTIFICATE',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},vertex_cuts=vertices,steps=[],
        latest_numeric_weights=None,changed_profiles=changed,
        scope='Only the targeted H3xA4-end coefficient is remapped by actualH3count. Every inherited forest is recounted literally and compared to its new DP component sum.'),separators=(',',':'))+'\n')

def normalization():
    records=[];ag=M.MODELS['A4'];hg=M.MODELS['H3'];ends=[u for u in ag if len(ag[u])==1]
    for h in range(1,5):
        for c in range(1,6):
            if 4*h+3*c>21:continue
            pairs=[(('A4',i,u),('H3',j,v))for i in range(h)for u in ends for j in range(c)for v in hg]
            assert len(pairs)==6*h*c
            for name,g in [('A4',ag),('H3',hg)]:
                for root in g:
                    retained=sum((name,0,root)not in pair for pair in pairs)
                    expected=F(2*h)- (F(1)if name=='A4'and root in ends else F(2*h,3*c)if name=='H3'else F(0))
                    assert F(retained,3*c)==expected
                    records.append(dict(A4_count=h,H3_count=c,facet_component=name,facet_root=root,
                        full_pairs=len(pairs),retained_pairs=retained,beta_coefficient=str(expected)))
    return records

def barrier_release():
    quad=json.loads((ROOT/'results/fresh_sqrt5_adjacent_discount_frontier_20261016.json').read_text())
    row=next(r['row']for r in quad['immediate_sealed_weight_obligations']if r['class_index']==7661)
    es=[(a,b,m)for(a,b),m in zip(combinations(range(row['nodes']),2),row['labels'])];white=row['white']
    corners=[M.compile_corner(white+[a,b],es,[a,b],[4])for a,b in product((0,2),(1,3))]
    p,tables=M.add_bonus(corners,(1,0,4,0,0,0,0))
    for table in tables:p=M.plus(p,table[4])
    old=next(iter(ORIGINAL_BONUS('h3cross',M.A4_ENDPOINT,(2,0,1,0,4,0,0,0))))
    new=SPLIT[2,0,5]
    assert p=={-1:F(1),1:F(-1),new:F(1,15)}and old!=new
    front=dict(context=[2,0,1,0,4,0,0,0],marked_hub='A1',marked_component='A1',facet_root=0,residual_components=[])
    v=M.literal_profile(front);assert v[old]==8 and new not in v
    return dict(face_old_coefficient=old,face_new_coefficient=new,face_profile=[[j,str(z)]for j,z in sorted(p.items())],
                count1_vertex_beta_coefficient='8',count5_vertex_beta_coefficient='0')

def run(saved=None):
    install();A.TAG=TAG;A.INITIAL=INITIAL;A.CHECKPOINT=CHECKPOINT
    S.TAG=TAG;S.CHECKPOINT=CHECKPOINT;S.RESUME=True;S.checkpoint=A.checkpoint
    audits=normalization();release=barrier_release()
    obligations,counts=A.faces();print('MODEL',N,'variables',len(obligations),'profiles',counts,flush=True)
    cert,steps=S.solve(obligations)if saved is None else(saved['certificate'],saved['steps'])
    proof=M.verify(obligations,cert)
    result=json.loads((ROOT/f'results/{BASE}.json').read_text())
    order=result['variable_order'];assert len(order)==OLD_N
    target_old={next(iter(ORIGINAL_BONUS('h3cross',M.A4_ENDPOINT,(a,b,1,0,1,0,0,0))))for a,b,c in KEYS}
    for j in target_old:order[j]=dict(order[j],actual_H3_count=1)
    for a,b,c in KEYS:
        base=next(iter(ORIGINAL_BONUS('h3cross',M.A4_ENDPOINT,(a,b,c,0,1,0,0,0))))
        order.append(dict(kind='actual_H3_count_split',source_component='A4',source_root_orbit='endpoint',
            source_root_key=M.A4_ENDPOINT,actual_A1_count=a,I25_context=b,actual_H3_count=c,
            inherited_count1_index=base,pair_normalization=f'1 / {3*c}'))
    assert len(order)==N
    paths=[f'scripts/{TAG}.py',INITIAL,PRIOR_INITIAL,f'scripts/{BASE}.py',f'results/{BASE}.json',
        'results/fresh_sqrt5_h3_multiplicity_barrier_20261016.json']
    result['sources'].update({p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths})
    result.update(status='PASS_EXACT_TARGETED_H3_A4_COUNT_SPLIT_PROPOSAL',variables=N,variable_order=order,
        face_counts=counts,unique_face_profiles=len(obligations),certificate=cert,verification=proof,steps=steps,
        added_actual_count_coefficients=len(KEYS),split_normalization=audits,multiplicity_barrier_release=release,
        scope='Exactly the H3xA4-end inherited family is split by actual H3 count; all3H3roots retain equal1/(3c) per-pair payments. Count1 uses old indices; counts2..5 use new independent nonnegative coefficients. All other bonuses unchanged, E6 two-edge exemption retained, distance2 discount0 and adjacent<=2/3. All996triangle states,37paid-quad stopping rows,oldzeroquads and complete marked-forest DP retained. Further paid-quad extension coverage remains open; no dimension claim.')
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    install()
    if not args.verify:
        rebuilt_initial()
        if not CHECKPOINT.exists():CHECKPOINT.write_bytes((ROOT/INITIAL).read_bytes())
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text())if args.verify else None
    result=json.loads(json.dumps(run(saved),default=str))
    if args.verify:assert result==saved
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    proof=result['verification'];number=proof.get('maximum',proof.get('lower_bound'))
    print(result['status'],result['certificate']['kind'],float(F(number)),flush=True)
