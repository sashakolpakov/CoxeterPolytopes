#!/usr/bin/env python3
"""Actual H3 count for every existing ordinary H3 cross-payment orbit.

The separate central-root beta/c family is unchanged. No new root-pair
positions are introduced; only existing coefficients are separated by count.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_h3_center_e6_orbit_exempt_20261016 as E
M=E.M;A=E.A;S=E.S;ROOT=M.ROOT
TAG='fresh_sqrt5_h3_all_actual_count_split_20261016'
BASE='fresh_sqrt5_h3_center_e6_orbit_exempt_20261016'
OLD_N=M.N;ORIGINAL_BONUS=M.bonus
INITIAL=f'data/enumerations/{TAG}_initial.json';PRIOR_INITIAL=f'data/enumerations/{BASE}_initial.json'
CHECKPOINT=Path('/private/tmp/'+TAG+'_checkpoint.json')
ROOT_TYPES={M.R.tree_key(g,[u]):(name,len(g))for name,g in M.MODELS.items()for u in g}
OLD_FAMILIES={};ORBITS={}
for(kind,key,a,b,c),j in M.P.RULE.index.items():
    if kind!='h3cross'or c!=1 or a<1:continue
    name,rank=ROOT_TYPES[key];assert name not in M.HUBS
    OLD_FAMILIES[M.SMALL+j]=(name,key,a-1,b)
    ORBITS[key]=(name,rank)
KEYS=sorted((j,c)for j,(name,key,a,b)in OLD_FAMILIES.items()
            for c in range(2,(21-len(M.MODELS[name]))//3+1)
            if len(M.MODELS[name])+a+2*b+3*c<=21)
SPLIT={key:OLD_N+i for i,key in enumerate(KEYS)};N=OLD_N+len(KEYS)
assert len(ORBITS)==6 and len(OLD_FAMILIES)==267 and len(KEYS)==392

def bonus(kind,key,ctx):
    old=ORIGINAL_BONUS(kind,key,ctx)
    if kind!='h3cross'or ctx[2]<=1 or not old:return old
    j=next(iter(old));pair=(j,ctx[2])
    if pair not in SPLIT:
        assert ROOT_TYPES[key][1]+ctx[0]+2*min(ctx[1],2)+3*ctx[2]>21
        return {}
    return {SPLIT[pair]:F(1)}

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
        expected={}
        for j,z in old.items():
            new=j
            if j in OLD_FAMILIES and ctx[2]>1:
                assert(j,ctx[2])in SPLIT;new=SPLIT[j,ctx[2]]
            expected[new]=expected.get(new,F(0))+z
        assert literal==expected,(k,meta);changed+=old!=literal
        vals,marks=M.components(ctx,meta['marked_hub']);assembled=marks[meta['marked_component'],meta['facet_root']]
        counts=M.Counter(dict(zip(M.RESERVED,ctx[:7])))
        if meta['marked_component']in M.RESERVED:counts[meta['marked_component']]-=1
        counts.update(meta['residual_components'])
        for name,count in counts.items():assembled=M.plus(assembled,vals[name],count)
        assert assembled==literal,(k,meta)
        row={j:z for j,z in literal.items()if j>=0};row[N]=F(-1)
        vertices.append(dict(row=[[j,str(z)]for j,z in sorted(row.items())],rhs=str(-literal.get(-1,F(0))),metadata=meta))
        if k%3000==0:print('ALL_COUNT_SPLIT_REBUILT',k,'changed',changed,flush=True)
    paths=[f'scripts/{TAG}.py',PRIOR_INITIAL,f'scripts/{BASE}.py']
    path.write_text(json.dumps(dict(status='LITERAL_ALL_H3_COUNT_SPLIT_HANDOFF_NOT_A_CERTIFICATE',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},vertex_cuts=vertices,steps=[],
        latest_numeric_weights=None,changed_profiles=changed,
        scope='All six ordinary H3 cross-payment root orbits are remapped by actualH3count. Every inherited forest is recounted literally and compared to its new DP component sum.'),separators=(',',':'))+'\n')

def normalization():
    records=[];hg=M.MODELS['H3']
    for key,(name,rank)in sorted(ORBITS.items(),key=lambda x:(x[1][0],repr(x[0]))):
        g=M.MODELS[name];roots=[u for u in g if M.R.tree_key(g,[u])==key];m=len(roots)
        for h in range(1,22//rank):
            for c in range(1,(21-rank*h)//3+1):
                pairs=[(('source',i,u),('H3',j,v))for i in range(h)for u in roots for j in range(c)for v in hg]
                assert len(pairs)==m*h*3*c
                for side,fg in [('source',g),('H3',hg)]:
                    for root in fg:
                        retained=sum((side,0,root)not in pair for pair in pairs)
                        loss=F(int(root in roots))if side=='source'else F(m*h,3*c)
                        assert F(retained,3*c)==m*h-loss
                        records.append(dict(source_component=name,source_root_key=key,source_count=h,H3_count=c,
                            facet_component=name if side=='source'else'H3',facet_root=root,
                            full_pairs=len(pairs),retained_pairs=retained,beta_coefficient=str(m*h-loss)))
    return records

def barrier_releases():
    q=json.loads((ROOT/'results/fresh_sqrt5_adjacent_discount_frontier_20261016.json').read_text());records=[]
    for index,actuala,source,copies,coefficient in[(7661,2,'A4',4,8),(7668,3,'A3',5,10)]:
        row=next(r['row']for r in q['immediate_sealed_weight_obligations']if r['class_index']==index)
        es=[(a,b,m)for(a,b),m in zip(combinations(range(row['nodes']),2),row['labels'])];white=row['white']
        corners=[M.compile_corner(white+[a,b],es,[a,b],[4])for a,b in product((0,2),(1,3))]
        p,tables=M.add_bonus(corners,(1,0,4,0,0,0,0))
        for table in tables:p=M.plus(p,table[4])
        key=M.R.tree_key(M.MODELS[source],[next(u for u in M.MODELS[source]if len(M.MODELS[source][u])==1)])
        old=next(iter(ORIGINAL_BONUS('h3cross',key,(actuala,0,1,0,0,0,0,0))));new=SPLIT[old,5]
        assert p=={-1:F(1),1:F(-1),new:F(1,15)}
        ctx=[actuala,0,1,0,copies if source=='A4'else 0,0,0,0]
        front=dict(context=ctx,marked_hub='A1',marked_component='A1',facet_root=0,residual_components=[source]*copies if source=='A3'else[])
        vertex=M.literal_profile(front);assert vertex[old]==coefficient and new not in vertex
        records.append(dict(face_index=index,count1_coefficient=old,count5_coefficient=new,
            face_profile=[[j,str(z)]for j,z in sorted(p.items())],count1_vertex_beta_coefficient=str(coefficient)))
    return records

def run(saved=None):
    install();A.TAG=TAG;A.INITIAL=INITIAL;A.CHECKPOINT=CHECKPOINT
    S.TAG=TAG;S.CHECKPOINT=CHECKPOINT;S.RESUME=True;S.checkpoint=A.checkpoint
    audits=normalization();releases=barrier_releases()
    obligations,counts=A.faces();print('MODEL',N,'variables',len(obligations),'profiles',counts,flush=True)
    cert,steps=S.solve(obligations)if saved is None else(saved['certificate'],saved['steps'])
    proof=M.verify(obligations,cert)
    result=json.loads((ROOT/f'results/{BASE}.json').read_text());order=result['variable_order'];assert len(order)==OLD_N
    for j in OLD_FAMILIES:order[j]=dict(order[j],actual_H3_count=1)
    for j,c in KEYS:
        name,key,a,b=OLD_FAMILIES[j]
        order.append(dict(kind='actual_H3_count_split',source_component=name,source_root_key=key,
            actual_A1_count=a,I25_context=b,actual_H3_count=c,inherited_count1_index=j,
            inherited_orbit=order[j]['orbit'],pair_normalization=str(F(1,3*c))))
    assert len(order)==N
    paths=[f'scripts/{TAG}.py',INITIAL,PRIOR_INITIAL,f'scripts/{BASE}.py',f'results/{BASE}.json',
        'results/fresh_sqrt5_h3_multiplicity_barrier_20261016.json','results/fresh_sqrt5_h3_a3_multiplicity_barrier_20261016.json']
    result['sources'].update({p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths})
    result.update(status='PASS_EXACT_ALL_H3_CROSS_ACTUAL_COUNT_PROPOSAL',variables=N,variable_order=order,
        face_counts=counts,unique_face_profiles=len(obligations),certificate=cert,verification=proof,steps=steps,
        added_actual_count_coefficients=len(KEYS),split_normalization=audits,multiplicity_barrier_releases=releases,
        split_source_orbits=[dict(component=name,root_key=key,rank=rank)for key,(name,rank)in sorted(ORBITS.items(),key=lambda x:(x[1][0],repr(x[0])))],
        scope='Every existing ordinary H3 cross-payment orbit is split by actual H3 count: six A3/A4/D4/D5 source-root orbits. Count1 uses old indices; counts2..6 use392new independent nonnegative coefficients. Per-pair normalization stays1/(3c), and the separate H3-central-root beta/c family is unchanged. No new root-pair positions. E6two-edge exemption, distance2discount0,adjacent<=2/3,996triangle states,37paidquadstops,oldzeroquads and fullmarkedDP retained. Remaining paidquadextension coverage is open; no dimension claim.')
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();install()
    if not args.verify:
        rebuilt_initial()
        if not CHECKPOINT.exists():CHECKPOINT.write_bytes((ROOT/INITIAL).read_bytes())
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text())if args.verify else None
    result=json.loads(json.dumps(run(saved),default=str))
    if args.verify:assert result==saved
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    proof=result['verification'];number=proof.get('maximum',proof.get('lower_bound'))
    print(result['status'],result['certificate']['kind'],float(F(number)),flush=True)
