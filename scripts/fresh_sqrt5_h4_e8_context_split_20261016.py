#!/usr/bin/env python3
"""H4 compensation indexed by actual E8 presence in the full vertex forest."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_h3_all_actual_count_split_20261016 as E
M=E.M;A=E.A;S=E.S;ROOT=M.ROOT
TAG='fresh_sqrt5_h4_e8_context_split_20261016';BASE=E.TAG
PRIOR=json.loads((ROOT/f'results/{BASE}.json').read_text())
OLD_N=E.N
KEYS=[j for j,v in enumerate(PRIOR['variable_order'])if isinstance(v,dict)and v.get('type')=='H4'
      and(v.get('H4_count')in(1,2)or(v.get('H4_count')==3 and v.get('A1_presence')==1))]
SPLIT={j:OLD_N+i for i,j in enumerate(KEYS)};N=OLD_N+len(KEYS)
assert len(KEYS)==108
OLD_COMPONENTS=M.components;OLD_LITERAL=M.literal_profile
OLD_COMPILE=M.compile_corner;OLD_ADD=M.add_bonus
INITIAL=f'data/enumerations/{TAG}_initial.json'
CHECKPOINT=Path('/private/tmp/'+TAG+'_checkpoint.json')

def remap(profile,present):
    return {SPLIT.get(j,j):z for j,z in profile.items()}if present else profile.copy()

def actual_flag(front):
    return int(front['context'][3]in(1,2,3)and
               (front['marked_component']=='E8'or'E8'in front['residual_components']))

@lru_cache(None)
def components(ctx,mh):
    assert len(ctx)==9
    vals,marks=OLD_COMPONENTS(ctx[:8],mh)
    return {n:remap(p,ctx[8])for n,p in vals.items()},{n:remap(p,ctx[8])for n,p in marks.items()}

def literal_profile(front):
    flag=actual_flag(front)
    if len(front['context'])==9:assert front['context'][8]==flag
    p=OLD_LITERAL(dict(front,context=front['context'][:8]))
    return remap(p,flag)

def compile_corner(vertices,edges,pair,marks):
    old=OLD_COMPILE(vertices,edges,pair,marks)
    e8=sum(M.R.type_of(g)=='E8'for g in M.R.components(M.R.graph(vertices,edges)))
    return old+(e8,)

def add_bonus(compiled,ext):
    assert len(ext)==8
    out={};tables=[]
    for corner in compiled:
        common,one=OLD_ADD([corner[:9]],ext[:7])
        assert len(one)==1
        present=bool(corner[9]or ext[7])
        out=M.plus(out,remap(common,present))
        tables.append({mark:remap(p,present)for mark,p in one[0].items()})
    return out,tables

@lru_cache(None)
def external_contexts(rank):
    # E8 is counted exactly outside the connected face core. All remaining
    # components exclude the seven old reserved types and E8. Every residual
    # rank0 or >=2 is attainable using A2/A3; residual rank1 is impossible.
    out=[]
    for e8 in range(rank//8+1):
        for ctx in M.contexts(rank-8*e8):
            residual=rank-8*e8-sum(r*c for r,c in zip(M.RANKS,ctx))
            assert residual==0 or residual>=2
            out.append(ctx+(e8,))
    return tuple(out)

def oracle(x,exact=False):
    fronts=[]
    for base in M.contexts(21):
        oldctx=base+(int(base[3]==4 and(base[5]or base[6])),)
        remaining=21-sum(r*c for r,c in zip(M.RANKS,base))
        tracked=base[3]in(1,2,3)
        for flag in(range(2)if tracked and remaining>=8 else(0,)):
            ctx=oldctx+(flag,)
            modes=[None]+[mode for hub,count in zip(M.RESERVED,base)if count
                         for mode in(('H3_center','H3_leaf')if hub=='H3'else(hub,))]
            for mh in modes:
                vals,marks=components(ctx,mh)
                if exact:
                    costs={n:M.value(p,x)for n,p in vals.items()}
                    mc={n:M.value(p,x)for n,p in marks.items()}
                else:
                    nv,nm=M.numerical(ctx,mh)
                    evaluate=lambda p:p[0]+sum(z*x[j]for j,z in p[1])
                    costs={n:evaluate(p)for n,p in nv.items()};mc={n:evaluate(p)for n,p in nm.items()}
                dp={(0,0):(F(0)if exact else 0.0,())}
                for rank in range(1,remaining+1):
                    for name,g in M.MODELS.items():
                        if name in M.RESERVED:continue
                        for oldflag in(0,1):
                            entry=dp.get((rank-len(g),oldflag))
                            if entry is None:continue
                            newflag=int(bool(oldflag or(tracked and name=='E8')))
                            score=entry[0]+costs[name];state=(rank,newflag)
                            if state not in dp or score>dp[state][0]:dp[state]=(score,entry[1]+(name,))
                const=sum(c*costs[name]for name,c in zip(M.RESERVED,base));best=None;winner=None
                for(name,e),mv in marks.items():
                    if name in M.RESERVED:
                        if not base[M.RESERVED.index(name)]:continue
                        rank=remaining;score=const-costs[name]+mc[name,e]
                    else:rank=remaining-len(M.MODELS[name]);score=const+mc[name,e]
                    for residualflag in(0,1):
                        if int(bool(residualflag or(tracked and name=='E8')))!=flag:continue
                        if(rank,residualflag)not in dp:continue
                        total=score+dp[rank,residualflag][0]
                        if best is None or total>best:
                            best=total;winner=dict(context=ctx,marked_hub=mh,marked_component=name,
                                facet_root=e,residual_components=dp[rank,residualflag][1])
                if best is None:continue
                p=marks[winner['marked_component'],winner['facet_root']]
                counts=Counter(dict(zip(M.RESERVED,base)))
                if winner['marked_component']in M.RESERVED:counts[winner['marked_component']]-=1
                counts.update(winner['residual_components'])
                for name,count in counts.items():p=M.plus(p,vals[name],count)
                assert actual_flag(winner)==flag
                fronts.append(dict(**winner,maximum=best,profile=p))
    return fronts

def install():
    E.install()
    M.N=N;M.components=components;M.literal_profile=literal_profile;M.oracle=oracle
    M.compile_corner=compile_corner;M.add_bonus=add_bonus
    components.cache_clear();M.numerical.cache_clear()

def rebuild_initial():
    path=ROOT/INITIAL
    if path.exists():return
    prior=json.loads((ROOT/E.INITIAL).read_text());vertices=[];changed=0
    for k,rec in enumerate(prior['vertex_cuts']):
        meta=rec['metadata'];meta=dict(meta,context=list(meta['context'][:8])+[actual_flag(meta)])
        p=literal_profile(meta)
        old={int(j):F(z)for j,z in rec['row']if int(j)!=OLD_N}
        if F(rec['rhs']):old[-1]=-F(rec['rhs'])
        assert p==remap(old,meta['context'][8]);changed+=p!=old
        vals,marks=components(tuple(meta['context']),meta['marked_hub'])
        assembled=marks[meta['marked_component'],meta['facet_root']]
        counts=Counter(dict(zip(M.RESERVED,meta['context'][:7])))
        if meta['marked_component']in M.RESERVED:counts[meta['marked_component']]-=1
        counts.update(meta['residual_components'])
        for name,count in counts.items():assembled=M.plus(assembled,vals[name],count)
        assert assembled==p
        row={j:z for j,z in p.items()if j>=0};row[N]=F(-1)
        vertices.append(dict(row=[[j,str(z)]for j,z in sorted(row.items())],rhs=str(-p.get(-1,F(0))),metadata=meta))
        if k%3000==0:print('E8_CONTEXT_REBUILD',k,'changed',changed,flush=True)
    paths=[f'scripts/{TAG}.py',E.INITIAL,f'scripts/{BASE}.py']
    path.write_text(json.dumps(dict(status='LITERAL_E8_CONTEXT_HANDOFF_NOT_A_CERTIFICATE',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},vertex_cuts=vertices,steps=[],
        latest_numeric_weights=None,changed_profiles=changed),separators=(',',':'))+'\n')

def faces():
    obligations={};old=M.H.read('fresh_sqrt5_parameter_weight_stress_20261015')
    queue=M.H.read('fresh_sqrt5_nonsingular_rank29_fast_20261015')
    rows=[('prior',i,r['row'])for i,r in enumerate(old['faces'])if r['row']['nodes']<=23]
    rows += [('finite_queue',r['index'],r['row'])for r in queue['records']if r['result']['singular']and r['row']['nodes']<=23]
    oldcases=0
    for family,index,row in rows:
        es=[(a,b,m)for(a,b),m in zip(combinations(range(row['nodes']),2),row['labels'])]
        compiled=[OLD_COMPILE(row['white']+[a,b],es,[a,b],[None])for a,b in product((0,2),(1,3))]
        assert all(not any(v[3][None])for v in compiled)
        for ext in M.contexts(23-row['nodes']):
            p,tables=OLD_ADD(compiled,ext)
            assert not any(j in SPLIT for j in p)and all(not t[None]for t in tables)
            obligations.setdefault((M.key(p),F(1)),dict(kind='quad',family=family,index=index,external_counts=list(ext)))
            oldcases+=1
    tri=json.loads((ROOT/f'results/{A.TRI}_rank22.json').read_text())
    quad=json.loads((ROOT/f'results/{A.QUAD}.json').read_text())
    assert len(tri['records'])==996 and tri['complete_through_bound']and not tri['pending']
    todo=[('triangle',r['index'],r['row'],r['status'])for r in tri['records']]
    todo += [('quad',r['class_index'],r['row'],'SEALED_SINGULAR')for r in quad['immediate_sealed_weight_obligations']]
    todo += [('quad',r['class_index'],r['row'],'SINGULAR_EXTENSIBLE_STOP')for r in quad['records']if r['class_index']in(4446,4448)]
    counts=dict(inherited_baseline_zero_quad_rows=len(rows),inherited_quad_common_contexts=oldcases,
                triangle_rows=996,paid_quad_stopping_rows=37,triangle_marked_contexts=0,paid_quad_marked_contexts=0)
    for pos,(kind,index,row,status)in enumerate(todo):
        n=row['nodes'];white=row['white'];remaining=(22 if kind=='triangle'else 23)-n
        es=[(a,b,m)for(a,b),m in zip(combinations(range(n),2),row['labels'])]
        pairs=[[b for b in row['black']if b!=o]for o in row['black']]if kind=='triangle'else list(product((0,2),(1,3)))
        marks=white+([-1,-2,-3]if remaining else[])
        compiled=[compile_corner(white+list(pair),es,pair,marks)for pair in pairs]
        for ext in external_contexts(remaining):
            common,tables=add_bonus(compiled,ext)
            valid=white+([-1]if remaining-4*ext[3]-ext[0]>0 else[])+([-2]if ext[3]>0 else[])+([-3]if ext[0]>0 else[])
            for mark in valid:
                p=common
                for table in tables:p=M.plus(p,table[mark])
                obligations.setdefault((M.key(p),F(2 if kind=='triangle'else 1)),dict(kind=kind,
                    family='complete_triangles'if kind=='triangle'else'paid_quad_stopping',index=index,
                    producer_status=status,external_counts=list(ext),facet_root=mark))
                counts['triangle_marked_contexts'if kind=='triangle'else'paid_quad_marked_contexts']+=1
        if pos%100==0:print('E8_FACE_REBUILD',pos,'profiles',len(obligations),flush=True)
    obligations[M.key({0:F(-1)}),F(0)]=dict(kind='model_bound',name='distance2_discount_zero')
    obligations[M.key({-1:F(2,3),1:F(-1)}),F(0)]=dict(kind='model_bound',name='adjacent_discount_at_most_two_thirds')
    return [(dict(p),target,meta)for(p,target),meta in obligations.items()],counts

def barrier_release():
    old=json.loads((ROOT/'results/fresh_sqrt5_h4_triangle_compensation_barrier_20261016.json').read_text())
    records=[]
    for f in old['faces']:
        row=f['row'];mark=f['facet_root'];es=[(a,b,m)for(a,b),m in zip(combinations(range(row['nodes']),2),row['labels'])]
        compiled=[compile_corner(row['white']+[a,b],es,[a,b],[mark])for a,b in combinations(row['black'],2)]
        p,tables=add_bonus(compiled,tuple(f['external_counts'])+(0,))
        for t in tables:p=M.plus(p,t[mark])
        assert p=={int(j):F(z)for j,z in f['profile']}
        records.append(dict(index=f['index'],old_bonus=f['bonus_index'],E8_bonus=SPLIT[f['bonus_index']]))
    front=old['vertices'][1]['metadata'];p=literal_profile(front)
    for f in records:assert f['old_bonus']not in p and p[f['E8_bonus']]==3
    return records

def run(saved=None):
    install();A.TAG=TAG;A.INITIAL=INITIAL;A.CHECKPOINT=CHECKPOINT
    S.TAG=TAG;S.CHECKPOINT=CHECKPOINT;S.RESUME=True;S.checkpoint=A.checkpoint
    releases=barrier_release();obligations,counts=faces()
    print('MODEL',N,'variables',len(obligations),'profiles',counts,flush=True)
    cert,steps=S.solve(obligations)if saved is None else(saved['certificate'],saved['steps'])
    proof=M.verify(obligations,cert);result=json.loads(json.dumps(PRIOR));order=result['variable_order']
    for j in KEYS:order[j]=dict(order[j],non_H4_remainder_contains_E8=False)
    order += [dict(order[j],non_H4_remainder_contains_E8=True,inherited_E8_absent_index=j)for j in KEYS]
    assert len(order)==N
    paths=[f'scripts/{TAG}.py',INITIAL,f'scripts/{BASE}.py',f'results/{BASE}.json',E.INITIAL,
           'results/fresh_sqrt5_h4_triangle_compensation_barrier_20261016.json']
    result['sources'].update({p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths})
    result.update(status='PASS_EXACT_H4_E8_CONTEXT_PROPOSAL',variables=N,variable_order=order,
        face_counts=counts,unique_face_profiles=len(obligations),certificate=cert,verification=proof,steps=steps,
        added_E8_context_coefficients=len(KEYS),H4_triangle_barrier_release=releases,
        vertex_context_order=PRIOR['vertex_context_order']+['E8_present_if_H4_count1to3'],
        triangle_external_context_order=PRIOR['triangle_external_context_order']+['actual_E8_count'],
        quadrilateral_external_context_order=PRIOR['quadrilateral_external_context_order']+['actual_E8_count_for_new_paid_rows'],
        scope='All existing H4 compensation at feasible H4count1to3 now distinguishes actual E8 presence in the complete vertex forest.108new nonnegative coefficients; all ordinary H3 actual-count splits retained. Common external E8 counts are shared across every corner. Complete residual-forest DP has an E8-presence state, including marked E8 components. All996triangle states and37paidquadstops rebuilt; inherited zero-baseline quads are unchanged because their literal H4 compensation is zero. Full newly affected paid-quad extension coverage remains open. No dimension claim.')
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();install()
    if not args.verify:
        rebuild_initial()
        if not CHECKPOINT.exists():CHECKPOINT.write_bytes((ROOT/INITIAL).read_bytes())
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text())if args.verify else None
    result=json.loads(json.dumps(run(saved),default=str))
    if args.verify:assert result==saved
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    proof=result['verification'];number=proof.get('maximum',proof.get('lower_bound'))
    print(result['status'],result['certificate']['kind'],float(F(number)),flush=True)
