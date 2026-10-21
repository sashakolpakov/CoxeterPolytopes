#!/usr/bin/env python3
"""A parameter family for the existing rank13 contextual replacement.

The inherited arbitrary-support seals are unchanged. Their four local
totals and complete current output interfaces are checked for the entire
parameter interval by affine endpoint verification. No forest DP is run.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import fresh_rank13_contextual_coupon_rebalance_20261018 as R
import fresh_d14_a16_d15_joint_face_composition_20261018 as G
import fresh_d14_small_d15_ambient58_frontier_20261018 as C
import fresh_small_d15_fork_coupon_provenance_20261018 as I

ROOT=G.ROOT;B=G.B;S=G.S;T=R.T;A=R.A;Q=R.Q
TAG='fresh_rank13_variable_context_retention_20261018'
OLD=F(1,16);TARGET=F(9,64);ENDPOINTS=(OLD,F(1,2))
STEMS=tuple(dict.fromkeys((*R.STEMS,R.TAG,G.TAG,C.TAG,S.TAG,B.TAG,*B.STEMS,
 I.TAG,*(s for _,ss,_ in I.ROUTES for s in ss))))
enc=R.enc


def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def ordinary_coefficient(a,b):
    pair=frozenset((a,b))
    if pair in(frozenset((('A',2),('A',13))),frozenset((('A',2),('D',13)))):return -2
    return int(pair==frozenset((('A',1),('A',13))))
def actual_base(t):return F(C.actual_base(t),32)
def actual_cross(a,b,delta=TARGET):
    assert OLD<=delta<=F(1,2)
    return F(C.actual_cross(a,b),32)+(delta-OLD)*ordinary_coefficient(a,b)
def context(fs,delta=TARGET,upper=False):
    assert OLD<=delta<=F(1,2)
    co=Counter(fs);k=co['A',2]
    has10=any(n==10 for t,n in fs)if upper else co['A',10]>0
    return F(C.context(fs,upper),32)+(4*(delta-OLD)*k*(k-1)if has10 and co['A',1]>=2 else 0)
def value(fs,delta=TARGET):
    return sum(actual_base(t)for t in fs)+sum(actual_cross(a,b,delta)for i,a in enumerate(fs)for b in fs[i+1:])+context(fs,delta)


def coefficient(g,white,marks):
    out=R.changes(g,white,marks)
    x=out['delta']/OLD;assert x.denominator==1
    return dict(coefficient=int(x),roles=out['roles'],components=out['components'])


def seal_key(g,white):return json.dumps(enc([g,sorted(white)]),sort_keys=True)


def local(docs):
    inherited=R.local_seals(docs)
    assert enc(inherited)==docs[R.TAG]['local_seals']
    items=[('order16',inherited['order16'],inherited['order16']['gram'],inherited['order16']['white']),
           ('cycle14',inherited['cycle14'],inherited['cycle14']['gram'],docs[R.CY.TAG]['equality']['white'])]
    items += [('D15_'+str(r['index']),r,r['old_seal']['gram'],r['old_seal']['white'])for r in inherited['D15']]
    rows=[];negative_keys=set()
    for name,old,gg,white in items:
        g=A.matrix(gg);black=sorted({m for c in old['corners']for m in c['marks']})
        assert len(black)==4 and not set(black)&set(white)
        assert set(black)|set(white)==set(range(len(g)))
        assert Q.inertia(A.sub(g,white))==[len(white),0,0]
        assert all(Q.inertia(A.sub(g,white+[b]))==[len(white)+1,0,0]for b in black)
        entries=[]
        for c in old['corners']:
            ch=coefficient(g,white,c['marks']);before=F(c['payment'])
            values={str(t):before+(t-OLD)*ch['coefficient']for t in(*ENDPOINTS,TARGET)}
            assert all(v>=0 for v in values.values())
            entries.append(dict(marks=c['marks'],current_lower=before,**ch,parameter_values=values))
            if ch['coefficient']<0:negative_keys.add((seal_key(g,white),tuple(c['marks'])))
        totals={str(t):sum(c['parameter_values'][str(t)]for c in entries)for t in(*ENDPOINTS,TARGET)}
        assert all(v>=1 for v in totals.values())
        if name.startswith('D15_'):assert all(v==1+F(t)for t,v in totals.items())
        else:assert set(totals.values())=={F(1)}
        rows.append(dict(name=name,gram=g,white=white,corners=entries,totals=totals))
    assert len(rows)==4 and len(negative_keys)==5
    return rows,negative_keys


def literal(docs,negative_keys):
    old=docs[S.TAG]['literal_interfaces']['records'];rows=[];counts=Counter()
    for prev,(stem,case,g,white,cs)in zip(old,T.literal_inputs(docs),strict=True):
        assert(prev['source'],prev['case'])==(stem,enc(case))and len(prev['corners'])==len(cs)
        corners=[]
        for prior,(marks,*_)in zip(prev['corners'],cs):
            assert marks==prior['marks'];ch=coefficient(g,white,marks);counts.update(ch['roles'])
            if ch['coefficient']<0:assert(seal_key(g,white),tuple(marks))in negative_keys
            before=F(prior['retained_lower'])
            values={str(t):before+(t-OLD)*ch['coefficient']for t in(*ENDPOINTS,TARGET)}
            assert all(v>=0 for v in values.values())
            corners.append(dict(marks=marks,current_lower=before,**ch,parameter_values=values))
        totals={str(t):sum(c['parameter_values'][str(t)]for c in corners)for t in(*ENDPOINTS,TARGET)}
        assert all(v>=1 for v in totals.values())
        rows.append(dict(source=stem,case=case,corners=corners,totals=totals))
    assert len(rows)==216
    assert counts['reduce_A2_D13_middle_or_center']==2 and counts['reduce_A2_A13_middle_or_center']==0
    return dict(rows=rows,role_counts=dict(counts))


def later(docs,negative_keys):
    out=[];total=0;hits=[]
    for stem,keys in B.FINITE.items():
        rows=[]
        for section in keys:
            for path,c,gg,white,blocks in B.corners(docs[stem][section],section):
                g=A.matrix(gg);marks=c['marks'];ids=sorted(set(white)|set(marks))
                assert sorted(A.components(g,ids))==sorted(b['nodes']for b in blocks)
                ch=coefficient(g,white,marks)
                if ch['coefficient']<0:
                    assert(seal_key(g,white),tuple(marks))in negative_keys,(stem,path,ch)
                    hits.append(dict(source=stem,path=path,marks=marks,roles=ch['roles']))
                rows.append(dict(path=path,marks=marks,**ch));total+=1
        out.append(dict(source=stem,sections=list(keys),corners=rows))
    assert len(out)==22 and total==564
    # These 24 are newer than the baseline's original 564 finite interfaces.
    fresh=[]
    for p in docs[S.TAG]['selected_proofs']:
        for v in p.get('variants',[]):
            g=A.matrix(v['gram']);cs=[]
            for c in v['corners']:
                ch=coefficient(g,v['white'],c['marks']);assert ch['coefficient']>=0
                cs.append(dict(marks=c['marks'],**ch))
            fresh.append(dict(index=p['index'],extra_contact=v['extra_contact'],corners=cs))
    assert len(fresh)==24
    return dict(packets=out,corners=total,decreased_outputs=hits,new_D14_variants=fresh)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(p)==h,(s,p)
    assert R.DELTA==OLD and docs[R.TAG]['delta']=='1/16'
    local_rows,negative_keys=local(docs)
    li=literal(docs,negative_keys);late=later(docs,negative_keys)
    for d in(*ENDPOINTS,TARGET):
        for a in C.TYPES:
            for b in C.TYPES:assert 0<=actual_cross(a,b,d)==actual_cross(b,a,d)
    assert actual_cross(('A',2),('A',13))==actual_cross(('A',2),('D',13))==F(151,32)
    assert actual_cross(('A',1),('A',13))==F(137,64)
    a13=[('A',1)]+[('A',2)]*9+[('A',13)]*3
    quadratic=[('A',10)]+[('A',1)]*2+[('A',2)]*23
    examples=[]
    for name,fs in(('A13_obstruction',a13),('quadratic_context_test',quadratic)):
        assert sum(n for t,n in fs)==58
        values={str(t):value(fs,t)for t in(OLD,TARGET,F(5,32))}
        examples.append(dict(name=name,forest=fs,values=values))
    assert value(a13)==F(24117,64)<377
    assert value(quadratic,TARGET)==F(2821,8)<377<value(quadratic,F(5,32))==F(1537,4)
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in STEMS for folder,ext in
        [('scripts','py'),('results','json'),('results','md')]}
    return enc(dict(status='PASS_VARIABLE_RANK13_CONTEXT_FACE_RETENTION',
        sources={p:digest(p)for p in sorted(paths)},baseline=G.TAG,
        inherited_parameter=OLD,parameter_interval=list(ENDPOINTS),chosen_parameter=TARGET,
        rule=dict(A2_A13_middle='1/2-delta',A2_D13_center='1/2-delta',
            A1_A13_middle='1/2+delta',guarded_A2_A2_mark_pair='2*delta'),
        chosen_mark_weights=dict(A2_A13_middle='23/64',A2_D13_center='23/64',
            A1_A13_middle='41/64',guarded_A2_A2='9/32'),
        local_seals=local_rows,literal_interfaces=li,later_finite_interfaces=late,
        complete_actual_type_pairs=len(C.TYPES)**2,checked_parameter_values=list((*ENDPOINTS,TARGET)),
        illustrative_formal_forests=examples,triangles='2',quadrilaterals='1',
        same_component_weights_changed=False,ideal_weights_changed=False,
        other_context_rule_unchanged=True,forest_count_performed=False,new_dimension_bound=None,
        scope='Parameter extension of four completed arbitrary-support seals with complete current output inventory and216/564/24 exact interfaces. Affine endpoint checks cover the full delta interval. Two forest evaluations are diagnostics, not a complete changed count or dimension theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists()
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['chosen_parameter'],len(out['later_finite_interfaces']['decreased_outputs']),flush=True)
