#!/usr/bin/env python3
"""Independent SymPy-field cone clipping and polar-extension review."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sympy as s
import fresh_h9_h8_all_high_20261014_review as C

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_quartet_polar_extension_20261014_review'
BASE='results/fresh_h9_h8_disjoint_walk_20261014_exact.json'
POLAR='results/fresh_h9_quartet_polar_20261014.json'
EXT='results/fresh_h9_quartet_polar_extension_20261014.json'
SECOND='results/fresh_h9_quartet_polar_second_20261014.json'
RAD=6+s.Rational(16,5)*s.sqrt(5)
C.RR=s.sqrt(RAD);C.K=s.QQ.algebraic_field(s.sqrt(5),C.RR)
C.Z,C.O=C.K.zero,C.K.one
K,Z,O=C.K,C.Z,C.O
BASIS=[O,K.from_sympy(s.sqrt(5)),K.from_sympy(C.RR),K.from_sympy(s.sqrt(5)*C.RR)]
THETA=C.ev(K.ext.as_expr())

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
@lru_cache(None)
def decode4(t):return sum((K.convert(F(a))*b for a,b in zip(t,BASIS)),Z)
def dec(a):return decode4(tuple(a[0]+a[1]))
def dv(v):return list(map(dec,v))
def dm(v):return list(map(dv,v))
@lru_cache(None)
def sg(a):
    if not a:return 0
    out=C.I()
    for x in a.to_list():out=out*THETA+C.I(F(str(x)))
    assert out.l>0 or out.h<0,('unresolved interval',a)
    return 1 if out.l>0 else -1
def canon(v):return tuple(x/next(q for q in v if q)for x in v)
def encode(v):
    if isinstance(v,(list,tuple)):return[encode(x)for x in v]
    return list(map(str,v.to_list()))
def pair(v,w,Q):return C.pair(v,w,Q)
def pd(G):return all(sg(p)>0 for p in C.ldl(G))

def cone(normals,Q,H):
    """Clip a simplicial cone; detect edges by their complete face ray set.

    This does not call the producer's elimination or rank-adjacency engine.
    Coordinates are evaluations on ten independent starting walls. All ray
    coordinates stay nonpositive, so normalization always uses -first entry.
    """
    functionals=[C.mv(Q,n)for n in normals]
    seed=[]
    for i in range(len(normals)):
        candidate=seed+[i]
        if len(C.eliminate([functionals[j]for j in candidate])[1])==len(candidate):seed=candidate
        if len(seed)==10:break
    assert len(seed)==10
    BI=C.inverse([functionals[i]for i in seed])
    transformed=[[sum((row[k]*BI[k][j]for k in range(10)),Z)for j in range(10)]for row in functionals]
    def normalize(v):
        first=next(x for x in v if x);assert sg(first)<0
        return tuple(x/-first for x in v)
    rays=[]
    mask=sum(1<<i for i in seed)
    for j,i in enumerate(seed):
        v=[Z]*10;v[j]=-O;rays.append((tuple(v),mask^(1<<i)))
    stages=[]
    for i in range(len(normals)):
        if i in seed:continue
        values=[C.dot(v,transformed[i])for v,_ in rays]
        neg=[j for j,a in enumerate(values)if sg(a)<0]
        pos=[j for j,a in enumerate(values)if sg(a)>0]
        zero=[j for j,a in enumerate(values)if not a]
        new={rays[j][0]:rays[j][1]for j in neg}
        new.update({rays[j][0]:rays[j][1]|(1<<i)for j in zero})
        edges=0
        for j in neg:
            for k in pos:
                common=rays[j][1]&rays[k][1]
                if common.bit_count()<8:continue
                # The smallest common face is two-dimensional exactly when
                # these are its only extreme rays. This is independent of
                # the producer's cached rank-of-common-walls criterion.
                if any(h not in(j,k)and active&common==common for h,(_,active)in enumerate(rays)):continue
                edges+=1
                v=normalize([values[k]*a-values[j]*b for a,b in zip(rays[j][0],rays[k][0])])
                active=common|(1<<i)
                if v in new:assert new[v]==active
                new[v]=active
        stages.append(dict(wall=i,negative=len(neg),positive=len(pos),zero=len(zero),cut_edges=edges,rays=len(new)))
        rays=list(new.items())
        print('clip',len(normals),i,len(rays),flush=True)
    result={name:{}for name in('finite_vertices','cusps','hyperideal_rays','opposite_sheet_rays')}
    for y,active in rays:
        v=C.mv(BI,y);vals=[C.dot(row,v)for row in functionals]
        assert all(sg(a)<=0 for a in vals)
        ids=tuple(i for i,a in enumerate(vals)if not a)
        assert active==sum(1<<i for i in ids)
        assert len(C.eliminate([functionals[i]for i in ids])[1])==9
        q=pair(v,v,Q);time=pair(v,H,Q)
        kind='hyperideal_rays'if sg(q)>0 else'finite_vertices'if sg(q)<0 else'cusps'
        if sg(q)<=0 and sg(time)>=0:kind='opposite_sheet_rays'
        assert ids not in result[kind]
        result[kind][ids]=v
    return result,stages,seed

def relation(g,q,r):
    assert sg(g)<=0
    if not g:return 'right',2,None
    u=g*g/(q*r)
    if sg(u-O)>=0:return('parallel'if u==O else'dotted'),None,u
    allowed={3:K.convert(F(1,4)),4:K.convert(F(1,2)),5:(3*O+BASIS[1])/8,
             6:K.convert(F(3,4)),10:(5*O+BASIS[1])/8}
    for label,c in allowed.items():
        if u==c:return'finite',label,u
    return'forbidden',None,u

def run():
    seen={}
    def read(p):
        if p in seen:return seen[p]
        d=json.loads((ROOT/p).read_text());seen[p]=d
        for f,h in d.get('sources',{}).items():
            assert sha(f)==h,f
            if f.endswith('.json'):read(f)
        return d
    b,p,e,second=map(read,(BASE,POLAR,EXT,SECOND))
    assert p['status']=='EXACT_QUARTET_EXTERIOR_POLAR_CENSUS'
    assert e['status']=='PASS_EXACT_FOUR_POLAR_CAP_EXTENSION'
    assert second['status']=='PASS_EXACT_SECOND_POLAR_CAP_EXTENSION'
    assert len(K.mod.to_list())==5 and s.simplify(RAD*(6-s.Rational(16,5)*s.sqrt(5)))==-s.Rational(76,5)
    assert C.ev(RAD).l>0 and C.ev(6-s.Rational(16,5)*s.sqrt(5)).h<0
    Q,N,H=dm(second['form']),dm(second['normals']),dv(second['center'])
    assert len(N)==25 and Q==dm(b['form'])==dm(e['form'])and N[:19]==dm(b['normals'])and H==dv(b['center'])==dv(e['center'])
    assert N[:23]==dm(e['normals'])
    assert [sg(x)for x in C.ldl(Q)].count(-1)==1
    assert sg(pair(H,H,Q))<0 and all(sg(pair(H,n,Q))<0 for n in N)
    norms=[pair(n,n,Q)for n in N];assert norms==dv(second['norms'])and norms[:23]==dv(e['norms'])and all(sg(x)>0 for x in norms)
    G=[[pair(n,m,Q)for m in N]for n in N];assert G==dm(second['gram'])and[row[:23]for row in G[:23]]==dm(e['gram'])
    pc=Counter();pc23=Counter();newpairs={tuple(r['pair']):r for r in e['new_pairs']+second['new_pairs']}
    for i,j in combinations(range(25),2):
        typ,label,u=relation(G[i][j],norms[i],norms[j]);assert typ!='forbidden'
        pc[typ+(str(label)if label else'')]+=1
        if j<23:pc23[typ+(str(label)if label else'')]+=1
        if j>=19:
            row=newpairs.pop((i,j));assert(typ,label)==(row['kind'],row['label'])
            assert u==(None if row['invariant']is None else dec(row['invariant']))
    assert not newpairs and len(e['new_pairs'])==82 and len(second['new_pairs'])==47
    assert dict(pc23)==e['pair_counts']and dict(pc)==second['pair_counts']
    for number,doc in[(23,e),(25,second)]:
        for i,row in enumerate(doc['facet_points']):
            v=dv(row['point']);expected=[h-pair(H,N[i],Q)/norms[i]*n for h,n in zip(H,N[i])]
            assert v==expected and pair(v,v,Q)==dec(row['norm'])and sg(pair(v,v,Q))<0 and sg(pair(v,H,Q))<0
            assert all((not pair(v,n,Q))if j==i else sg(pair(v,n,Q))<0 for j,n in enumerate(N[:number]))
        assert len(doc['facet_points'])==number
    E=dv(second['strict_escape']['point']);assert E==dv(b['strict_escape']['point'])==dv(e['strict_escape']['point'])
    assert pair(E,E,Q)==dec(second['strict_escape']['norm'])==dec(e['strict_escape']['norm'])and sg(pair(E,E,Q))>0 and sg(pair(E,H,Q))<0
    assert [pair(E,n,Q)for n in N]==dv(second['strict_escape']['pairings'])
    assert [pair(E,n,Q)for n in N[:23]]==dv(e['strict_escape']['pairings'])
    assert all(sg(pair(E,n,Q))<0 for n in N)
    censuses=[];polar_results=[]
    for number,doc,records in[(19,p,p['polar_candidates']),(23,e,e['next_polar_candidates']),(25,second,second['next_polar_candidates'])]:
        current=N[:number];cn,stages,seed=cone(current,Q,H)
        counts={name:len(rows)for name,rows in cn.items()}
        assert counts==doc['counts']
        assert{key:sorted(map(list,rows))for key,rows in cn.items()}==doc['complete_census']
        accepted=[];statuses=Counter();functionals=[C.mv(Q,n)for n in current]
        assert len(records)==len(cn['hyperideal_rays'])
        for index,row in enumerate(records):
            assert row['index']==index
            ids=tuple(row['active']);v=dv(row['normal']);computed=cn['hyperideal_rays'][ids]
            assert canon(v)==canon(computed)
            q=pair(v,v,Q);vals=[C.dot(v,f)for f in functionals]
            assert q==dec(row['norm'])and sg(q)>0 and vals==dv(row['pairings'])
            assert tuple(j for j,x in enumerate(vals)if not x)==ids
            supp=row['independent_support'];assert set(supp)<=set(ids)and len(supp)==9
            assert len(C.eliminate([functionals[j]for j in supp])[1])==9
            bad=[]
            for j,g in enumerate(vals):
                typ,label,u=relation(g,q,norms[j]);saved=row['relations'][j]
                assert(typ,label)==(saved['kind'],saved['label'])
                assert u==(None if saved['invariant']is None else dec(saved['invariant']))
                if typ=='forbidden':bad.append(j)
            assert bad==row['forbidden_indices']and(not bad)==row['coxeter_compatible']
            assert pair(v,H,Q)==dec(row['center_pairing'])and sg(pair(v,H,Q))<0
            assert pair(v,E,Q)==dec(row['escape_pairing'])
            if not bad:accepted.append(index)
            statuses['compatible'if not bad else'forbidden']+=1
        assert accepted=={19:[37,38,95,186],23:[82,85],25:[62,95,96,128]}[number]
        if number==19:
            assert [dv(records[i]['normal'])for i in accepted]==N[19:23]
            assert set(cn['finite_vertices'])=={tuple(r['incident'])for r in b['complete_finite_vertex_census']['vertices']}
        else:
            saved={tuple(r['incident']):r for r in doc['finite_vertices']}
            assert set(saved)==set(cn['finite_vertices'])
            for ids,v in cn['finite_vertices'].items():
                row=saved[ids];sv=dv(row['point'])
                assert canon(v)==canon(sv)and len(ids)==9
                assert pd([[G[i][j]for j in ids]for i in ids])
                assert pair(sv,sv,Q)==dec(row['norm'])and [pair(sv,n,Q)for n in current]==dv(row['pairings'])
            vertices=doc['finite_vertices'];edges=[]
            for i,j in combinations(range(len(vertices)),2):
                support=sorted(set(vertices[i]['incident'])&set(vertices[j]['incident']))
                if len(support)==8:edges.append(dict(vertices=[i,j],support=support))
            assert edges==doc['finite_edges']and len(edges)=={23:213,25:247}[number]
            closed={(i,tuple(r['support']))for r in edges for i in r['vertices']}
            germs=[dict(vertex=i,support=[j for j in r['incident']if j!=release],released_wall=release)
                   for i,r in enumerate(vertices)for release in r['incident']
                   if(i,tuple(j for j in r['incident']if j!=release))not in closed]
            assert germs==doc['open_finite_germs']and len(germs)=={23:195,25:235}[number]
            old={tuple(v['incident'])for v in(b['complete_finite_vertex_census']['vertices']if number==23 else e['finite_vertices'])}
            assert old<set(saved)and len(saved)-len(old)==doc['new_finite_vertex_count']=={23:16,25:12}[number]
            if number==23:assert[dv(records[i]['normal'])for i in accepted]==N[23:]
            else:
                nextnormals=[dv(records[i]['normal'])for i in accepted]
                assert all(sg(pair(n,E,Q))<0 for n in nextnormals)
                assert all(relation(pair(a,b,Q),pair(a,a,Q),pair(b,b,Q))[0]=='dotted'for a,b in combinations(nextnormals,2))
        censuses.append(dict(walls=number,seed=seed,counts=counts,clipping_stages=stages))
        polar_results.append(dict(walls=number,tested=len(records),accepted=accepted,statuses=dict(statuses)))
    paths=[BASE,POLAR,EXT,SECOND,f'scripts/{TAG}.py','scripts/fresh_h9_h8_all_high_20261014_review.py',
           'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return dict(status='PASS_INDEPENDENT_QUARTET_POLAR_CONES_AND_EXTENSION',sources={p:sha(p)for p in paths},
      arithmetic='Independent SymPy algebraic number field, 192-bit outward rational intervals, and cone clipping with combinatorial face-adjacency; no producer or chamber-engine imports.',
      interval_bits=192,recursive_manifests=len(seen),cone_censuses=censuses,polar_tests=polar_results,
      pair_counts_23=dict(pc23),pair_counts_25=dict(pc),new_pairs=[82,47],actual_facets=[23,25],finite_vertices=[69,81],new_finite_vertices=[16,12],
      finite_edges=[213,247],open_finite_germs=[195,235],strict_escape_norm=encode(pair(E,E,Q)),
      all_label_alphabet=[2,3,4,5,6,10],
      scope='Complete independent 19/23/25 cone census; every285/413/513 exterior polar test; all300 final chamber pairs,25actual facets,81finite vertices,247edges,235open germs and strict escape. Allfour unadded next polars are mutually dotted and miss the escape. Scope is literal polar additions; arbitrary other caps/completions remain open. Both extended chambers have infinite volume.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['finite_vertices'],d['finite_edges'])
