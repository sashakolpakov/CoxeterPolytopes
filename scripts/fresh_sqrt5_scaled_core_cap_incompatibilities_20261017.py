#!/usr/bin/env python3
"""Two complete selected cap families discovered through exact scaled cores.

No optimization or witness repair.  All represented copies of the three
literal raw profiles are matched to their actual provenance and maximal edge.
"""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import fresh_sqrt5_multispectator_cap_capacity_20261017 as M

V,K=M.V,M.K
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_scaled_core_cap_incompatibilities_20261017'
MODEL='fresh_sqrt5_core23_25_spectator_incidence_20261017'
FULL='fresh_sqrt5_cross_coface_cap_full_model_20261017'
FAMILIES={555:(254,),371:(259,262)}


def sub(G,ids):return [[G[i][j]for j in ids]for i in ids]
def out(G):return [[z.out()for z in row]for row in G]
def dot(p,J,q):return sum((p[i]*J[i][j]*q[j]for i in range(len(p))for j in range(len(q))),K())
def conjugate(G):return [[K(z.a,-z.b)for z in row]for row in G]


def gauge(cg):
    n=len(cg);roots={}
    for i,j in combinations(range(n),2):
        if cg[i][j]!=V.RIGHT:roots[i,j]=M.ksqrt(V.squared(cg[i][j]))
    if all(z is not None for z in roots.values()):
        scales=[K(1)]*n;G=[[K(2)if i==j else K()for j in range(n)]for i in range(n)]
        for (i,j),z in roots.items():G[i][j]=G[j][i]=-z
    else:
        scales={0:K(1)};todo=[0]
        for i in todo:
            for j in range(n):
                if i!=j and cg[i][j]!=V.RIGHT and j not in scales:
                    scales[j]=scales[i]*V.squared(cg[i][j]);todo.append(j)
        assert len(scales)==n;scales=[scales[i]for i in range(n)]
        G=[[2*scales[i]if i==j else K()for j in range(n)]for i in range(n)]
        for i,j in roots:
            z=M.ksqrt(scales[i]*scales[j]*V.squared(cg[i][j]));assert z is not None
            G[i][j]=G[j][i]=-z
    assert all(s.sign()>0 and K(s.a,-s.b).sign()>0 for s in scales)
    for i,j in combinations(range(n),2):
        assert G[i][j].sign()<=0 and G[i][j]*G[i][j]/(scales[i]*scales[j])==V.squared(cg[i][j])
    return G,scales


def cap_gauge(es,scales):
    ports=[i for i,e in enumerate(es)if e!=V.RIGHT];assert ports
    s=scales[ports[0]]*V.squared(es[ports[0]])
    assert s.sign()>0 and K(s.a,-s.b).sign()>0
    p=[K()for _ in es]
    for i in ports:
        z=M.ksqrt(s*scales[i]*V.squared(es[i]));assert z is not None
        p[i]=-z
    assert all(p[i]*p[i]/(s*scales[i])==V.squared(es[i])for i in range(len(es)))
    return s,p


def orientation_orbit(es,vertices,aut):
    ans=set()
    for perm in aut:
        image=[V.RIGHT for _ in es]
        for i,e in enumerate(es):image[perm[i]]=e
        ans.add((tuple(image),tuple(sorted(perm[i]for i in vertices))))
    return sorted(ans)


def run():
    D={key:M.read(tag)for key,tag in M.TAGS.items()if key!='model'}
    S,I,A,Q,R=[D[k]for k in('slices','incidence','actual','quad','raw')]
    model,full=M.read(MODEL),M.read(FULL)
    alpha={i:F(a)for i,a in model['face_multipliers']};assert all(a>=0 for a in alpha.values())
    pmap=dict(full['prepared']['provenance_to_face']);assert len(set(pmap.values()))==len(pmap)
    sfaces={}
    for p in S['coface_slice_demands']:
        if p['base_type']in sfaces:assert sfaces[p['base_type']]==p['simplex_face']
        sfaces[p['base_type']]=p['simplex_face']
    templates={t['canonical_core']:t for t in I['prism_templates']}
    contexts={r[0]:r for r in I['contexts']};assert len(contexts)==len(I['contexts'])
    groups=defaultdict(list);occurrences=[];seed={};coredata={};seen=set()
    for B in S['actual_bases']:
        t=B['base_type'];typ=I['base_types'][t];core=typ['core']
        if core not in FAMILIES or B['coface_location']!=['finite','A1',0]:continue
        prov=A['provenance_masses'][B['provenance']];ci=prov['canonical_prism_core'];C=Q['canonical_cores'][ci]
        raw=C['source_core']
        if raw not in FAMILIES[core]:continue
        T=next(x for x in templates[ci]['bases']if x['deleted_normal']==B['deleted'])
        if T['finite_extra']!=['A1']:continue
        b,q=B['deleted'],B['coface_root'];mp=dict(T['core_isomorphism'])
        r=R['records'][raw];g=V.raw_graph(r);cg=V.core_graph(I['cores'][core]);n=len(cg)
        assert set(range(r['nodes']))==set(mp)|{b,q} and len(mp)==n
        assert sorted(mp.values())==list(range(n))
        assert all(g[x][y]==cg[mp[x]][mp[y]]for x,y in combinations(mp,2))
        assert all(g[q][x]==V.RIGHT for x in mp)
        assert T['base_core']==core and T['opposite_base_facet']==q
        assert typ['simplex_dimension']==T['simplex_dimension']==1
        assert B['prism_dimension']==prov['prism_dimension']==2
        assert B['canonical_vertices']==I['cores'][core]['spherical_deletions']==[0,1]
        assert sorted(mp[x]for x in T['full_spherical_deletions'])==[0,1]
        assert Counter(typ['spectator_forest'])['A1']>=1
        assert n+sum(map(V.rank,typ['spectator_forest']))==14
        assert B['face']==pmap[B['provenance']]
        context=contexts[prov['source_context']]
        assert context[2:4]==[ci,2] and dict(context[4])[t]==1
        assert (B['provenance'],b)not in seen;seen.add((B['provenance'],b))
        vertices=set()
        for left,right in C['quad_partitions']:
            if b in left:side,other=left,right
            elif b in right:side,other=right,left
            else:continue
            assert set(side)=={b,q};vertices.update(other)
        assert sorted(mp[x]for x in vertices)==[0,1]
        es=tuple(g[b][x]for x in sorted(mp,key=mp.get));tau=V.squared(g[b][q])
        assert tau.sign()>0
        if core not in coredata:
            G,scales=gauge(cg);J,det=V.inverse(G);aut=V.automorphisms(cg)
            assert det.sign()<0 and V.pd(conjugate(G))
            assert all({perm[i]for i in(0,1)}=={0,1}for perm in aut)
            assert [i for i in range(n)if V.pd(sub(G,[j for j in range(n)if j!=i]))]==[0,1]
            coredata[core]=(G,scales,J,aut)
        G,scales,J,aut=coredata[core];orbit=orientation_orbit(es,[0,1],aut)
        if raw in seed:assert seed[raw]['orbit']==orbit and seed[raw]['tau']==tau
        else:seed[raw]=dict(core=core,tau=tau,orbit=orbit)
        rec=dict(index=len(occurrences),base_index=B['index'],base_type=t,simplex_face=sfaces[t],
            provenance=B['provenance'],face=B['face'],raw_source=raw,canonical_prism_core=ci,
            source_context=prov['source_context'],
            deleted=b,coface=q,core=core,core_isomorphism=T['core_isomorphism'],
            canonical_base_vertices=[0,1],spectator_forest=typ['spectator_forest'])
        occurrences.append(rec);groups[core,t].append(rec)
    assert set(seed)=={254,259,262}
    families=[]
    for core,raws in sorted(FAMILIES.items()):
        G,scales,J,aut=coredata[core];orients=[]
        for raw in raws:
            for es,vertices in seed[raw]['orbit']:
                s,p=cap_gauge(es,scales);tau=seed[raw]['tau']
                assert dot(p,J,p)+s*tau/2==2*s
                orients.append(dict(raw_source=raw,contacts=es,vertices=vertices,scale=s,covector=p,tau=tau))
        tests=[]
        for ia,ib in product(range(len(orients)),repeat=2):
            a,b=orients[ia],orients[ib];h=dot(a['covector'],J,b['covector'])
            sq=h*h/(a['scale']*b['scale'])
            same_q_sign=M.signrad(h,K(F(1,2)),a['scale']*b['scale']*a['tau']*b['tau'])
            assert same_q_sign>0
            if core==555:
                assert h.sign()>0
                expected=K(F(5,4),F(-1,4))if ia==ib else K(F(3,4),F(-1,4))
                assert a['scale']==b['scale']==1 and h==expected
                reason='positive_product'
            elif not M.integral(sq):reason='nonintegral_two_cycle'
            else:
                # Every remaining product lies strictly between -1 and 0.
                assert h.sign()<0 and sq.sign()>0 and (sq-1).sign()<0
                reason='open_coxeter_angle_gap'
            tests.append(dict(orientations=[ia,ib],raw_sources=[a['raw_source'],b['raw_source']],
                scaled_product=h.out(),normalized_two_cycle=sq.out(),different_q_forbidden_reason=reason,
                same_q_product=dict(rational_part=h.out(),radical_coefficient='1/2',
                    radicand=(a['scale']*b['scale']*a['tau']*b['tau']).out(),sign=same_q_sign)))
        if core==555:assert len(aut)==2 and len(orients)==2 and seed[254]['tau']==K(F(3,2),F(1,2))
        else:
            assert len(orients)==2 and all(r['vertices']==(0,1)for r in orients)
            assert tests[1]['normalized_two_cycle']==K(F(1,2)).out()
        families.append(dict(core=core,raw_sources=list(raws),capacity=1,
            core_gram=out(G),core_inverse=out(J),scales=[x.out()for x in scales],automorphisms=aut,
            orientations=[dict(raw_source=o['raw_source'],contacts=o['contacts'],base_vertices=o['vertices'],
                cap_scale=o['scale'].out(),covector=[x.out()for x in o['covector']],tau=o['tau'].out())for o in orients],
            ordered_pair_tests=tests))
    columns=[];violations=[];required={}
    for (core,t),rs in sorted(groups.items()):
        si=sfaces[t];positive=Counter(r['face']for r in rs)
        load=sum((w*alpha.get(f,F())for f,w in positive.items()),F());mass=alpha.get(si,F())
        rec=dict(index=len(columns),core=core,base_type=t,simplex_face=si,capacity=1,
            source_bases=[r['base_index']for r in rs],occurrences=[r['index']for r in rs],
            positive_face_coefficients=sorted(positive.items()),
            provenance_coefficients=sorted(Counter(r['provenance']for r in rs).items()),
            spectator_forest=I['base_types'][t]['spectator_forest'],load=str(load),
            initial_simplex_mass=str(mass),gap=str(load-mass))
        columns.append(rec)
        if load>mass:violations.append(rec['index']);required[si]=max(required.get(si,mass),load)
    assert len(violations)==2
    assert [(columns[i]['base_type'],F(columns[i]['gap']))for i in violations]==[(48730,F(9,2500000000)),(114227,F(2991,200000000))]
    paths=[Path(__file__),Path(M.__file__),Path(V.__file__)]
    paths +=[ROOT/f'results/{t}.json'for k,t in M.TAGS.items()if k!='model']
    paths +=[ROOT/f'results/{t}.json'for t in(MODEL,FULL)]
    return dict(status='PASS_EXACT_SCALED_CORE_CAP_INCOMPATIBILITY_FAMILIES',dimension=13,ground_field='Q(sqrt5)',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        families=families,occurrences=occurrences,columns=columns,violating_columns=violations,
        required_simplex_masses=[[si,str(alpha.get(si,F())),str(v)]for si,v in sorted(required.items())],
        counts=dict(families=len(families),literal_profiles=3,actual_base_occurrences=len(occurrences),
            columns=len(columns),positive_columns=sum(F(c['load'])>0 for c in columns),violations=len(violations)),
        input_bound=model['bound'],alpha_changed=False,mu_changed=False,optimizer_invocations=0,
        scope='Two complete selected raw-profile families over all represented actual-provenance copies. Capacity one holds across every isolated A1 coface and every colored core orientation. Raw254 is excluded by acuteness alone; raw259/262 also use necessary two-cycle integrality or a finite-angle gap. No repaired witness, optimality, full geometric coverage or dimension theorem is asserted.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));f=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(f.read_text())
    else:f.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],d['required_simplex_masses'],flush=True)
