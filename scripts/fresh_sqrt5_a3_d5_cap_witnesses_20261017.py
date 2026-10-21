#!/usr/bin/env python3
"""Two exact finite-extra cap witnesses; no optimizer or witness repair."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import fresh_sqrt5_multispectator_cap_capacity_20261017 as M
V,K=M.V,M.K
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_a3_d5_cap_witnesses_20261017'
MODEL='fresh_sqrt5_triple_groundfield_cap_joint_ports_20261017'


def gate(a,b,t):
    sg=M.signrad(a,b,t);sroot=M.ksqrt(t)
    if sroot is not None:a,b=a+b*sroot,K()
    inK=not b or not a
    squared=a*a+t*b*b
    gap=sg<0 and M.signrad(a+1,b,t)>0
    integral=inK and M.integral(squared)
    return dict(a=a.out(),b=b.out(),sign=sg,two_cycle_in_K=inK,
        integral_two_cycle=integral,forbidden_finite_angle_gap=gap,
        allowed=sg<=0 and not gap and integral)


def run():
    D={k:M.read(v)for k,v in M.TAGS.items()if k!='model'}
    S,I,A,Q,R=[D[k]for k in ('slices','incidence','actual','quad','raw')]
    model=M.read(MODEL);alpha={i:F(a)for i,a in model['face_multipliers']}
    profiles={x['base_type']:x for x in S['coface_slice_demands']}
    templates={x['canonical_core']:x for x in I['prism_templates']};out=[]
    for bi,finite_name,expected_core,capacity in [(4465,'A3',494,2),(5885,'D5',317,1)]:
        B=S['actual_bases'][bi];assert B['index']==bi and B['coface_location']==['finite','A1',0]
        ci=A['provenance_masses'][B['provenance']]['canonical_prism_core'];C=Q['canonical_cores'][ci]
        raw=R['records'][C['source_core']];g=V.raw_graph(raw)
        T=next(t for t in templates[ci]['bases']if t['deleted_normal']==B['deleted'])
        assert T['finite_extra']==['A1',finite_name] and T['base_core']==expected_core
        b,q=B['deleted'],B['coface_root'];mp=dict(T['core_isomorphism']);fin=sorted(set(range(raw['nodes']))-set(mp)-{b,q})
        assert T['opposite_base_facet']==q and all(g[q][i]==V.RIGHT for i in set(range(raw['nodes']))-{q,b})
        assert all(g[x][y]==V.RIGHT for x in fin for y in mp)
        H=[[K(2)if x==y else K(-1)if g[x][y][0]==3 else K()for y in fin]for x in fin]
        assert all(g[x][y][0]in(2,3)for x in fin for y in fin)
        J,det=V.inverse(H);assert V.pd(H)
        contacts=[int(g[b][x][0]==3)for x in fin];assert all(g[b][x][0]in(2,3)for x in fin)and sum(contacts)==1
        cg=V.core_graph(I['cores'][expected_core]);cc=[g[b][x]for x in sorted(mp,key=mp.get)]
        assert all(g[x][y]==cg[mp[x]][mp[y]]for x,y in combinations(mp,2))
        assert cc==[V.RIGHT,V.edge(3)]+[V.RIGHT]*(len(cg)-2)
        assert list(B['canonical_vertices'])==[0,1]
        aut=V.automorphisms(cg);assert all({p[i]for i in I['cores'][expected_core]['spherical_deletions']}=={0,1}for p in aut)
        rootorbit=sorted({p[1]for p in aut});assert rootorbit==[0,1]
        fg=[[V.RIGHT if i==j else g[x][y]for j,y in enumerate(fin)]for i,x in enumerate(fin)]
        froot=contacts.index(1);forbit=sorted({p[froot]for p in V.automorphisms(fg)})
        typ=I['base_types'][B['base_type']];counter=Counter(typ['spectator_forest']);assert counter[finite_name]==1
        k=counter['A1'];tau=V.squared(g[b][q]);h=2-tau/2-J[froot][froot]
        if finite_name=='A3':
            assert len(forbit)==2 and tau==K(F(3,2),F(1,2)) and h==K(F(1,2),F(-1,4)) and k==3
            root=M.ksqrt(V.squared(cg[0][1]));assert root==K(1,1)
            G=[[K(2)if i==j else -M.ksqrt(V.squared(cg[i][j]))if cg[i][j]!=V.RIGHT else K()for j in range(4)]for i in range(4)]
            L,det=V.inverse(G);assert det.sign()<0 and L[0][0]==L[1][1]==h and L[0][1]==K(F(-1,4))
            assert all(J[i][i]==K(F(3,4))for i in forbit)and J[forbit[0]][forbit[1]]==K(F(1,4))
            corepairs=lambda i,j:(L[i][j],K())
        else:
            assert len(forbit)==1 and tau==K(3,1)and h==K(F(-1,2),F(-1,2))and k==2
            assert V.squared(cg[0][1])==tau and M.ksqrt(tau) is None and J[froot][froot]==1
            assert K(2)/(4-tau)==h
            corepairs=lambda i,j:(h,K())if i==j else(K(),K(1)/(4-tau))
        positions=list(product(rootorbit,range(k),forbit));pairs=[];adj=[0]*len(positions)
        for i,j in combinations(range(len(positions)),2):
            o,q0,f0=positions[i];p,q1,f1=positions[j];a,z=corepairs(o,p)
            a+=J[f0][f1]+(tau/2 if q0==q1 else K())
            r=gate(a,z,tau);pairs.append(dict(pair=[i,j],**r))
            if r['allowed']:adj[i]|=1<<j;adj[j]|=1<<i
        memo={0:(0,())}
        def best(mask):
            if mask in memo:return memo[mask]
            bit=mask&-mask;i=bit.bit_length()-1;rest=mask^bit
            a,x=best(rest);b,y=best(rest&adj[i]);b+=1
            memo[mask]=(b,(i,)+y)if b>a else(a,x);return memo[mask]
        maximum,witness=best((1<<len(positions))-1);assert maximum==capacity
        for mask,(val,_)in memo.items():
            if mask:
                bit=mask&-mask;i=bit.bit_length()-1;rest=mask^bit
                assert val==max(memo[rest][0],1+memo[rest&adj[i]][0])
        pf=profiles[B['base_type']];load=alpha[B['face']];mass=alpha[pf['simplex_face']]
        assert load>capacity*mass
        out.append(dict(base_index=bi,base_type=B['base_type'],provenance=B['provenance'],Qface=B['face'],
            simplex_face=pf['simplex_face'],core=expected_core,raw_source=C['source_core'],core_isomorphism=T['core_isomorphism'],
            finite_extra=finite_name,finite_roots=fin,finite_gram=[[x.out()for x in row]for row in H],
            finite_inverse=[[x.out()for x in row]for row in J],finite_contact=contacts,finite_profile_roots=forbit,
            core_orientations=rootorbit,core_automorphisms=aut,spectator_forest=typ['spectator_forest'],A1_count=k,
            tau=tau.out(),core_projection_norm=h.out(),positions=positions,pairs=pairs,capacity=capacity,
            witness=witness,recurrence=[[m,v,list(w)]for m,(v,w)in sorted(memo.items())],
            positive_face_coefficients=[[B['face'],1]],load=str(load),simplex_mass=str(mass),slack=str(capacity*mass-load)))
    paths=[Path(__file__),Path(M.__file__),Path(V.__file__),ROOT/f'results/{MODEL}.json']
    paths += [ROOT/f'results/{tag}.json'for key,tag in M.TAGS.items()if key!='model']
    return dict(status='PASS_TWO_EXACT_A3_D5_CAP_WITNESSES',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        ground_field='Q(sqrt5)',witnesses=out,input_bound=model['bound'],new_bound=None,optimizer_invocations=0,
        scope='Two selected actual cap-profile families at fixed full simplex stars. Complete identical spectator placements and core/finite contact orientations checked. Capacity2 for base110693 and capacity1 for base8351. Exact local necessary compatibility, not a complete79-family census, repaired witness, ambient realization or dimension theorem.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],[(r['base_type'],r['capacity'],r['slack'])for r in d['witnesses']],flush=True)
