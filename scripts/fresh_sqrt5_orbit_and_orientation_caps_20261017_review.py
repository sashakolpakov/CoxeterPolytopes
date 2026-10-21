#!/usr/bin/env python3
"""Independent phi-basis review of orbit and cross-orientation cap capacities.

No producer/project imports. Raw incidence is inherited but literal maps,
all selected core automorphisms/PD deletions, projections and loads are rebuilt.
"""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_orbit_and_orientation_caps_20261017_review'
TAGS=dict(orbit='fresh_sqrt5_orbit_integral_cap_capacity_20261017',
 cross='fresh_sqrt5_cross_orientation_cap_capacity_20261017',
 slices='fresh_sqrt5_labelled_coface_base_packing_20261017',
 incidence='fresh_sqrt5_prism_simplex_base_incidence_20261017',
 actual='fresh_sqrt5_actual_provenance_dual_20261017',
 quad='fresh_sqrt5_rank13_packing_quad_separation_20261017',
 raw='fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017',
 previous='fresh_sqrt5_cross_coface_cap_capacity_20261017',
 model='fresh_sqrt5_cross_coface_cap_full_model_20261017',
 cliques='fresh_sqrt5_cap_profile_cliques_20261017')

class Phi:
    """u+v*phi, phi^2=phi+1; independent of producer sqrt5 engine."""
    def __init__(self,u=0,v=0):
        if isinstance(u,Phi):self.u,self.v=u.u,u.v
        else:self.u,self.v=F(u),F(v)
    @classmethod
    def radical(cls,z):
        a,b=map(F,z);return cls(a-b,2*b)
    def __add__(self,x):
        x=Phi(x);return Phi(self.u+x.u,self.v+x.v)
    __radd__=__add__
    def __neg__(self):return Phi(-self.u,-self.v)
    def __sub__(self,x):return self+-Phi(x)
    def __rsub__(self,x):return Phi(x)+-self
    def __mul__(self,x):
        x=Phi(x);return Phi(self.u*x.u+self.v*x.v,self.u*x.v+self.v*x.u+self.v*x.v)
    __rmul__=__mul__
    def __truediv__(self,x):
        x=Phi(x);d=x.u*x.u+x.u*x.v-x.v*x.v
        assert d;return self*Phi((x.u+x.v)/d,-x.v/d)
    def __rtruediv__(self,x):return Phi(x)/self
    def __bool__(self):return bool(self.u or self.v)
    def __eq__(self,x):
        x=Phi(x);return self.u==x.u and self.v==x.v
    def sign(self):
        a,b=self.u+self.v/2,self.v/2
        if not b:return (a>0)-(a<0)
        if not a:return (b>0)-(b<0)
        if a*b>0:return (a>0)-(a<0)
        d=a*a-5*b*b;return ((a>0)-(a<0))*((d>0)-(d<0))
    def out(self):return [str(self.u+self.v/2),str(self.v/2)]
    def integral(self):
        trace=2*self.u+self.v;norm=self.u*self.u+self.u*self.v-self.v*self.v
        test=self.u.denominator==self.v.denominator==1
        assert test==(trace.denominator==norm.denominator==1)
        return test
    def trace_norm(self):return [str(2*self.u+self.v),str(self.u*self.u+self.u*self.v-self.v*self.v)]


def invert(A):
    """Reverse-column Gaussian elimination, then independent back substitution."""
    n=len(A);order=list(reversed(range(n)))
    M=[[Phi(A[i][j]) for j in order]+[Phi(i==j) for j in range(n)] for i in order]
    det=Phi(1)
    for k in range(n):
        p=next(i for i in range(k,n) if M[i][k])
        if p!=k:M[p],M[k]=M[k],M[p];det=-det
        d=M[k][k];det*=d
        for i in range(k+1,n):
            c=M[i][k]/d
            for j in range(k,n+n):M[i][j]-=c*M[k][j]
    Z=[[Phi() for _ in range(n)] for _ in range(n)]
    for rhs in range(n):
        for i in reversed(range(n)):
            Z[i][rhs]=(M[i][n+rhs]-sum((M[i][j]*Z[j][rhs] for j in range(i+1,n)),Phi()))/M[i][i]
    J=[None]*n
    for i,j in enumerate(order):J[j]=Z[i]
    assert all(sum((A[i][k]*J[k][j] for k in range(n)),Phi())==int(i==j)
               for i in range(n) for j in range(n))
    return J,det


def positive(A):
    M=[[Phi(x) for x in row] for row in A]
    while M:
        d=M[-1][-1]
        if d.sign()<=0:return False
        M=[[M[i][j]-M[i][-1]*M[-1][j]/d for j in range(len(M)-1)] for i in range(len(M)-1)]
    return True


def edge(m,a='0',b='0'):
    return (int(m),str(F(a)),str(F(b))) if m==0 else (int(m),'0','0')
R=edge(2)

def graph(c):
    n=c['nodes'];G=[[R]*n for _ in range(n)]
    for i,j,e in c['edges']:G[i][j]=G[j][i]=edge(*e)
    return G

def raw_graph(r):
    n=r['nodes'];G=[[R]*n for _ in range(n)]
    for (i,j),m in zip(combinations(range(n),2),r['row']['labels']):
        G[i][j]=G[j][i]=edge(m,*r['parameters'][f'{i},{j}']) if m==0 else edge(m)
    return G

def symmetries(G):
    n=len(G);degrees=[Counter(G[i][j] for j in range(n) if j!=i) for i in range(n)]
    out=[]
    def fill(mp,used):
        if len(mp)==n:out.append(tuple(mp[i] for i in range(n)));return
        choices={i:[j for j in range(n) if j not in used and degrees[i]==degrees[j]
             and all(G[i][a]==G[j][b] for a,b in mp.items())] for i in range(n) if i not in mp}
        i=min(choices,key=lambda i:(len(choices[i]),-i))
        for j in choices[i]:fill({**mp,i:j},used|{j})
    fill({},set());return sorted(out)

def gram(G):
    n=len(G);H=[[Phi() for _ in range(n)] for _ in range(n)]
    if all(G[i][j][0] in (2,3,5) for i in range(n) for j in range(i)):
        for i in range(n):
            for j in range(n):H[i][j]=Phi(2) if i==j else {2:Phi(),3:Phi(-1),5:Phi(0,-1)}[G[i][j][0]]
        return H,None
    # Opposite rooting from the older checker. Positive gauge handles radicals.
    E=[(i,j) for i in range(n) for j in range(i) if G[i][j]!=R]
    assert len(E)==n-1
    s={n-1:Phi(1)};todo=[n-1]
    for i in todo:
        for j in reversed(range(n)):
            if G[i][j]==R or j in s:continue
            m,a,b=G[i][j]
            w=Phi.radical((a,b)) if m==0 else {3:Phi(1),4:Phi(2),5:Phi(1,1),6:Phi(3),10:Phi(2,1)}[m]
            assert w.sign()>0
            s[j]=s[i]*w;H[i][j]=H[j][i]=-s[j];todo.append(j)
    assert len(s)==n
    for i in range(n):H[i][i]=2*s[i]
    return H,[s[i] for i in range(n)]

def rank(s):return 2 if s.startswith('I2(') else int(s[1:])
def clean(x):return json.loads(json.dumps(x))

def run():
    D={k:json.loads((ROOT/f'results/{v}.json').read_text()) for k,v in TAGS.items()}
    P,X,S,I,A,Q,W,O,M=[D[k] for k in ('orbit','cross','slices','incidence','actual','quad','raw','previous','model')]
    assert P['status']=='PASS_EXACT_ORBIT_AND_INTEGRALITY_CAP_CAPACITIES'
    assert X['status']=='PASS_EXACT_CROSS_ORIENTATION_CAP_CAPACITIES'
    dual=next(c for c in M['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={i:F(a) for i,a in dual['face_multipliers']};assert all(a>=0 for a in alpha.values())
    profiles={p['base_type']:p for p in S['coface_slice_demands']}
    templates={t['canonical_core']:t for t in I['prism_templates']}
    savedcores={c['core']:c for c in S['core_symmetry_and_subset_orbits']}
    saved={o['base_index']:o for o in P['occurrences']};assert len(saved)==len(P['occurrences'])
    cache={};checks=[];groups=defaultdict(list);seen=set();selected=[];candidate_counts=Counter()
    for B in S['actual_bases']:
        if B['coface_location']!=['finite','A1',0]:continue
        candidate_counts['A1_base_occurrences']+=1
        t=B['base_type'];p=profiles[t];core=p['core'];typ=I['base_types'][t]
        j,b,q=B['provenance'],B['deleted'],B['coface_root'];prov=A['provenance_masses'][j]
        ci=prov['canonical_prism_core'];C=Q['canonical_cores'][ci];r=W['records'][C['source_core']]
        T=next(z for z in templates[ci]['bases'] if z['deleted_normal']==b);mp=dict(T['core_isomorphism'])
        if set(mp)!=(set(range(r['nodes']))-{b,q}):continue
        candidate_counts['sole_A1_extra']+=1
        G=raw_graph(r)
        if G[b][q][0]!=0:continue
        tau=Phi.radical(G[b][q][1:]);assert (tau-4).sign()>0 and tau.integral()
        h=2-tau/2;w=h*h
        finite=(h+2).sign()>0 and h not in (Phi(-1),Phi(0,-1))
        nonintegral=not w.integral()
        if not(finite or nonintegral):continue
        selected.append(B['index']);o=saved[B['index']]
        assert (j,b) not in seen;seen.add((j,b))
        assert b!=q and b not in mp and q not in mp
        assert T['opposite_base_facet']==q and T['base_core']==core==typ['core']
        assert T['finite_extra']==['A1'] and Counter(typ['spectator_forest'])['A1']>=1
        cg=graph(I['cores'][core]);n=len(cg)
        assert sorted(mp.values())==list(range(n))
        assert all(G[q][i]==R for i in mp)
        assert all(G[i][k]==cg[mp[i]][mp[k]] for i,k in combinations(mp,2))
        assert n+sum(map(rank,typ['spectator_forest']))==14
        assert r['nodes']+sum(map(rank,typ['spectator_forest']))-1==15
        vertices=set();pairs=set()
        for left,right in C['quad_partitions']:
            if b in left:side,other=left,right
            elif b in right:side,other=right,left
            else:continue
            assert set(side)=={b,q};vertices.update(other);pairs.update(tuple(sorted((b,k))) for k in other)
        assert vertices and vertices<=set(mp)
        assert pairs<={tuple(sorted(z)) for z in C['pd_deletions']}
        vb=tuple(sorted(mp[k] for k in vertices))
        assert list(vb)==B['canonical_vertices']==o['base_vertices']
        assert len(vb)==B['prism_dimension']==templates[ci]['prism_dimension']
        if core not in cache:
            aut=symmetries(cg);H,scales=gram(cg);J,det=invert(H);assert det.sign()<0
            pdroots=[k for k in range(n) if positive([[H[i][j] for j in range(n) if j!=k] for i in range(n) if i!=k])]
            assert pdroots==I['cores'][core]['spherical_deletions'] and pdroots
            assert aut==sorted(map(tuple,savedcores[core]['automorphisms']))
            assert all({a[i] for i in pdroots}==set(pdroots) for a in aut)
            cache[core]=(aut,H,scales,J,det,pdroots)
        aut,H,scales,J,det,pdroots=cache[core]
        assert set(vb)<=set(pdroots)
        assert sorted(mp[k] for k in T['full_spherical_deletions'])==pdroots
        contacts=tuple(G[b][i] for i in sorted(mp,key=mp.get))
        assert clean(contacts)==o['contacts'] and all(e[0] in (2,3) for e in contacts)
        ports=[i for i,e in enumerate(contacts) if e[0]==3]
        if scales is None:norm=sum((J[i][k] for i in ports for k in ports),Phi())
        else:
            assert len(ports)==1;norm=scales[ports[0]]*J[ports[0]][ports[0]]
        assert norm==h and h.sign()<0
        orbit=set()
        for a in aut:
            moved=[None]*n
            for k in range(n):moved[a[k]]=contacts[k]
            orbit.add((tuple(moved),tuple(sorted(a[k] for k in vb))))
        orbit=tuple(sorted(orbit));capacity=len(orbit)
        assert capacity in (1,2) and len(aut)%capacity==0
        assert clean(orbit)==o['orientation_orbit'] and capacity==o['capacity']
        assert tau.out()==o['dotted_squared_product'] and h.out()==o['repeated_wall_product'] and w.out()==o['two_cycle']
        assert w.trace_norm()==[o['two_cycle_trace'],o['two_cycle_norm']]
        assert finite==o['forbidden_finite_angle'] and nonintegral==o['nonintegral_two_cycle']
        assert (h+1).sign()>0 if o['old_open_gap'] else (h+1).sign()<=0
        assert o['arithmetic_needed']==(not finite)
        assert (j,B['face'],ci,C['source_core'],b,q,core)==tuple(o[k] for k in ('provenance','face','canonical_prism_core','source_core','deleted_normal','coface_root','core'))
        assert o['core_isomorphism']==sorted(map(list,mp.items())) and o['automorphism_count']==len(aut)
        key=(t,orbit,tuple(tau.out()));groups[key].append((B,p['simplex_face']))
        checks.append(dict(base_index=B['index'],core=core,capacity=capacity,canonical_vertices=list(vb),raw_omitted_pairs=sorted(map(list,pairs)),projection_norm=h.out(),two_cycle=w.out(),finite_gap=finite,nonintegral_cycle=nonintegral))
    assert sorted(selected)==sorted(saved) and len(selected)==815
    assert dict(candidate_counts)=={k:P['counts'][k] for k in candidate_counts}
    expected={tuple(c['source_bases']):c for c in P['columns']};xcols={c['index']:c for c in X['columns']}
    pairchecks={z['column']:z for z in X['pair_checks']};loads=[];viol=[];xviol=[];tightened=[];crosschecks=[]
    required=defaultdict(F);xrequired=defaultdict(F)
    for (t,orbit,tau),rows in sorted(groups.items()):
        ids=tuple(sorted(b['index'] for b,_ in rows));c=expected[ids];xc=xcols[c['index']]
        si=rows[0][1];assert all(s==si for _,s in rows)
        co=Counter(b['face'] for b,_ in rows);pr=Counter(b['provenance'] for b,_ in rows)
        cap=len(orbit);load=sum((m*alpha.get(i,F()) for i,m in co.items()),F());mass=alpha.get(si,F())
        assert c['simplex_face']==si and c['capacity']==cap
        assert c['positive_face_coefficients']==sorted(map(list,co.items())) and c['provenance_coefficients']==sorted(map(list,pr.items()))
        assert (F(c['load']),F(c['initial_simplex_mass']),F(c['slack']))==(load,mass,cap*mass-load)
        newcap=cap
        if cap==2:
            core=c['core'];aut,H,scales,J,det,pdroots=cache[core]
            ports=[[i for i,e in enumerate(o[0]) if e[0]==3] for o in orbit]
            overlap=sorted(set(orbit[0][1])&set(orbit[1][1]))
            if scales is None:
                cross=sum((J[i][j] for i in ports[0] for j in ports[1]),Phi());csq=cross*cross;product=cross.out()
            else:
                # Derive square from the independently built positive gauge,
                # not the producer's closed formula for a two-node core.
                assert len(ports[0])==len(ports[1])==1
                i,j=ports[0][0],ports[1][0]
                csq=scales[i]*scales[j]*J[i][j]*J[i][j];product=None
            bad=not csq.integral();take=bool(overlap) and bad
            if take:newcap=1;tightened.append(c['index'])
            z=pairchecks[c['index']]
            assert z['ports']==ports and z['common_base_vertices']==overlap
            assert z['cross_product']==product and z['cross_two_cycle']==csq.out()
            assert [z['cross_two_cycle_trace'],z['cross_two_cycle_norm']]==csq.trace_norm()
            assert z['projection_norms']==[c['repeated_wall_product']]*2
            assert z['nonintegral_two_cycle']==bad and z['tightened']==take
            crosschecks.append(dict(column=c['index'],common_base_vertices=overlap,cycle=csq.out(),integral=not bad,capacity=newcap))
        for key in ('source_bases','positive_face_coefficients','provenance_coefficients','orientation_orbit','simplex_face','load','initial_simplex_mass'):
            assert xc[key]==c[key]
        assert xc['capacity']==newcap and xc['inherited_capacity']==cap and F(xc['slack'])==newcap*mass-load
        if load>cap*mass:viol.append(c['index'])
        if load>newcap*mass:xviol.append(c['index'])
        required[si]=max(required[si],load/cap);xrequired[si]=max(xrequired[si],load/newcap)
        loads.append((c['index'],load,cap,newcap))
    assert len(groups)==len(expected)==815
    assert sorted(viol)==P['violating_columns'] and len(viol)==19
    assert sorted(xviol)==X['violating_columns'] and len(xviol)==42
    assert sorted(tightened)==X['strengthened_columns'] and len(tightened)==215
    assert len(crosschecks)==243 and sum(v>0 for _,v,_,_ in loads)==167
    assert sorted(set(xviol)-set(viol))==X['new_violating_columns']
    for output,req in [(P,required),(X,xrequired)]:
        assert [[i,str(alpha.get(i,F())),str(v)] for i,v in sorted(req.items())]==output['required_simplex_masses']
    for c in O['columns']:
        n=expected[tuple(c['source_bases'])]
        for key in ('capacity','simplex_face','positive_face_coefficients'):assert n[key]==c[key]
    # Independent small degree census; the prime-factor bound is in the report.
    from math import gcd
    orders=[n for n in range(1,121) if sum(gcd(i,n)==1 for i in range(1,n+1))<=4]
    assert orders==P['cyclotomic_orders_degree_at_most_two']==[1,2,3,4,5,6,8,10,12]
    # Distinct-profile successor. Reconstruct ALL selected-domain orientation
    # products, then maximal cliques by graph recursion rather than the
    # producer's all-subsets scan.
    CL=D['cliques'];bytype=defaultdict(list)
    for c in X['columns']:
        if c['capacity']==1:bytype[c['base_type']].append(c)
    savedpairs={tuple(z['columns']):z for z in CL['pair_checks']}
    edgepairs=set();cpchecks=[];tested=0
    for bt,cs in sorted(bytype.items()):
        for a,b in combinations(cs,2):
            tested+=1;key=(a['index'],b['index']);record=savedpairs[key]
            assert a['core']==b['core'] and a['simplex_face']==b['simplex_face']
            assert not(set(a['source_bases'])&set(b['source_bases']))
            aut,H,scales,J,det,pdroots=cache[a['core']]
            if scales is not None:
                assert record['reason']=='outside_selected_physical_Qsqrt5_core_domain'
                cpchecks.append(dict(columns=list(key),excluded=False,reason='outside_selected_domain'));continue
            allgood=True;crosses=[]
            for ia,oa in enumerate(a['orientation_orbit']):
                for ib,ob in enumerate(b['orientation_orbit']):
                    overlap=sorted(set(oa[1])&set(ob[1]))
                    assert all(e[0] in (2,3) for e in oa[0]+ob[0])
                    pa=[i for i,e in enumerate(oa[0]) if e[0]==3];pb=[i for i,e in enumerate(ob[0]) if e[0]==3]
                    assert sum((J[i][j] for i in pa for j in pa),Phi()).out()==a['repeated_wall_product']
                    assert sum((J[i][j] for i in pb for j in pb),Phi()).out()==b['repeated_wall_product']
                    cross=sum((J[i][j] for i in pa for j in pb),Phi());sq=cross*cross
                    good=bool(overlap) and not sq.integral();allgood &= good
                    crosses.append(dict(orientations=[ia,ib],common_base_vertices=overlap,
                        cross_product=cross.out(),cross_two_cycle=sq.out(),trace=sq.trace_norm()[0],norm=sq.trace_norm()[1]))
            if allgood:
                edgepairs.add(key)
                assert crosses==record['checks']
            else:
                assert record['reason']=='integral_two_cycle_not_excluded'
            cpchecks.append(dict(columns=list(key),excluded=allgood,checks=crosses))
    assert tested==len(savedpairs)==80 and len(edgepairs)==70
    assert edgepairs=={tuple(z['columns']) for z in CL['incompatibility_edges']}
    maximal=[]
    for bt,cs in sorted(bytype.items()):
        verts={c['index'] for c in cs}
        adj={v:{w for w in verts if tuple(sorted((v,w))) in edgepairs} for v in verts}
        found=[]
        def bk(r,p,x):
            if not p and not x:
                if len(r)>1:found.append(tuple(sorted(r)))
                return
            for v in sorted(p.copy()):
                bk(r|{v},p&adj[v],x&adj[v]);p.remove(v);x.add(v)
        bk(set(),verts.copy(),set())
        maximal.extend((bt,c) for c in sorted(found))
    expectedcliques={(c['base_type'],tuple(c['profile_columns'])):c for c in CL['columns'][815:]}
    assert set(maximal)==set(expectedcliques) and len(maximal)==36
    cr=defaultdict(F);cv=[];cloads=[]
    for a,c in zip(X['columns'],CL['columns'][:815]):assert a==c
    for (bt,vs),c in expectedcliques.items():
        co=Counter();pr=Counter();bases=[]
        for v in vs:
            z=xcols[v];co.update(dict(z['positive_face_coefficients']));pr.update(dict(z['provenance_coefficients']));bases+=z['source_bases']
            assert z['base_type']==bt and z['simplex_face']==c['simplex_face'] and z['capacity']==1
        assert len(bases)==len(set(bases)) and sorted(bases)==c['source_bases']
        assert sorted(map(list,co.items()))==c['positive_face_coefficients']
        assert sorted(map(list,pr.items()))==c['provenance_coefficients'] and c['capacity']==1
    for c in CL['columns']:
        si=c['simplex_face'];load=sum((alpha.get(i,F())*m for i,m in c['positive_face_coefficients']),F());mass=alpha.get(si,F())
        assert (F(c['load']),F(c['initial_simplex_mass']),F(c['slack']))==(load,mass,c['capacity']*mass-load)
        cr[si]=max(cr[si],load/c['capacity'])
        if load>c['capacity']*mass:cv.append(c['index'])
        cloads.append((c['index'],load))
    assert len(cv)==46 and cv==CL['violating_columns']
    assert len([v for v in cv if v>=815])==4
    assert [[i,str(alpha.get(i,F())),str(v)] for i,v in sorted(cr.items())]==CL['required_simplex_masses']
    paths=[f'scripts/{TAG}.py']+[f'results/{t}.json' for t in TAGS.values()]
    paths += [f'scripts/{TAGS[k]}.py' for k in ('orbit','cross','cliques')]
    return dict(status='PASS_INDEPENDENT_ORBIT_AND_CROSS_ORIENTATION_CAP_CAPACITIES',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        selected_occurrences=815,core_types=len(cache),orbit_columns=815,inherited_columns=241,
        orbit_violations=19,cross_pair_checks=243,tightened_capacities=215,cross_violations=42,
        candidate_counts=dict(candidate_counts),positive_loads=167,checks=checks,cross_checks=sorted(crosschecks,key=lambda z:z['column']),
        profile_pair_checks=cpchecks,clique_columns=len(maximal),final_columns=851,final_violations=46,
        maximal_cliques=[dict(base_type=bt,profile_columns=list(vs)) for bt,vs in maximal],
        core_checks=[dict(core=i,nodes=len(cache[i][1]),automorphism_count=len(cache[i][0]),
          determinant=cache[i][4].out(),complete_PD_deletions=cache[i][5],positive_gauge=cache[i][2] is not None) for i in sorted(cache)],
        scope='Independent phi-basis Fraction arithmetic, complete selected-core symmetry/PD checks, raw incidence maps,815cap selection and loads,243cross-cycle checks,80profile pairs and36maximal clique columns. Inherited arithmetic core realization and maximal-simplex/prism shared-base theorem assumed; no optimizer or dimension claim.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    r=clean(run());out=ROOT/f'results/{TAG}.json'
    if args.verify:assert r==json.loads(out.read_text())
    else:out.write_text(json.dumps(r,indent=2)+'\n')
    print(r['status'],r['selected_occurrences'],r['core_types'],r['tightened_capacities'],flush=True)
