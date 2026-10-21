#!/usr/bin/env python3
"""Independent literal mapping and Fraction projection checks for selected caps.

No producer or project arithmetic imports. This reviews the selected valid
capacity family; omitted cap types are harmless and are not re-enumerated.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_cross_coface_cap_capacity_20261017_review'
PRIMARY = 'fresh_sqrt5_cross_coface_cap_capacity_20261017'
TAGS = dict(primary=PRIMARY,
    slices='fresh_sqrt5_labelled_coface_base_packing_20261017',
    incidence='fresh_sqrt5_prism_simplex_base_incidence_20261017',
    actual='fresh_sqrt5_actual_provenance_dual_20261017',
    quad='fresh_sqrt5_rank13_packing_quad_separation_20261017',
    raw='fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017')


class K:
    def __init__(self, a=0, b=0):
        if isinstance(a, K): self.a, self.b = a.a, a.b
        else: self.a, self.b = F(a), F(b)
    def __add__(self, x):
        x=K(x); return K(self.a+x.a,self.b+x.b)
    __radd__=__add__
    def __neg__(self): return K(-self.a,-self.b)
    def __sub__(self,x): return self+-K(x)
    def __rsub__(self,x): return K(x)+-self
    def __mul__(self,x):
        x=K(x); return K(self.a*x.a+5*self.b*x.b,self.a*x.b+self.b*x.a)
    __rmul__=__mul__
    def __truediv__(self,x):
        x=K(x); d=x.a*x.a-5*x.b*x.b
        assert d; return self*K(x.a/d,-x.b/d)
    def __bool__(self): return bool(self.a or self.b)
    def __eq__(self,x):
        x=K(x); return (self.a,self.b)==(x.a,x.b)
    def sign(self):
        a,b=self.a,self.b
        if not b:return (a>0)-(a<0)
        if not a:return (b>0)-(b<0)
        if a*b>0:return (a>0)-(a<0)
        delta=a*a-5*b*b
        return ((a>0)-(a<0))*((delta>0)-(delta<0))
    def out(self):return [str(self.a),str(self.b)]


def inverse(A):
    n=len(A); M=[[K(x) for x in r]+[K(i==j) for j in range(n)] for i,r in enumerate(A)]
    det=K(1)
    for j in range(n):
        p=next(i for i in range(j,n) if M[i][j])
        if p!=j:M[p],M[j]=M[j],M[p];det=-det
        z=M[j][j];det=det*z;M[j]=[x/z for x in M[j]]
        for i in range(n):
            if i!=j:
                z=M[i][j]
                if z:M[i]=[a-z*b for a,b in zip(M[i],M[j])]
    assert all(M[i][j]==int(i==j) for i in range(n) for j in range(n))
    return [r[n:] for r in M],det


def pd(A):
    M=[[K(x) for x in r] for r in A]
    for k in range(len(M)):
        d=M[k][k]
        if d.sign()<=0:return False
        for i in range(k+1,len(M)):
            for j in range(k+1,len(M)):
                M[i][j]=M[i][j]-M[i][k]*M[k][j]/d
    return True


def edge(m,a='0',b='0'):
    return (int(m),str(F(a)),str(F(b))) if m==0 else (int(m),'0','0')


RIGHT=edge(2)


def core_graph(c):
    n=c['nodes'];G=[[RIGHT for _ in range(n)] for _ in range(n)]
    for i,j,(m,a,b) in c['edges']:G[i][j]=G[j][i]=edge(m,a,b)
    return G


def raw_graph(r):
    n=r['nodes'];G=[[RIGHT for _ in range(n)] for _ in range(n)]
    for (i,j),m in zip(combinations(range(n),2),r['row']['labels']):
        e=edge(m,*r['parameters'][f'{i},{j}']) if m==0 else edge(m)
        G[i][j]=G[j][i]=e
    return G


def automorphisms(G):
    """Exhaustive adjacency-preserving backtracking, without producer colours."""
    n=len(G); sig=[tuple(sorted(G[i][j] for j in range(n) if j!=i)) for i in range(n)]
    domains={i:[j for j in range(n) if sig[i]==sig[j]] for i in range(n)}
    order=sorted(range(n),key=lambda i:(len(domains[i]),i));out=[]
    def visit(k,mp,used):
        if k==n:out.append([mp[i] for i in range(n)]);return
        i=order[k]
        for j in domains[i]:
            if j in used or any(G[i][a]!=G[j][b] for a,b in mp.items()):continue
            mp[i]=j;visit(k+1,mp,used|{j});del mp[i]
    visit(0,{},set())
    return out


def squared(e):
    m,a,b=e
    if m==0:return K(a,b)
    return {2:K(),3:K(1),4:K(2),5:K(F(3,2),F(1,2)),
            6:K(3),10:K(F(5,2),F(1,2))}[m]


def physical_core(G):
    """Literal norm-2 K matrix, or independent tree gauge over K.

    The latter has Hii=2*s_i, Hij=-s_i*w_ij along a parent-child
    edge, with s_child=s_parent*w_ij. It is positively congruent to
    the physical normalized matrix. One-contact projection norms then
    equal s_i*(H^-1)ii and require no adjoined square roots.
    """
    n=len(G)
    if all(G[i][j][0] in (2,3,5) for i in range(n) for j in range(i)):
        vals={2:K(),3:K(-1),5:K(F(-1,2),F(-1,2))}
        return [[K(2) if i==j else vals[G[i][j][0]] for j in range(n)] for i in range(n)],None
    E=[(i,j) for i in range(n) for j in range(i) if G[i][j]!=RIGHT]
    assert len(E)==n-1
    H=[[K() for _ in range(n)] for _ in range(n)];s={0:K(1)};todo=[0]
    for i in todo:
        for j in range(n):
            if i==j or G[i][j]==RIGHT or j in s:continue
            w=squared(G[i][j]);assert w.sign()>0
            s[j]=s[i]*w;H[i][j]=H[j][i]=-s[j];todo.append(j)
    assert len(s)==n and all(v.sign()>0 for v in s.values())
    for i in range(n):H[i][i]=2*s[i]
    return H,[s[i] for i in range(n)]


def rank(name):
    if name.startswith('I2('):return 2
    assert name[0] in 'ABCDEH' and name[1:].isdigit(),name
    return int(name[1:])


def run():
    D={k:json.loads((ROOT/f'results/{v}.json').read_text()) for k,v in TAGS.items()}
    P,S,I,A,Q,R=[D[k] for k in ('primary','slices','incidence','actual','quad','raw')]
    assert P['status']=='PASS_EXACT_CROSS_COFACE_A1_CAP_CAPACITIES'
    profiles={}
    for p in S['coface_slice_demands']:
        key=(p['simplex_face'],p['simplex_dimension'],p['core'])
        if p['base_type'] in profiles:assert profiles[p['base_type']]==key
        profiles[p['base_type']]=key
    alpha={i:F(a) for i,a in S['face_multipliers']};assert all(a>=0 for a in alpha.values())
    templates={t['canonical_core']:t for t in I['prism_templates']}
    savedcores={c['core']:c for c in S['core_symmetry_and_subset_orbits']}
    caches={};checks=[];groups=defaultdict(list);seen=set()
    for o in P['occurrences']:
        B=S['actual_bases'][o['base_index']]
        assert B['index']==o['base_index']
        j,b,q=B['provenance'],B['deleted'],B['coface_root']
        assert (j,b) not in seen;seen.add((j,b))
        assert B['coface_location']==['finite','A1',0]
        prov=A['provenance_masses'][j];ci=prov['canonical_prism_core']
        C=Q['canonical_cores'][ci];r=R['records'][C['source_core']];G=raw_graph(r)
        t=B['base_type'];typ=I['base_types'][t];si,dim,core=profiles[t]
        assert core==typ['core'] and dim==typ['simplex_dimension']
        T=next(x for x in templates[ci]['bases'] if x['deleted_normal']==b)
        assert T['opposite_base_facet']==q and T['base_core']==core
        mp=dict(T['core_isomorphism']);cg=core_graph(I['cores'][core])
        assert set(mp)==set(range(r['nodes']))-{b,q}
        assert sorted(mp.values())==list(range(len(cg)))
        assert T['finite_extra']==['A1']
        assert G[b][q][0]==0 and all(G[q][i]==RIGHT for i in mp)
        assert all(G[i][j]==cg[mp[i]][mp[j]] for i,j in combinations(mp,2))
        # All corners of this particular base, reconstructed from the
        # canonical K2,2 union, are labelled by the OTHER omitted root.
        vertices=set();base_pairs=set()
        for left,right in C['quad_partitions']:
            if b in left:side,other=left,right
            elif b in right:side,other=right,left
            else:continue
            assert set(side)=={b,q}
            vertices.update(other)
            base_pairs.update(tuple(sorted((b,p))) for p in other)
        assert vertices and vertices<=set(mp)
        assert base_pairs<={tuple(sorted(z)) for z in C['pd_deletions']}
        canonical=sorted(mp[p] for p in vertices)
        assert canonical==B['canonical_vertices']==o['base_vertices']
        assert len(canonical)==B['prism_dimension']==templates[ci]['prism_dimension']
        assert set(canonical)<=set(I['cores'][core]['spherical_deletions'])
        assert sorted(mp[p] for p in T['full_spherical_deletions'])==I['cores'][core]['spherical_deletions']
        copies=Counter(typ['spectator_forest'])['A1']
        assert copies==o['A1_copies'] and copies>=1
        assert len(cg)+sum(map(rank,typ['spectator_forest']))==14
        assert r['nodes']+sum(map(rank,typ['spectator_forest']))-1==15
        if core not in caches:
            aut=automorphisms(cg)
            assert aut==[list(range(len(cg)))]==savedcores[core]['automorphisms']
            H,scales=physical_core(cg);J,det=inverse(H)
            assert det.sign()<0
            pd_roots=[]
            for k in range(len(cg)):
                ids=[i for i in range(len(cg)) if i!=k]
                if pd([[H[i][j] for j in ids] for i in ids]):pd_roots.append(k)
            assert pd_roots==I['cores'][core]['spherical_deletions']
            assert pd_roots
            caches[core]=(J,scales,det,pd_roots)
        J,scales,det,pd_roots=caches[core]
        contact=[]
        for i in sorted(mp,key=mp.get):
            e=G[b][i];contact.append([mp[i],e[0],list(e[1:]) if e[0]==0 else None])
        assert contact==o['core_contacts']
        # All selected cap-to-core contacts happen to be right or pi/3;
        # this is rechecked rather than assumed by the geometric theorem.
        assert all(m in (2,3) for _,m,_ in contact)
        ports=[i for i,m,_ in contact if m==3]
        if scales is None:norm=sum((J[i][j] for i in ports for j in ports),K())
        else:
            assert len(ports)==1
            k=ports[0];norm=scales[k]*J[k][k]
        tau=K(*G[b][q][1:]);assert (tau-4).sign()>0 and (6-tau).sign()>0
        assert norm==2-tau/2 and norm.sign()<0 and (norm+1).sign()>0
        assert tau.out()==o['dotted_squared_norm2_product'] and norm.out()==o['common_projection_norm']
        assert (j,B['face'],ci,C['source_core'],b,q,core)==tuple(o[k] for k in
            ('provenance','face','canonical_prism_core','source_core','deleted_normal','original_coface_root','core'))
        assert sorted(map(list,mp.items()))==o['core_isomorphism']
        key=(t,json.dumps(contact,separators=(',',':')),tuple(canonical),tuple(tau.out()))
        groups[key].append((B,si))
        checks.append(dict(base_index=B['index'],core=core,canonical_vertices=canonical,
            raw_vertex_omitted_pairs=sorted(map(list,base_pairs)),projection_norm=norm.out(),
            dotted_squared_product=tau.out(),A1_copies=copies))
    expected={tuple(c['source_bases']):c for c in P['columns']}
    loads=[];violations=[];required=defaultdict(F)
    for bs in groups.values():
        ids=tuple(sorted(b['index'] for b,_ in bs));c=expected[ids]
        assert len({si for _,si in bs})==1
        si=bs[0][1];f=Counter(b['face'] for b,_ in bs);pr=Counter(b['provenance'] for b,_ in bs)
        load=sum((m*alpha.get(i,F()) for i,m in f.items()),F());old=alpha.get(si,F())
        assert c['simplex_face']==si and c['capacity']==1
        assert c['positive_face_coefficients']==sorted(map(list,f.items()))
        assert c['provenance_coefficients']==sorted(map(list,pr.items()))
        assert (F(c['load']),F(c['initial_simplex_mass']),F(c['slack']))==(load,old,old-load)
        if load>old:violations.append(c['index'])
        required[si]=max(required[si],load)
        loads.append((c['index'],load,old))
    assert len(checks)==len(groups)==len(expected)==241
    assert sorted(violations)==P['violating_columns'] and len(violations)==53
    assert sum(load>0 for _,load,_ in loads)==70
    assert [[i,str(alpha.get(i,F())),str(v)] for i,v in sorted(required.items())]==P['necessary_simplex_mass_lower_bounds']
    paths=[f'scripts/{TAG}.py',f'scripts/{PRIMARY}.py']+[f'results/{tag}.json' for tag in TAGS.values()]
    return dict(status='PASS_INDEPENDENT_SELECTED_CROSS_COFACE_CAP_CAPACITIES',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        selected_occurrences=len(checks),capacity_columns=len(groups),core_types=len(caches),
        positive_loads=70,violations=53,checks=checks,
        core_checks=[dict(core=i,determinant_of_positive_gauge=caches[i][2].out(),
            complete_PD_deletions=caches[i][3],complete_automorphism_group_size=1)
            for i in sorted(caches)],
        scope='Independent stdlib Fraction geometry/mapping/load review of every selected241capacity. Direct13core projection and completePDdeletion/automorphism checks. Inherited raw arithmetic/local-realization and whole-star incidence theorems are not recertified. No census of omitted cap types, no optimizer, no repaired witness or dimension bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['selected_occurrences'],d['core_types'],d['violations'],flush=True)
