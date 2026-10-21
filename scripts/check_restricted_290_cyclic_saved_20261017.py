#!/usr/bin/env python3
"""No-project-import exact cycle/cover accounting check; not physical replay."""
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_restricted_290_cyclic_cover_20261017'

def q(x):return tuple(map(F,x))
def mul(x,y):return(x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def sign(x):
    a,b=x
    if not b:return(a>0)-(a<0)
    if not a:return(b>0)-(b<0)
    if a*b>0:return(a>0)-(a<0)
    n=a*a-5*b*b
    return((a>0)-(a<0))*((n>0)-(n<0))
def integral(x):return(2*x[1]).denominator==1 and(x[0]-x[1]).denominator==1
def native_square(e):
    if e[0]==0:return q(e[1:])
    if e[0]==2:return(F(),F())
    assert e[0]in(3,5)
    return(F(1),F())if e[0]==3 else(F(3,2),F(1,2))

def run():
    data=json.loads((ROOT/f'results/{TAG}.json').read_text())
    base=json.loads((ROOT/f'results/{data["baseline_source"]}.json').read_text())
    old={r['vertex']:r for r in base['forest_records']};alpha={i:F(x)for i,x in base['face_multipliers']}
    mu={i:F(x)for i,x in base['baseline_restricted_repair']['vertex_multipliers']};counts=Counter()
    for path,h in data['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h;counts['bindings']+=1
    for r in data['forest_records']:
        b=old[r['vertex']];d=list(map(F,r['demands']));groups=b['corner_orbits']
        assert r['forest']==b['forest'] and r['pair_domain_sha256']==b['all_unowned_deletion']['pair_domain_sha256']
        assert d==[sum((alpha.get(i,F())*F(x)for i,x in g['positive_face_coefficients']),F())for g in groups]
        owner={i:j for j,g in enumerate(groups)for i in g['literal_stars']}
        wallcontacts=r['native_wall_contacts'];starwalls={c['star']:c['walls']for c in r['native_star_checks']}
        pairs={tuple(c['walls']):c for c in r['cross_wall_checks']}
        for c in pairs.values():
            s=q(c['squared_product']);assert sign(s)>=0 and integral(s)
            assert c['product_sign']==(-1 if s!=(0,0)else 0)
        for c in r['native_star_checks']:
            key=tuple(sorted(c['walls']));assert len(set(key))==2
            assert q(pairs[key]['squared_product'])==native_square(c['native_omitted_product'])
        qrow=b['all_unowned_deletion']
        expected=Counter()
        source=[(qrow['literal_representatives'][i],x)for i,x in qrow['certificate']['mixture']]if 'certificate'in qrow else qrow['split_cover']['terms']
        for cl,x in source:
            cl=tuple(i for i in cl if d[owner[i]]>0);assert cl;expected[cl]+=F(x)
        assert expected==Counter({tuple(t['literal_stars']):F(t['amount'])for t in r['cover_terms']})
        marginal=[F()for _ in d];mass=F();baseedges=None;usedpairs=set()
        for t in r['cover_terms']:
            amount=F(t['amount']);assert amount>0;mass+=amount
            for i in t['literal_stars']:marginal[owner[i]]+=amount
            ids=t['walls'];assert ids==sorted({w for i in t['literal_stars']for w in starwalls[i]})
            local={w:13+j for j,w in enumerate(ids)};edges={(i,j):q(s)for i,j,s in t['nonzero_squared_edges']}
            assert len(edges)==len(t['nonzero_squared_edges'])
            be={k:s for k,s in edges.items()if k[1]<13}
            if baseedges is None:baseedges=be
            else:assert be==baseedges
            expectededges=dict(baseedges)
            for w in ids:
                for i,e in wallcontacts[w]:expectededges[i,local[w]]=native_square(e)
            for a,bw in combinations(ids,2):
                usedpairs.add((a,bw));s=q(pairs[a,bw]['squared_product'])
                if s!=(0,0):expectededges[local[a],local[bw]]=s
            assert edges==expectededges
            pots=list(map(q,t['tree_squareclass_potentials']));parents=t['spanning_parents'];n=13+len(ids)
            assert len(pots)==len(parents)==n and all(sign(x)>0 for x in pots)
            assert t['spanning_roots']==[i for i,x in enumerate(parents)if x is None]
            for i,pa in enumerate(parents):
                if pa is None:assert pots[i]==(1,0)
                else:assert pots[i]==mul(pots[pa],edges[tuple(sorted((i,pa)))])
                seen=set();j=i
                while j is not None:assert j not in seen;seen.add(j);j=parents[j]
            roots={(i,j):q(z)for i,j,z in t['edge_membership_roots']};assert roots.keys()==edges.keys()
            for key,s in edges.items():
                i,j=key;z=roots[key];assert i<j and sign(s)>0 and integral(s)
                assert sign(z)>0 and integral(z) and mul(z,z)==mul(mul(pots[i],pots[j]),s)
            rank=len(edges)-n+len(t['spanning_roots']);assert rank==t['cycle_rank']>=0
            counts.update(cover_terms=1,edge_membership_checks=len(edges),cycle_rank=rank,
                external_wall_occurrences=len(ids),singleton_star_terms=int(len(t['literal_stars'])==1),
                multi_star_terms=int(len(t['literal_stars'])>1))
        assert usedpairs==pairs.keys()
        assert marginal==list(map(F,r['covered_marginals'])) and all(x>=y for x,y in zip(marginal,d))
        assert mass==F(r['cover_mass'])==F(qrow['optimal_mass_upper_bound'])
        m=max(mu.get(r['vertex'],F()),mass);assert m==F(r['preserved_multiplier']);mu[r['vertex']]=m
        counts.update(forests=1,distinct_external_walls=len(wallcontacts),
            distinct_cross_wall_pairs=len(pairs),native_star_checks=len(starwalls))
    assert all(counts[k]==v for k,v in data['counts'].items())
    h=data['hypothetical_restricted_witness'];assert sum(mu.values(),F())==F(h['normalizer'])
    assert F(h['objective_numerator'])/F(h['normalizer'])==F(h['bound'])>12
    assert [[i,str(x)]for i,x in sorted(mu.items())if x]==h['vertex_multipliers']
    print('PASS_NO_PROJECT_IMPORT_CYCLIC_SAVED_CHECK',dict(counts),flush=True)

if __name__=='__main__':run()
