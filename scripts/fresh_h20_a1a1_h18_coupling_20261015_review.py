#!/usr/bin/env python3
"""Independent physical168-bank reconstruction and all-label proof checks."""
import argparse,json
from itertools import combinations,product
from pathlib import Path
import sympy as s
import fresh_h20_a2_h18_coupling_20261014_review as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_a1a1_h18_coupling_20261015_review'
PRIMARY='results/fresh_h20_a1a1_h18_coupling_20261015.json'
SOURCE='data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json'
B2='results/fresh_h20_b2_h18_coupling_20261014.json'
K,Z,O,rt2,rt3,rt5,phi=A.K,A.Z,A.O,A.rt2,A.rt3,A.rt5,A.phi

def run():
    d,raw,b2=map(A.read,(PRIMARY,SOURCE,B2))
    G=[[A.dec(v)for v in row]for row in raw['gram']]
    basis=list(range(19));assert basis==d['source_basis']
    Q0=[[G[i][j]for j in basis]for i in basis];J=A.inverse(Q0)
    roots=[A.mv(J,[G[j][i]for j in basis])for i in range(40)]
    norms=[G[i][i]for i in range(40)]
    assert all(A.pair(roots[i],roots[j],Q0)==G[i][j]for i in range(40)for j in range(40))
    assert roots==[[A.dec(v)for v in row]for row in d['source_normals']]
    edges=[(i,j)for i,j in combinations(range(40),2)if A.sign(G[i][j])<0 and 4*G[i][j]**2==norms[i]*norms[j]]
    component={0}
    while True:
        nxt=component|{j for i,j in edges if i in component}|{i for i,j in edges if j in component}
        if nxt==component:break
        component=nxt
    assert len(edges)==27 and len(component)==24 and sorted(component)==d['forced_ordinary']
    assert set(map(frozenset,edges))==set(map(frozenset,d['pi3_edges']))
    assert all(norms[i]==2*O for i in component)
    assert all(not G[i][j]for i,j in combinations([0,2,4],2))
    subsets=[{0},{1},{0,1}]
    impossible_count=sum(all(not(a&b)for a,b in combinations(triple,2))for triple in product(subsets,repeat=3))
    assert impossible_count==0
    eligible=sorted(set(range(40))-component);core=[17,23,24]
    assert eligible==d['eligible_special']and core==d['core']
    assert sum(norms[i]==O for i in eligible if i not in core)==11
    assert sum(norms[i]==2*O for i in eligible if i not in core)==2
    # All nonzero contact magnitudes are at least1. Two positive contacts
    # have h>0 only at33,34,43,35,53; labels>=6 leave h<=0.
    angle={0:Z,3:O,4:rt2,5:phi,6:rt3}
    positive_doubles=[]
    for m,n in product((3,4,5,6),repeat=2):
        h=2*O-(angle[m]**2+angle[n]**2)/2
        if A.sign(h)>0:positive_doubles.append([m,n])
    assert positive_doubles==[[3,3],[3,4],[3,5],[4,3],[5,3]]
    assert sorted(positive_doubles)==sorted(d['initial_positive_double_profiles'])
    # Independent exact gap endpoints. At m5 the products lie strictly
    # between consecutive Coxeter values; monotonic h handles every tail.
    h5=(5*O-rt5)/4;h8=O-rt2/2
    assert A.sign(4*h5-phi**2)>0 and A.sign(3*O-4*h5)>0
    assert A.sign(2*h5-O)>0 and A.sign(2*O-2*h5)>0
    assert A.sign(h5)>0 and A.sign(O-h5)>0
    h35=(3*O-rt5)/4
    assert A.sign(h35)>0 and A.sign(O-4*h35)>0
    assert A.sign(h8-O/4)>0 and A.sign(O/2-h8)>0
    # For u^2=7/8, f(u)=4u^3-3u-1/2=(u-1)/2<0.
    # f increases on (1/2,1), so cos(pi/9)>u and h9<1/4.
    u2=s.Rational(7,8)
    assert 4*u2-3==s.Rational(1,2)and s.Rational(3,4)<u2<1 and 12*u2-3>0
    profile_labels=[(0,0),(3,0),(0,3),(4,0),(0,4),(6,0),(0,6),(3,3),(3,4),(4,3)]
    factors={2*O:rt2,3*O/2:rt2*rt3/2,O:O,O/2:rt2/2}
    profiles={}
    for labels in profile_labels:
        aa=[angle[m]for m in labels];h=2*O-sum((t*t for t in aa),Z)/2
        f=factors[h];assert f*f==h and A.sign(f)>0
        name='O'if labels==(0,0)else ''.join(map(str,labels))
        profiles[name]=(labels,aa,h,f)
        assert d['allowed_profiles'][name]['labels']==list(labels)
        assert K.convert(s.Rational(d['allowed_profiles'][name]['h']))==h
        assert A.dec(d['allowed_profiles'][name]['contact_sum_square'])==sum(aa,Z)**2
    core_names=['O','40','04','33'];assignments=[];anchors=[]
    for i in range(40):
        names=['O']if i in component else core_names if i in core else list(profiles)
        assignments.append(dict(source_wall=i,profiles=names))
        if i in component:continue
        saved=next(row for row in d['ordinary_anchors']if row['source_wall']==i)
        j=saved['ordinary_anchor'];g=-O if i in core else-2*O
        assert j in component and G[i][j]==g and A.dec(saved['source_pairing'])==g
        assert A.dec(saved['source_norm'])==norms[i]
        for name in names:
            h=profiles[name][2];prod2=g*g*h/norms[i]
            assert prod2 in (O,2*O,3*O)or A.sign(prod2-4*O)>=0
            relation={O:'pi3',2*O:'pi4',3*O:'pi6',4*O:'parallel'}.get(prod2,'dotted')
            anchors.append(dict(source_wall=i,profile=name,squared_product=A.enc(prod2),relation=relation))
            row=next(row for row in d['complete_anchor_relations']if row['source_wall']==i and row['profile']==name)
            assert A.dec(row['squared_product'])==prod2 and row['relation']==relation
    assert assignments==d['per_source_profiles']and len(anchors)==142
    Q=[row+[Z,Z]for row in Q0]+[[Z]*19+[2*O,Z],[Z]*19+[Z,2*O]]
    bank=[[Z]*19+[O,Z],[Z]*19+[Z,O]];records=[]
    for row in assignments:
        i=row['source_wall'];rootnorm=O if norms[i]==O else rt2
        for name in row['profiles']:
            labels,aa,h,f=profiles[name]
            normal=[f*u/rootnorm for u in roots[i]]+[-aa[0]/2,-aa[1]/2]
            assert A.pair(normal,normal,Q)==2*O
            assert A.pair(normal,bank[0],Q)==-aa[0]and A.pair(normal,bank[1],Q)==-aa[1]
            bank.append(normal);records.append(dict(source=i,profile=name))
    assert len(bank)==168==d['relaxed_bank_count']
    x=[A.dec(v)for v in b2['source_point']];h=[sum((r[j]for r in roots),Z)for j in range(19)]
    xn=A.pair(x,x,Q0);hn=A.pair(h,h,Q0);xh=A.pair(x,h,Q0)
    assert xn==(3857*O-2800*rt2)/529 and hn==-2952*O-1280*rt2 and xh==-16*O-40*rt2
    gaps=[A.pair(x,r,Q0)for r in roots]
    assert all(A.sign(g)<=0 for g in gaps)and all(A.sign(A.pair(h,r,Q0))<0 for r in roots)
    assert all(A.sign(gaps[i]+(O if norms[i]==O else rt2))<=0 for i in eligible)
    eps=O/100000
    E=[2*u+eps*v for u,v in zip(x,h)]+[-O/2,-O/2]
    H=[100*u for u in h]+[-O/2,-O/2]
    assert E==[A.dec(v)for v in d['strict_escape']]and H==[A.dec(v)for v in d['strict_center']]
    margins={name:4*p[2]-sum(p[1],Z)**2/4 for name,p in profiles.items()}
    assert all(A.sign(v)>0 for v in margins.values())and margins['34']==(5*O-2*rt2)/4
    en=A.pair(E,E,Q);cn=A.pair(H,H,Q);ec=A.pair(E,H,Q)
    assert en==A.dec(d['escape_norm'])and cn==A.dec(d['center_norm'])and ec==A.dec(d['escape_center_pairing'])
    assert A.sign(en)>0 and A.sign(cn)<0 and A.sign(ec)<0
    ep=[A.pair(E,n,Q)for n in bank];cp=[A.pair(H,n,Q)for n in bank]
    assert all(A.sign(v)<0 for v in ep+cp)
    paths=[PRIMARY,SOURCE,B2,f'scripts/{TAG}.py','scripts/fresh_h20_a2_h18_coupling_20261014_review.py']
    return dict(status='PASS_INDEPENDENT_ALL_DIRECT_A1A1_H18_ESCAPE',sources={p:A.sha(p)for p in paths},
        arithmetic='Separate SymPy QQ(sqrt2,sqrt3,sqrt5) reconstruction using the independent A2 reviewer; no A1A1 producer or Fraction engine imports.',
        source_gram_entries=1600,pi3_edges=[list(e)for e in edges],forced_ordinary=sorted(component),
        orthogonal_triple=[0,2,4],disjoint_three_nonempty_color_patterns=impossible_count,
        eligible_special=eligible,anchor_relations=anchors,bank_count=168,bank_profiles=records,
        initial_doubles=positive_doubles,all_label_tail=dict(cos_pi9_comparison_square='7/8',monotonicity_proof_checked=True),
        profile_margins={k:A.enc(v)for k,v in margins.items()},escape_norm=A.enc(en),
        center_norm=A.enc(cn),escape_center_pairing=A.enc(ec),escape_norm_interval=list(map(str,A.interval(en))),
        escape=[A.enc(v)for v in E],center=[A.enc(v)for v in H],
        bank_escape_pairings=[A.enc(v)for v in ep],bank_center_pairings=[A.enc(v)for v in cp],
        scope='Complete necessary all-label inducing profiles and independent physical168-normal bank/336 strict signs. Every coherent direct A1xA1 assignment has a common strict spacelike escape. Additional noninducing walls, changed source metrics and arbitrary larger completions are outside. No full ambient incidence or assignment census claimed.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved independent A1A1 bank replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['bank_count'])
