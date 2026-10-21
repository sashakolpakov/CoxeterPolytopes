#!/usr/bin/env python3
"""All-label necessary contact profiles for a wall crossing the S75 plane."""
import argparse,json
from fractions import Fraction as F
from itertools import combinations
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_s75_projection_gate_20261014'
SOURCE='results/fresh_h20_bridge_faces_20261014_review.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def sign(a,b):
    if not b:return(a>0)-(a<0)
    if not a:return(b>0)-(b<0)
    if a*b>0:return(a>0)-(a<0)
    v=a*a-2*b*b;return((a>0)-(a<0))*((v>0)-(v<0))
def mul(x,y):return(x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return[enc(y)for y in x]
    if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
    return x
def run():
    doc=json.loads((ROOT/SOURCE).read_text())
    for p,h in doc['sources'].items():assert sha(p)==h,p
    r=doc['faces'][1];S=r['support'];G=[list(map(F,row))for row in r['rational_gram']];J=[list(map(F,row))for row in r['inverse']]
    n=len(S);assert n==18 and 75 in S
    assert all(G[i][i]==2 for i in range(n))and min(x for row in J for x in row)>=0
    assert all(sum(G[i][k]*J[k][j]for k in range(n))==int(i==j)for i in range(n)for j in range(n))
    # Independent exact LDL of the saved rational Gram.
    B=[row[:]for row in G];piv=[]
    for k in range(n):
        p=B[k][k];assert p>0;piv.append(p)
        for i in range(k+1,n):
            for j in range(k+1,n):B[i][j]-=B[i][k]*B[k][j]/p
    diag=[J[i][i]for i in range(n)];assert min(diag)==F(1,2)
    assert J[S.index(75)][S.index(75)]==5
    # c_i=cos(pi/m_i) for finite contacts, with c_i>=1 for dotted/parallel.
    # Positive projection requires c^T J c <1/2.
    min3=sum(sorted(diag)[:3])/4;assert min3==F(51,88)>F(1,2)
    forced=[S[i]for i,d in enumerate(diag)if d>=2]
    domains={}
    for i,d in enumerate(diag):
        wall=S[i]
        if d>=2:domains[wall]='right only'
        elif d==F(1,2):domains[wall]='any finite label m>=2; no parallel/dotted'
        elif d==F(10,11):
            assert F(1,2)<1/(2*d)<F(5,8)
            domains[wall]=[2,3,4]
        else:
            assert F(1,4)<1/(2*d)<=F(1,2)
            domains[wall]=[2,3]
    assert domains[17].startswith('any')
    # cos^2(pi/5)=(3+sqrt5)/8>5/8. This rules out all >=5 labels
    # in the paired17 tail, even with the cheapest other pi3 contact.
    paired_tail_bound=F(1,2)*F(5,8)+F(10,11)/4
    assert paired_tail_bound==F(95,176)>F(1,2)
    options=[]
    for wall,dom in domains.items():
        if wall==17:options.extend([(wall,3),(wall,4)])
        elif isinstance(dom,list):options.extend((wall,m)for m in dom if m!=2)
    cos={3:(F(1,2),F(0)),4:(F(0),F(1,2))}
    def cost(profile):
        a=b=F(0)
        for i,m in profile:
            for j,l in profile:
                u,v=mul(cos[m],cos[l]);q=J[S.index(i)][S.index(j)];a+=q*u;b+=q*v
        return a,b
    singles=[];pairs=[];rejected=[]
    for row in options:
        x=cost([row]);assert sign(x[0]-F(1,2),x[1])<0
        if row[0]!=17:singles.append(dict(contacts=[row],cost=x))
    for p in combinations(options,2):
        if p[0][0]==p[1][0]:continue
        x=cost(p);ok=sign(x[0]-F(1,2),x[1])<0
        (pairs if ok else rejected).append(dict(contacts=p,cost=x))
    assert len(singles)==9 and len(pairs)==7 and len(rejected)==45
    expected={frozenset([(1,3),(12,3)]),frozenset([(10,3),(12,3)]),frozenset([(17,3),(12,3)])}
    for i in[1,10]:
        for m in[3,4]:expected.add(frozenset([(i,3),(17,m)]))
    assert {frozenset(map(tuple,row['contacts']))for row in pairs}==expected
    return enc(dict(status='PASS_ALL_LABEL_S75_PROJECTION_CONTACT_GATE',sources={SOURCE:sha(SOURCE),f'scripts/{TAG}.py':sha(f'scripts/{TAG}.py')},
        support=S,gram=G,inverse=J,PD_pivots=piv,inverse_diagonals=dict(zip(S,diag)),
        exact_normalized_projection_condition='c^T J c <1/2, where c_i=-(r dot s_i)/sqrt(2q)>=0',
        diagonal_cosine_bound='Every nonright finite contact must satisfy cos^2(pi/m)<1/(2J_ii). Any J_ii>=2 forces right.',
        forced_right_walls=forced,wall75_diagonal=5,wall75_must_be_right=True,
        individual_label_domains=domains,three_contact_cost_lower_bound=min3,maximum_nonright_support_contacts=2,
        no_dotted_or_parallel_support_contacts=True,
        paired17_all_labels_ge5_cost_lower_bound=paired_tail_bound,
        zero_contact_profile=True,single_finite_profiles=singles,two_contact_profiles=pairs,rejected_two_contact_candidates=rejected,
        isolated_label_tail=dict(wall=17,labels='every finite m>=3',cost='cos^2(pi/m)/2 <1/2'),
        finite_profiles_including_zero=17,unbounded_single_label_tails=1,
        scope='Complete necessary support-contact classification for an acute Coxeter normal whose projection onto the specified retainedS75 Lorentz3 complement is spacelike. Does not impose other98-wall gates or claim any profile has an admissible ambient realization.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
