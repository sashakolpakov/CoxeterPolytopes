#!/usr/bin/env python3
"""A single pi/4 third contact with shared or opposite single pi/3 anchors."""
import argparse
import json
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_star_third_four_20261015'
PRE='fresh_h9_star_single_third_opposite_20261015'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
BATCH=f'results/{PRE}_search_batch_014.json'
S=[3,4,6,7,9,10,12,13]
CORE=sorted(S+[11,14,17])
R2,R3,R5,R10=map(s.sqrt,(2,3,5,10))

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(e):
    e=s.expand(e);b,c,d=[e.coeff(r)for r in(R2,R5,R10)]
    a=s.expand(e-b*R2-c*R5-d*R10)
    assert all(v.is_Rational for v in(a,b,c,d)),e
    return list(map(str,(a,b,c,d)))
def sign5(a,b):
    if not a:return(b>0)-(b<0)
    if not b or(a>0)==(b>0):return(a>0)-(a<0)
    z=a*a-5*b*b;return((a>0)-(a<0))*((z>0)-(z<0))
def sign(e):
    a,b,c,d=map(Fraction,enc(e));sa,sb=sign5(a,c),sign5(b,d)
    if not sa:return sb
    if not sb or sa==sb:return sa
    return sa*sign5(a*a+5*c*c-2*b*b-10*d*d,2*a*c-4*b*d)

def run():
    seen={}
    def read(p):
        if p in seen:return seen[p]
        d=json.loads((ROOT/p).read_text());seen[p]=d
        for q,h in d.get('sources',{}).items():
            assert sha(q)==h,q
            if q.endswith('.json'):read(q)
        return d
    prior=read(f'results/{PRE}.json');review=read(f'results/{PRE}_review.json')
    assert prior['status']=='PASS_ALL_LABEL_SINGLE_THIRD_OPPOSITE_ANCHOR_CORE'
    assert review['status']=='PASS_INDEPENDENT_OPPOSITE_SINGLE_THIRD_DOMAINS_AND_POLYNOMIALS'
    assert review['sources'][f'results/{PRE}.json']==sha(f'results/{PRE}.json')
    geom=read(GEOM)['representatives'][0]['geometry'];batch=read(BATCH)
    assert batch['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR'and batch['cumulative_additional_nodes']==60015
    V=[set(v['incident'])for v in geom['vertices']];MN=set(map(frozenset,geom['minimal_nonfaces']))
    assert len(V)==378 and any(set(S)<=v for v in V)
    assert all(any(set(p)<=v for v in V)for p in combinations(CORE,2))
    maximal={frozenset(v&set(CORE))for v in V}
    maximal=sorted([sorted(v)for v in maximal if not any(v<w for w in maximal)])
    assert maximal==prior['complete_maximal_core_faces']
    def face(nodes):
        return dict(support=sorted(nodes),actual_vertex=next(i for i,v in enumerate(V)if set(nodes)<=v))
    all_guards=[];budgets=[];domains=[];families=[];orientations=[]
    for shared in(False,True):
        own={11:13,17:13 if shared else 12}
        base={tuple(p):m for p,m in zip(((3,4),(6,7),(9,10),(12,13)),(5,7,7,7))}
        base[4,14]=4
        for a in(11,17):base[tuple(sorted((a,own[a])))]=3
        guards=[]
        for old in prior['positive_kernel_guards']:
            nodes=list(old['ordered_support'])
            if shared and old['kind']in('mutual_and_third_contact_forbidden','mutual_and_owner_foreign_forbidden','nonright_owner_foreign_forces_other_third_right'):
                nodes=[13 if x==12 else x for x in nodes]
            assert len(set(nodes))==len(nodes)
            facts=dict(base)
            for name,value in old['conditional'].items():
                if name=='trial_owner_internal_at_least':facts[3,4]=value
                elif name=='trial_mutual_label_at_least':facts[11,17]=value
                elif name in('foreign_contacts','third_contacts'):
                    for p in value:facts[tuple(sorted(p))]=3
                else:
                    assert name in('trial_nonright','trial_label_at_least_four','owner_foreign','third_contact','trial_mutual_nonright','foreign_contact','trial_third_nonright')
                    facts[tuple(sorted(value))]=4 if name=='trial_label_at_least_four'else 3
            G=s.eye(len(nodes));edges=old['edges'];cos={3:s.Rational(1,2),4:R2/2,6:R3/2}
            for i,j,m in edges:
                assert facts.get(tuple(sorted((nodes[i],nodes[j]))),0)>=m
                G[i,j]=G[j,i]=-cos[m]
            w=s.Matrix([s.sympify(x)for x in old['positive_kernel']])
            assert all(x.is_positive for x in w)and all(s.expand(x)==0 for x in G*w)
            guards.append(dict(anchor_mode='shared'if shared else'opposite',kind=old['kind'],ordered_support=nodes,
                edges=edges,positive_kernel=list(map(str,w)),conditional=old['conditional'],
                available_edge_lower_bounds=[list(p)+[m]for p,m in sorted(facts.items())],**face(nodes)))
        assert len(guards)==54;all_guards+=guards
        for a in(11,17):
            budgets.append(dict(anchor_mode='shared'if shared else'opposite',leaf=a,removed_own_root=25-own[a],
                retained_own_root=own[a],own_cost='1/4',**face((set(S)-{25-own[a]})|{a})))
        budgets.append(dict(anchor_mode='shared'if shared else'opposite',leaf=14,removed_own_root=3,retained_own_root=4,
            own_cost='1/2',foreign_pi3_cost_lower='1/2',sum_lower='1',**face((set(S)-{3})|{14})))
        literal=[]
        for fx,fy,dx,dy,mu in product((None,3,4),(None,3,4),(2,3),(2,3),(2,3,4,5)):
            if fx is not None and fy is not None:continue
            if dx==dy==3:continue
            if mu>2 and(fx is not None or fy is not None or dx==3 or dy==3):continue
            if(fx==4 or fy==4)and(dx==3 or dy==3):continue
            literal.append([fx,fy,dx,dy,mu,None if fx is None and fy is None else 5])
        assert len(literal)==14
        reduced=set()
        for fx,fy,dx,dy,mu,m in literal:
            if fy is not None:fx,fy,dx,dy=fy,fx,dy,dx
            reduced.add((None if fx is None else fx-3,int(dx==3),int(dy==3),mu,m))
        expected={(None,0,0,m,None)for m in(2,3,4,5)}|{(None,1,0,2,None),(None,0,1,2,None)}
        expected|={(0,a,b,2,5)for a,b in((0,0),(1,0),(0,1))}|{(1,0,0,2,5)}
        assert reduced==expected and len(reduced)==10
        domains.append(dict(anchor_mode='shared'if shared else'opposite',literal_patterns=literal,
            canonical_patterns=sorted([list(z)for z in reduced],key=str)))
        for orev,lrev in product((False,True),repeat=2):
            perm=list(range(18))
            if orev:perm[3],perm[4]=4,3
            if lrev:perm[12],perm[13]=13,12
            assert{frozenset(perm[x]for x in v)for v in V}==set(map(frozenset,V))
            assert{frozenset(perm[x]for x in m)for m in MN}==MN
            fixed=[[3,14,2],[4,14,4],[11,12,2],[11,13,3],
                   [12,17,2 if shared else 3],[13,17,3 if shared else 2]]
            bounds=[[3,4,5],[6,7,7],[9,10,7],[12,13,7]]
            image={name:sorted([sorted((perm[a],perm[b]))+[m]for a,b,m in rows])for name,rows in(('fixed',fixed),('bounds',bounds))}
            families.append(dict(id=('shared'if shared else'opposite')+'_third_four'+('_owner_swap'if orev else'')+('_last_swap'if lrev else''),origin='single_pi4_third_rank_core',**image))
            orientations.append(dict(anchor_mode='shared'if shared else'opposite',permutation=perm))
    assert len(all_guards)==108 and len(budgets)==6 and all(len(r['support'])==8 for r in budgets)
    # Exact no-foreign sign identity. Both third-anchor contacts cannot be nonzero.
    D,mu,a,b,h,v=s.symbols('D mu a b h v',real=True)
    # h is the anchor last-block projection factor, with h=1 or h=t.
    aa=1-D/4;off=-mu-h*D/4
    M=s.Matrix([[aa,off,-a],[off,aa,-b],[-a,-b,v]])
    identity=s.expand(M.det()-(v*(aa*aa-off*off)-aa*(a*a+b*b)+2*off*a*b))
    assert identity==0
    t,c=s.symbols('t c',real=True)
    anchor_shared=s.expand((aa*aa-off*off).subs(h,1))
    anchor_opposite=s.cancel((aa*aa-off*off).subs({h:t,D:1/(1-t*t)}))
    assert s.expand(anchor_shared-(1-D/2-mu*D/2-mu*mu))==0
    assert s.cancel(anchor_opposite-(1-7/(16*(1-t*t))-mu*t/(2*(1-t*t))-mu*mu))==0
    c5=(1+R5)/4;beta=R2/2
    assert sign(c5*c5-s.Rational(1,2))>0
    assert s.Rational(100,19)>4 and 1-s.Rational(7,16)*s.Rational(100,19)<0
    # Reprove the exact integer tail intervals independently of saved conclusions.
    u=s.symbols('u');C=8*u**3+4*u*u-4*u-1
    assert s.expand(s.chebyshevt(7,u)-1-(u-1)*C*C)==0
    lo=s.Rational(31,50);hi=2*s.Rational(901,1000)**2-1
    assert C.subs(u,lo)<0<C.subs(u,hi)
    assert s.diff(C,u).subs(u,s.Rational(1,2))>0 and s.diff(C,u,2).subs(u,s.Rational(1,2))>0
    assert(1+lo)/2==s.Rational(81,100)
    assert sign((2+R2)/4-s.Rational(23,25)**2)>0
    # Four foreign formulas for each anchor pattern; a foreign contact pins owner5.
    q=s.expand(1-c5*c5);D=1/(1-t*t);cases=[];positive=negative=0
    for shared in(False,True):
        for foreign,dx,dy in((0,0,0),(0,0,1),(0,1,0),(1,0,0)):
            F=1/(4*q);T=beta*(c5 if foreign==0 else 1)/(2*q)
            R=s.Matrix([[1-D/4-F,-(1 if shared else t)*D/4,-s.Rational(dx,2)-T],
                        [-(1 if shared else t)*D/4,1-D/4,-s.Rational(dy,2)],
                        [-s.Rational(dx,2)-T,-s.Rational(dy,2),1-beta*beta/q]])
            P=s.simplify(s.cancel(q*(1-t*t)*R.det()));P=s.expand(P)
            assert s.Poly(P,t).degree()<=2
            exceptional=not shared and(foreign,dx,dy)==(0,0,1)
            if exceptional:
                old=[r for r in prior['continuous_cases']if r['branch']=='third5_last7'];assert len(old)==1
                assert s.expand(P-s.sympify(old[0]['scaled_determinant'],locals={'beta':c5,'c':R2/2,'t':t}))==0
            for branch in(['last7','last_ge8']if exceptional else['all_last_ge7']):
                low,high,target=s.Rational(9,10),s.Integer(1),1
                if branch=='last7':high,target=s.Rational(901,1000),-1
                if branch=='last_ge8':low=s.Rational(23,25)
                z=s.symbols('z');Q=s.Poly(s.expand(P.subs(t,low+(high-low)*z)),z)
                values=[s.expand(sum(Q.nth(j)*s.binomial(k,j)/s.binomial(2,j)for j in range(k+1)))for k in range(3)]
                assert all(sign(x)==target for x in values)
                if target==1:positive+=3
                else:negative+=3
                cases.append(dict(anchor_mode='shared'if shared else'opposite',foreign_owner_root=foreign,
                    third_anchor_cosines=[str(s.Rational(dx,2)),str(s.Rational(dy,2))],owner_internal=5,
                    scaled_determinant=str(P),branch=branch,interval=[str(low),str(high)],target_sign=target,
                    bernstein_coefficients=[enc(x)for x in values]))
    assert len(cases)==9 and positive==24 and negative==3
    labels={tuple(r[:2]):r[2]for r in batch['solution']['pair_labels']}
    assert any(all(labels[a,b]==m for a,b,m in f['fixed'])and all(labels[a,b]>=m for a,b,m in f['bounds'])for f in families)
    paths=[f'scripts/{TAG}.py',GEOM,BATCH,f'results/{PRE}.json',f'results/{PRE}_review.json']
    return dict(status='PASS_SINGLE_PI4_THIRD_SHARED_OR_OPPOSITE_ANCHORS',sources={p:sha(p)for p in paths},case=0,
        rank_core=CORE,complete_maximal_core_faces=maximal,positive_kernel_guards=all_guards,projection_budgets=budgets,
        contact_domains=domains,no_foreign_identity=dict(determinant='v*anchor_minor - aa*(delta_x^2+delta_y^2) + 2*off*delta_x*delta_y',
            anchor_diagonal='aa=1-D/4<0',shared_minor=str(anchor_shared),opposite_minor='1-7D/16-mu*t*D/2-mu^2',
            third_diagonal='v=1-1/(2*(1-c^2))<0 for c>=cos(pi/5)',delta_product='0',
            sign='anchor minor<0, v<0, aa<0: determinant strictly positive'),
        integer_tail=dict(cubic=str(C),cos7_interval=['9/10','901/1000'],cos8_lower='23/25',D_lower='100/19'),
        foreign_cases=cases,positive_coefficients=positive,negative_coefficients=negative,
        exceptional_formula_matches_previous=True,literal_orientations=orientations,families=families,covered_saved_cursor=60015,
        scope='Compact star0 core: owner internal>=5; other three internals>=7; third own(2,4); anchors shared or opposite singlepi3. All noncore labels and distances unrestricted. No whole-incidence exclusion or search continuation.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],len(data['positive_kernel_guards']),len(data['foreign_cases']),data['positive_coefficients'],data['negative_coefficients'])
