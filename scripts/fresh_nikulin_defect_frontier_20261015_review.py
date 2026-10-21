#!/usr/bin/env python3
"""Independent physical radical Grams, Fraction KKT and subunit-cost audit.

No producer/previous-review arithmetic imports. Published/global geometry and
previously reviewed all-rank coordinate-domain coverage are inherited explicitly.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_nikulin_defect_frontier_20261015';TAG=BASE+'_review'
OLD='fresh_nikulin_height_charge_20261015';UNIT='fresh_nikulin_unit_charge_20261015'
B=s.Rational(21,20)

def star(lengths):
    n=sum(lengths)+1;G=2*s.eye(n);start=1
    for length in lengths:
        path=[0]+list(range(start,start+length));start+=length
        for a,b in zip(path,path[1:]):G[a,b]=G[b,a]=-1
    return G

def grams():
    return {'A1':s.Matrix([[2,-2],[-2,2]]),'E6':star([2,2,2]),'E7':star([1,3,3]),'E8':star([1,2,5]),
      'F4':s.Matrix([[2,-1,0,0,0],[-1,2,-1,0,0],[0,-1,2,-s.sqrt(2),0],[0,0,-s.sqrt(2),2,-1],[0,0,0,-1,2]]),
      'G2':s.Matrix([[2,-1,0],[-1,2,-s.sqrt(3)],[0,-s.sqrt(3),2]])}

def omissions(name,forbidden):
    if forbidden is None:return {'A1':(0,1),'E6':(2,4),'E7':(3,6),'E8':(1,6),'F4':(0,4),'G2':(0,2)}[name]
    if name=='A1':return (1-forbidden,1-forbidden)
    if name=='E6':return tuple(v for v in (2,4,6)if v!=forbidden)[:2]
    if name=='E7':return (2,5)if forbidden in(4,7)else(4,7)
    if name=='E8':return (1,6)if forbidden in(3,7)else(3,7)
    if name=='F4':return (0,4)if forbidden in(1,3)else(1,3)
    return tuple(v for v in range(3)if v!=forbidden)

def subsets(nodes):
    for k in range(len(nodes)+1):yield from combinations(nodes,k)

def height(a,b):return s.simplify(max(min(a,b)**2,max(a,b)**2/B**2))
def qminimum(coeff):
    a,b,c=coeff;assert c>=0
    values=[a,a+b+c]
    if c and 0<=-b/(2*c)<=1:values.append(a-b*b/(4*c))
    return s.simplify(min(values))

def physical_exceptional(old,primary):
    rows=[];counts=Counter();defects=[]
    for name,G in grams().items():
        n=G.rows;d=G.nullspace()[0]
        if d[0]<0:d=-d
        d=s.simplify(d/min(d));assert all(x>=1 for x in d)
        assert G*d==s.zeros(n,1)and G.rank()==n-1
        P=s.simplify((G+d*d.T).inv())
        def R(i,j):
            v=s.zeros(n,1);v[i]+=1/d[i];v[j]-=1/d[j]
            assert (v.T*d)[0]==0
            return s.simplify((v.T*P*v)[0])
        inverses={}
        for omit in range(n):
            ids=[k for k in range(n)if k!=omit];Q=G.extract(ids,ids);J=s.simplify(Q.inv())
            assert all(Q[:k,:k].det()>0 for k in range(1,n))and all(x>=0 for x in J)
            assert all(J[k,k]>=s.Rational(1,2)for k in range(n-1))
            inverses[omit]=(ids,J)
        for forbidden in [None]+list(range(n)):
            om=omissions(name,forbidden);assert forbidden not in om
            T=s.Rational(2*(n-1-int(forbidden is not None)),9)
            saved=next(r for r in old['exceptional_cases']if r['type']==name and r['forbidden']==forbidden)
            widened=next(r for r in primary['widened_165_checks']['exceptional_records']if r['type']==name and r['forbidden']==forbidden)
            residual=next(r for r in primary['residual_exceptional_cases']if r['type']==name and r['forbidden']==forbidden)
            assert saved['omissions']==list(om)and s.Rational(saved['target'])==T
            assert saved['gram']==[[str(v)for v in row]for row in G.tolist()]
            assert list(map(s.sympify,saved['kernel']))==list(d)
            point=s.simplify(min(d[i]for i in om)**2*R(*om));assert point>=T
            cheap=[];subunit=[]
            for side,i in enumerate(om):
                ids,J=inverses[i]
                for k in ids:
                    c=J[ids.index(k),ids.index(k)]
                    if c<T:cheap.append((side,k,c))
                    if c<1:subunit.append((side,k,c))
            assert len(cheap)==saved['all_cheap_contact_count']
            def one(side,k):
                i,j=om[side],om[1-side]
                coeff=[R(i,j),-R(i,j)+R(k,j)-R(i,k),R(i,k)]
                return [s.simplify(height(d[k],d[j])*R(k,j)),s.simplify(height(d[i]+d[k],d[j])*qminimum(coeff))],coeff
            singles=[];multiple=[]
            for group in subsets(cheap):
                if not group:continue
                cost=sum((c for side,k,c in group),s.S(0))
                for (a,k,_),(b,l,_)in combinations(group,2):
                    if a==b:
                        ids,J=inverses[om[a]];cost+=2*J[ids.index(k),ids.index(l)]
                cost=s.simplify(cost)
                if cost>=T:continue
                if len(group)==1:
                    side,k,_=group[0];lo,coeff=one(side,k);assert min(lo)+cost>=T
                    raw=next(z for z in widened['singletons']if(z['side'],z['node'])==(side,k))
                    assert raw['energy_lower']==list(map(str,lo))and s.Rational(raw['cost'])==cost
                    singles.append(dict(side=side,node=k,cost=str(cost),energy_coefficients=list(map(str,coeff)),energy_lower=list(map(str,lo))))
                else:
                    assert len(group)==2
                    nodes=[(side,k)for side,k,_ in group]
                    special=name=='E7'and forbidden is None and nodes==[(0,4),(1,7)]
                    if special:
                        u,v=s.symbols('u v',real=True);z=s.zeros(n,1)
                        for i,a in [(3,1-u),(4,u),(6,v-1),(7,-v)]:z[i]=a/d[i]
                        energy=s.expand((z.T*P*z)[0]);assert energy==(1+u*u+v*v)/2 and cost==1
                        assert d[4]==d[7]==1 and d[3]+d[4]>=3 and d[6]+d[7]>=3
                        assert min(s.Rational(5,2),s.Rational(9,2)/B**2+1)>T
                        counts['special_double']+=1
                    else:
                        assert name in('E7','E8')and cost>=s.Rational(3,2)and T-cost<=s.Rational(2,9)
                        kind='ordinary_E8'if name=='E8'and forbidden is None else 'E7_type'
                        assert T==(s.Rational(16,9)if kind=='ordinary_E8'else s.Rational(14,9))
                        assert cost>=(s.Rational(14,9)if kind=='ordinary_E8'else s.Rational(3,2))
                        assert T-cost<=(s.Rational(2,9)if kind=='ordinary_E8'else s.Rational(1,18))
                        defects.append(dict(type=name,forbidden=forbidden,nodes=[list(z)for z in nodes],kind=kind,cost=str(cost),target=str(T),defect=str(T-cost)))
                    multiple.append(dict(nodes=[list(z)for z in nodes],cost=str(cost),special=special))
            assert len(singles)==len(saved['single_cases'])and len(multiple)==len(saved['multiple_cases'])
            subrecords=[]
            for side,k,cost in subunit:
                if T<=s.Rational(2,3):lower=s.S(0);reason='small target'
                else:
                    lower=min(one(side,k)[0]);assert lower>=T;reason='energy at least target'
                raw=next(z for z in residual['all_diagonals_below_one']if(z['side'],z['node'])==(side,k))
                assert s.Rational(raw['diagonal'])==cost and s.Rational(raw['energy_lower'])==lower
                subrecords.append(dict(side=side,node=k,cost=str(cost),energy_lower=str(lower),loss_upper=str(T-lower),reason=reason))
            assert len(subrecords)==len(residual['all_diagonals_below_one'])
            rows.append(dict(type=name,forbidden=forbidden,omissions=list(om),target=str(T),point_energy=str(point),singletons=singles,multiple=multiple,subunit_coordinates=subrecords))
            counts.update(records=1,singletons=len(singles),multiple=len(multiple),subunit_coordinates=len(subrecords))
    assert counts['records']==40 and counts['singletons']==71 and counts['multiple']==16 and len(defects)==15
    return dict(records=rows,defect_occurrences=defects,counts=dict(counts))

# Separate rational cycle Gram/dual replay, independent of all producer helpers.
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F())
def cycle_replay(old,primary):
    out=[];count=0;Bf=F(21,20)
    for m in range(2,18):
        h=m//2;T=F(2*(m-1),9)
        resistance=lambda i,j:F(((i-j)%m)*(m-(i-j)%m),m)
        K=lambda i,j:F(min(i,j))-F(i*j,m)
        assert all(K(i,i)+K(j,j)-2*K(i,j)==resistance(i,j)for i,j in product(range(m),repeat=2))
        potentials=[[F(int(k<i))-F(i,m)for k in range(m)]for i in range(m)]
        assert all(dot(potentials[i],potentials[j])==K(i,j)for i,j in product(range(m),repeat=2))
        domains=[]
        for omit in(0,h):
            nodes=[i for i in range(m)if i!=omit and resistance(i,omit)<=2];rows=[]
            for ids in subsets(nodes):
                cost=sum((resistance(i,omit)for i in ids),F())+sum((resistance(i,omit)+resistance(j,omit)-resistance(i,j)for i,j in combinations(ids,2)),F())
                if cost<=2:rows.append((ids,cost))
            domains.append(rows)
        raw=next(z for z in old['classical_A']['small_cycles']if z['cycle_size']==m)
        expected={}
        for (left,a),(right,b)in product(*domains):
            if a+b>=T:continue
            for i,j in product((False,True),repeat=2):
                U=left+((0,)if i else());V=right+((h,)if j else())
                if U and V:expected[(U,V)]=(a,b)
        assert len(expected)==len(raw['cases']);records=[]
        for rec in raw['cases']:
            U,V=tuple(rec['left_support']),tuple(rec['right_support']);a,b=expected.pop((U,V));n=len(U);ids=U+V
            signs=[1]*n+[-1]*len(V);G=[[signs[i]*signs[j]*K(k,l)for j,l in enumerate(ids)]for i,k in enumerate(ids)]
            d=rec['convex_dual'];x=list(map(F,d['barycentric_weights']));lm=list(map(F,d['equality_multipliers']))
            assert min(x)>=0 and sum(x[:n])==sum(x[n:])==1
            gx=[dot(row,x)for row in G];slack=[z+lm[int(i>=n)]for i,z in enumerate(gx)]
            assert min(slack)>=0 and all(z*t==0 for z,t in zip(x,slack))
            energy=dot(x,gx);assert energy==F(d['minimum_energy'])and energy==-sum(lm)
            H=max(min(len(U),len(V))**2,F(max(len(U),len(V))**2)/Bf**2)
            lower=H*energy+a+b;assert lower>=T
            stored=next(z for z in primary['widened_165_checks']['cycle_cases']if z['cycle_size']==m and z['left']==list(U)and z['right']==list(V))
            assert F(stored['unit_lower'])==lower
            records.append(dict(left=list(U),right=list(V),energy=str(energy),height_square=str(H),cost=str(a+b),lower=str(lower)))
        assert not expected;count+=len(records);out.append(dict(cycle_size=m,domain_counts=list(map(len,domains)),cases=records))
    assert count==580
    return dict(total=count,cycles=out)

def all_rank_A():
    m,t,u=s.symbols('m t u',positive=True);records=[]
    # For d>=2, concavity in d gives its minimum at2 (d<=m/2).
    d=s.symbols('d',real=True)
    assert s.expand(d*(m-d)-2*(m-2)-(d-2)*(m-d-2))==0
    for parity,start in ((0,2),(1,3)):
        h=(m-parity)/2;R=lambda z:z*(m-z)/m
        zero_margin=s.factor((R(h)-2*(m-1)/9).subs(m,start+2*t));zn,zd=s.fraction(zero_margin)
        assert all(v>=0 for v in s.Poly(zn,t).all_coeffs())and all(v>=0 for v in s.Poly(zd,t).all_coeffs())and zd.subs(t,0)>0
        point=s.factor(R(h-1));point_margin=s.factor(point-2*(m-1)/9+s.Rational(4,9))
        num,den=s.fraction(s.factor(point_margin.subs(m,start+2*t)))
        assert all(v>=0 for v in s.Poly(num,t).all_coeffs())and all(v>=0 for v in s.Poly(den,t).all_coeffs())and den.subs(t,0)>0
        ell=(m+parity)/m;poly=R(h)-ell*u+R(1)*u*u
        minimizer=s.factor(ell/(2*R(1)));minimum=s.factor(poly.subs(u,minimizer))
        difference=s.factor(minimum-(m-3)/4)
        # m>=4 even and m>=5 odd; minimizer belongs to[0,1].
        base=4+parity;sub={m:base+2*t}
        for z in (minimizer,1-minimizer,difference):
            a,b=s.fraction(s.factor(z.subs(sub)))
            assert all(v>=0 for v in s.Poly(a,t).all_coeffs())and all(v>=0 for v in s.Poly(b,t).all_coeffs())and b.subs(t,0)>0
        records.append(dict(parity=parity,point_energy=str(point),point_loss_margin=str(point_margin),segment_polynomial=str(s.expand(poly)),segment_minimizer=str(minimizer),minimum=str(minimum),difference_from_uniform=str(difference)))
    margin=s.expand(((m-3)/B**2-2*(m-1)/9).subs(m,4+t));assert margin==s.Rational(106,441)+s.Rational(302,441)*t
    # Recheck the inherited large-rank unit-charge tails without balance input.
    tails=[]
    for parity,start in((0,18),(1,19)):
        mm=start+2*t;R=mm/4 if not parity else (mm*mm-1)/(4*mm)
        margin2=s.factor(R-2*(mm-1)/9-4*(2+parity)/(mm-2));a,b=s.fraction(margin2)
        assert all(v>=0 for v in s.Poly(a,t).all_coeffs())and a.subs(t,0)>0
        assert all(v>=0 for v in s.Poly(b,t).all_coeffs())and b.subs(t,0)>0
        tails.append(dict(parity=parity,start=start,unit_margin=str(margin2)))
    return dict(single_contact_cases=records,positive_omitted_margin=str(margin),inherited_unit_tail_recheck=tails,avoid_node='Same opposite-omission cycle calculation by rotation; coincident avoid-A1 target0.')

def all_rank_BCD(unit):
    r,k=s.symbols('r k',integer=True,positive=True);t,u,v=s.symbols('t u v',integer=True,nonnegative=True)
    syms={str(z):z for z in(r,k,t,u,v)};rows=[]
    for raw in unit['classical_all_rank_checks']['prefix_coordinate_cases']:
        ds=s.sympify(raw['null_coefficient_squared'])
        # Parse symbols consistently, and check the inverse/loss identities.
        R=[s.sympify(z,locals=syms)for z in raw['point_resistances']]
        cost=s.sympify(raw['inverse_diagonal'],locals=syms);loss=s.sympify(raw['linear_loss'],locals=syms)
        assert s.simplify(cost-ds*R[1])==0 and s.simplify(loss-R[0]-R[1]+R[2])==0
        T=2*(r-int(raw['modified']))/9
        sub={syms[a]:s.sympify(b,locals=syms)for a,b in raw['domain_substitution'].items()}
        cost0=s.expand(cost.subs(sub,simultaneous=True));loss0=s.expand(loss.subs(sub,simultaneous=True));target=s.expand(T.subs(sub,simultaneous=True))
        point_margin=s.expand((R[0]-T).subs(sub,simultaneous=True));assert all(z>=0 for z in s.Poly(point_margin,t,u,v).coeffs())
        # The unchanged unit-charge tail has either ell<=cost or cost>=T.
        margins=[s.Poly(s.expand((cost-loss).subs(sub,simultaneous=True)),t,u,v),s.Poly(s.expand((cost-T).subs(sub,simultaneous=True)),t,u,v)]
        assert any(all(z>=0 for z in p.coeffs())for p in margins)
        coeff=[cost0.subs({t:0,u:0,v:0})]+[cost0.coeff(z)for z in(t,u,v)]
        assert cost0==coeff[0]+sum(a*z for a,z in zip(coeff[1:],(t,u,v)))and min(coeff)>=0
        domains=[];small=[]
        if coeff[0]>=1:reason='all diagonals at least1'
        elif loss0==0:reason='zero linear energy loss'
        else:
            for z,c in zip((t,u,v),coeff[1:]):
                if c==0:assert target.coeff(z)==0;domains.append([0])
                else:
                    bound=int(s.ceiling((1-coeff[0])/c));domains.append(range(bound))
                    assert coeff[0]+c*bound>=1
            for point in product(*domains):
                subst=dict(zip((t,u,v),point))
                if cost0.subs(subst)<1:
                    val=target.subs(subst);assert val<=s.Rational(2,3);small.append(dict(parameters=list(point),target=str(val),cost=str(cost0.subs(subst))))
            assert small;reason='integer rank bound gives target at most2/3'
        rows.append(dict(type=raw['type'],modified=raw['modified'],side=raw['side'],node_class=raw['node_class'],cost=str(cost0),linear_loss=str(loss0),reason=reason,small_cases=small))
    assert len(rows)==24
    # Independent norm-two physical B3 sharpness and modified B2 check.
    G=s.Matrix([[2,0,-1,0],[0,2,-1,0],[-1,-1,2,-s.sqrt(2)],[0,0,-s.sqrt(2),2]])
    J=G[:3,:3].inv();assert J[0,0]==s.Rational(3,4)
    d=G.nullspace()[0];d=s.simplify(d/min(d));assert d[0]==1
    smallB=s.Matrix([[2,0,-s.sqrt(2)],[0,2,-s.sqrt(2)],[-s.sqrt(2),-s.sqrt(2),2]])
    for i in range(3):
        ids=[j for j in range(3)if i!=j];Q=smallB.extract(ids,ids).inv();assert min(Q[a,a]for a in range(2))>=s.Rational(1,2)>s.Rational(2,9)
    return dict(classes=rows,separate_modified_B2_target='2/9',sharp_B3=dict(physical_gram=[[str(z)for z in row]for row in G.tolist()],omissions=[0,3],both_rows_node=0,sole_cost='3/4',height_ratio=1,energy=0,target='2/3'))

def aggregation(primary,exceptional):
    group=Counter(z['kind']for z in exceptional['defect_occurrences']);assert group==Counter(ordinary_E8=4,E7_type=11)
    rows=[]
    for large in range(3):
        minimum=large*F(14,9)+(2-large)*F(3,2);target=large*F(16,9)+(2-large)*F(14,9);residual=4-minimum
        if large:assert residual<1;loss=target+F(2,3)
        else:loss=target+residual
        assert loss<=F(38,9)
        rows.append(dict(ordinary_E8=large,cost_lower=str(minimum),residual_budget=str(residual),bad_target=str(target),loss=str(loss)))
    assert 3*F(3,2)>4 and 4+F(2,9)==F(38,9)
    t=s.symbols('t',nonnegative=True);near=s.expand((2*(164+t)/9-s.Rational(2,9)-s.Rational(38,9))/2-2);far=s.expand(B*(164+t-43)/9)
    assert near==14+t/9 and far==14+s.Rational(7,60)+s.Rational(7,60)*t
    assert primary['total_loss_at_most']=='38/9'and primary['threshold_checks']['combined_candidate_bound']==164
    return dict(defect_occurrences=dict(group),two_defect_cases=rows,zero_or_one_defect_loss='38/9',near_at_N164_plus_t=str(near),far_at_N164_plus_t=str(far),ideal_branch_bound=164,inherited_compact_branch=143,combined_candidate=164)

def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def run():
    old=json.loads((ROOT/f'results/{OLD}.json').read_text());primary=json.loads((ROOT/f'results/{BASE}.json').read_text());unit=json.loads((ROOT/f'results/{UNIT}.json').read_text())
    exceptional=physical_exceptional(old,primary);cycles=cycle_replay(old,primary);A=all_rank_A();BCD=all_rank_BCD(unit);agg=aggregation(primary,exceptional)
    paths=[f'scripts/{TAG}.py']
    for name in(BASE,OLD,OLD+'_review',UNIT,UNIT+'_review',UNIT+'_dimension169','fresh_nikulin_coupled_charge_20261015_dimension166'):
        paths.extend([f'scripts/{name}.py',f'results/{name}.json'])
    visited=set()
    def bind(p):
        if p in visited:return
        visited.add(p)
        if p.endswith('.json'):
            for q,h in json.loads((ROOT/p).read_text()).get('sources',{}).items():assert digest(q)==h,(p,q);bind(q)
    for p in paths:bind(p)
    return dict(status='PASS_INDEPENDENT_PHYSICAL_GRAM_KKT_SUBUNIT_ENVELOPE_AND_164_SCOPE_REVIEW',sources={p:digest(p)for p in paths},recursive_manifest_files=len(visited),
                physical_exceptional=exceptional,widened_cycle_KKT=cycles,all_rank_A=A,all_rank_BCD=BCD,aggregate=agg,
                scope='Independent new arithmetic and hostile proof read. Complete literal exceptional support lists, widened580KKT, all-rank subunit A algebra, all24 inherited coordinate-class implications and new aggregation checked. Previously reviewed coordinate-domain coverage and global geometric169/166/165 hypotheses inherited; no publication/current-best claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text());print('PASS saved independent164 replay')
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['physical_exceptional']['counts'],out['widened_cycle_KKT']['total'])
