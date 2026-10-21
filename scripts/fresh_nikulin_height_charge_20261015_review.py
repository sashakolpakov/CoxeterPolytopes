#!/usr/bin/env python3
"""Independent rational grounded matrices/KKT checks for the 165 candidate.

No SymPy, producer, optimizer, or project arithmetic imports.  Complete support
lists are rebuilt; supplied convex duals are checked, not re-optimized.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from math import isqrt
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_height_charge_20261015_review'
PRIMARY='fresh_nikulin_height_charge_20261015'
B=F(31,30)

def dot(a,b):return sum((x*y for x,y in zip(a,b)),F())
def mv(A,x):return [dot(row,x)for row in A]
def sub(A,ids):return [[F(A[i][j])for j in ids]for i in ids]
def inverse(A):
    n=len(A);R=[[F(v)for v in row]+[F(i==j)for j in range(n)]for i,row in enumerate(A)]
    for k in range(n):
        h=next(i for i in range(k,n)if R[i][k]);R[k],R[h]=R[h],R[k]
        z=R[k][k];R[k]=[x/z for x in R[k]]
        for i in range(n):
            if i!=k:
                z=R[i][k];R[i]=[a-z*b for a,b in zip(R[i],R[k])]
    return [r[n:]for r in R]
def positive(A):
    A=[[F(x)for x in row]for row in A]
    while A:
        assert A[0][0]>0
        A=[[A[i][j]-A[i][0]*A[0][j]/A[0][0]for j in range(1,len(A))]for i in range(1,len(A))]
def sqrtq(x):
    a,b=isqrt(x.numerator),isqrt(x.denominator)
    assert F(a*a,b*b)==x
    return F(a,b)
def subsets(a):
    for n in range(len(a)+1):yield from combinations(a,n)
def star(arms):
    n=1+sum(arms);A=[[F(2*(i==j))for j in range(n)]for i in range(n)];k=1
    for length in arms:
        old=0
        for _ in range(length):A[old][k]=A[k][old]=-1;old=k;k+=1
    return A

def model(A,scales):
    A=[[F(x)for x in row]for row in A];scales=list(map(F,scales));n=len(A)
    I=list(range(1,n));J=inverse(sub(A,I));kernel=[F(1)]+mv(J,[-A[i][0]for i in I])
    assert min(kernel)>0 and mv(A,kernel)==[0]*n
    normalization=min(d*k*k for d,k in zip(scales,kernel))
    delta2=[d*k*k/normalization for d,k in zip(scales,kernel)]
    grams=[];diags=[]
    for omit in range(n):
        keep=[i for i in range(n)if i!=omit];Q=sub(A,keep);positive(Q);K=inverse(Q)
        assert all(x>=0 for row in K for x in row)
        grams.append((keep,K));diags.append({j:scales[j]*K[a][a]for a,j in enumerate(keep)})
    def resistance(i,j):
        f=[F()for _ in range(n)];f[i]+=1/kernel[i];f[j]-=1/kernel[j]
        v=[F()]+mv(J,f[1:]);assert mv(A,v)==f and dot(f,kernel)==0
        return normalization*dot(f,v)
    def cross(omit,i,j):
        keep,K=grams[omit]
        return sqrtq(scales[i]*scales[j])*K[keep.index(i)][keep.index(j)]
    return delta2,diags,resistance,cross

def physical_guard(A,scales,saved):
    for i in range(len(A)):
        for j in range(len(A)):
            v=saved[i][j]
            if v in ('-sqrt(2)','-sqrt(3)'):sg=-1;sq=F(2 if v=='-sqrt(2)'else 3)
            else:z=F(v);sg=(z>0)-(z<0);sq=z*z
            assert sg==(A[i][j]>0)-(A[i][j]<0) and sq==F(A[i][j])**2/(scales[i]*scales[j])
def kernel_square(text):
    if 'sqrt'not in text:return F(text)**2
    coeff,rad=text.split('sqrt(');coeff=coeff.rstrip('*')or'1'
    return F(coeff)**2*int(rad.rstrip(')'))
def hbound(a2,b2):return max(min(a2,b2),max(a2,b2)/B**2)
def minimum(c):
    a,b,d=c;assert d>=0;vals=[a,a+b+d]
    if d and 0<=-b/(2*d)<=1:vals.append(a-b*b/(4*d))
    return min(vals)
def choices(name,rank,forbidden):
    ordinary={'A1':(0,1),'E6':(2,4),'E7':(3,6),'E8':(1,6),'F4':(0,4),'G2':(0,2)}
    if forbidden is None:return ordinary[name]
    if name=='A1':return(1-forbidden,)*2
    if name=='E6':return tuple(i for i in(2,4,6)if i!=forbidden)[:2]
    if name=='E7':return(2,5)if forbidden in(4,7)else(4,7)
    if name=='E8':return(1,6)if forbidden in(3,7)else(3,7)
    if name=='F4':return(0,4)if forbidden in(1,3)else(1,3)
    return tuple(i for i in range(3)if i!=forbidden)

def exceptions(saved):
    literal={'A1':([[2,-2],[-2,2]],[1,1]),'E6':(star([2,2,2]),[1]*7),
      'E7':(star([1,3,3]),[1]*8),'E8':(star([1,2,5]),[1]*9),
      'F4':([[2,-1,0,0,0],[-1,2,-1,0,0],[0,-1,2,-2,0],[0,0,-2,4,-2],[0,0,0,-2,4]],[1,1,1,2,2]),
      'G2':([[2,-1,0],[-1,2,-3],[0,-3,6]],[1,1,3])}
    rows=[];single_count=multiple_count=0;defects=[]
    for name,(A,scales)in literal.items():
        ds,diags,R,cross=model(A,scales);r=len(A)-1
        for forbidden in[None]+list(range(r+1)):
            om=choices(name,r,forbidden);T=F(2*(r-int(forbidden is not None)),9)
            old=next(z for z in saved['exceptional_cases']if z['type']==name and z['forbidden']==forbidden)
            assert old['rank']==r and old['omissions']==list(om)and F(old['target'])==T
            assert forbidden not in om
            physical_guard(A,scales,old['gram']);assert list(map(kernel_square,old['kernel']))==ds
            assert old['inverse_diagonals']==[{str(k):str(v)for k,v in diags[i].items()}for i in om]
            point=min(ds[i]for i in om)*R(*om);assert point>=T and F(old['point_energy'])==point
            cheap=[(side,k,c)for side,i in enumerate(om)for k,c in diags[i].items()if c<T]
            assert old['all_cheap_contact_count']==len(cheap)and old['all_cheap_subsets_checked']==2**len(cheap)-1
            supports=[]
            for group in subsets(cheap):
                if not group:continue
                cost=sum((c for _,_,c in group),F())
                cost+=sum((2*cross(om[a],k,l)for(a,k,_),(b,l,_)in combinations(group,2)if a==b),F())
                if cost<T:supports.append((group,cost))
            assert len(supports)==old['low_cost_supports']
            singles=[];multiples=[]
            for group,cost in supports:
                if len(group)==1:
                    side,k,_=group[0];i,j=om[side],om[1-side]
                    coeff=[R(i,j),-R(i,j)+R(k,j)-R(i,k),R(i,k)];emin=minimum(coeff)
                    plus_square=ds[i]+ds[k]+2*sqrtq(ds[i]*ds[k])
                    lower0=hbound(ds[k],ds[j])*R(k,j)+cost
                    lower1=hbound(plus_square,ds[j])*emin+cost
                    assert min(lower0,lower1)>=T
                    raw=next(z for z in old['single_cases']if(z['side'],z['node'])==(side,k))
                    assert F(raw['cost'])==cost and F(raw['energy_min'])==emin
                    assert F(raw['zero_omitted_lower'])==lower0 and F(raw['positive_omitted_lower'])==lower1
                    singles.append(dict(side=side,node=k,cost=str(cost),energy_coefficients=list(map(str,coeff)),
                        exact_minimum=str(emin),lower_bounds=list(map(str,(lower0,lower1)))))
                else:
                    assert len(group)==2
                    nodes=[(side,k)for side,k,_ in group]
                    special=name=='E7'and forbidden is None and nodes==[(0,4),(1,7)]
                    raw=next(z for z in old['multiple_cases']if z['nodes']==list(map(list,nodes)))
                    assert F(raw['cost'])==cost
                    if special:
                        # Polarization at six rational points recovers a degree2
                        # polynomial, known from a squared Euclidean difference.
                        def energy(u,v):
                            ids=[3,4,6,7];ws=[1-u,u,-(1-v),-v]
                            return -sum((ws[a]*ws[b]*R(ids[a],ids[b])for a,b in combinations(range(4),2)),F())
                        e00=energy(F(),F());e10=energy(F(1),F());e01=energy(F(),F(1))
                        ex=2*(energy(F(1),F())-2*energy(F(1,2),F())+e00)
                        ey=2*(energy(F(),F(1))-2*energy(F(),F(1,2))+e00)
                        coeff=[e00,e10-e00-ex,e01-e00-ey,ex,energy(F(1),F(1))-e10-e01+e00,ey]
                        assert coeff==[F(1,2),0,0,F(1,2),0,F(1,2)]and cost==1
                        lower0=F(3,2)+cost;lower1=F(9,2)/B**2+cost
                        assert min(lower0,lower1)>T
                        assert F(raw['neither_omitted_lower'])==lower0 and F(raw['some_omitted_lower'])==lower1
                        assert raw['status']=='UNIT_CHARGE_VIA_E7_HEIGHT'
                        multiples.append(dict(nodes=nodes,cost=str(cost),energy_coefficients=list(map(str,coeff)),defect=0))
                    else:
                        assert cost>=F(3,2)and T-cost<=F(2,9)
                        assert raw['status']=='POSSIBLE_BUDGETED_DEFECT'and F(raw['defect_upper'])==T-cost
                        item=dict(type=name,forbidden=forbidden,nodes=nodes,cost=str(cost),defect_upper=str(T-cost))
                        defects.append(item);multiples.append(item)
            assert len(singles)==len(old['single_cases'])and len(multiples)==len(old['multiple_cases'])
            single_count+=len(singles);multiple_count+=len(multiples)
            rows.append(dict(type=name,forbidden=forbidden,omissions=om,point_lower=str(point),
                complete_cheap_subsets=2**len(cheap)-1,single_cases=singles,multiple_cases=multiples))
    assert len(rows)==40 and single_count==71 and multiple_count==16 and len(defects)==15
    return rows,defects

def cycles(saved):
    rows=[];total=0
    for m in range(2,18):
        h=m//2;T=F(2*(m-1),9)
        # Grounding node0 turns the cycle (including double-bond A1) into A_(m-1).
        J=inverse([[F(2 if i==j else -1 if abs(i-j)==1 else 0)for j in range(m-1)]for i in range(m-1)])
        K=[[F()]*m]+[[F()]+row for row in J]
        def R(a,b):return K[a][a]+K[b][b]-2*K[a][b]
        assert all(K[i][j]==min(i,j)-F(i*j,m)for i in range(m)for j in range(m))
        assert all(R(i,j)==F(((i-j)%m)*(m-(i-j)%m),m)for i in range(m)for j in range(m))
        domains=[]
        for omit in(0,h):
            nodes=[k for k in range(m)if k!=omit and R(k,omit)<=2];group=[]
            for support in subsets(nodes):
                cost=sum((R(k,omit)for k in support),F())
                cost+=sum((R(k,omit)+R(l,omit)-R(k,l)for k,l in combinations(support,2)),F())
                if cost<=2:group.append((support,cost))
            domains.append(group)
        old=next(r for r in saved['classical_A']['small_cycles']if r['cycle_size']==m)
        assert old['retained_domains']==[[[list(ids),str(c)]for ids,c in d]for d in domains]
        expected={}
        for (left,cp),(right,cq)in product(*domains):
            if cp+cq>=T:continue
            for a,b in product((False,True),repeat=2):
                U=left+((0,)if a else());V=right+((h,)if b else())
                if U and V:expected[(U,V)]=(cp,cq)
        assert len(expected)==len(old['cases'])
        records=[]
        for raw in old['cases']:
            U,V=tuple(raw['left_support']),tuple(raw['right_support']);cp,cq=expected.pop((U,V))
            assert F(raw['left_min_cost'])==cp and F(raw['right_min_cost'])==cq
            w=raw['convex_dual'];x=list(map(F,w['barycentric_weights']));lag=list(map(F,w['equality_multipliers']))
            nodes=U+V;sign=[1]*len(U)+[-1]*len(V)
            H=[[sign[i]*sign[j]*K[a][b]for j,b in enumerate(nodes)]for i,a in enumerate(nodes)]
            assert len(x)==len(nodes)and len(lag)==2 and min(x)>=0
            assert sum(x[:len(U)])==sum(x[len(U):])==1
            res=[v+lag[int(i>=len(U))]for i,v in enumerate(mv(H,x))]
            assert min(res)>=0 and all(a*b==0 for a,b in zip(x,res))
            assert list(map(F,w['dual_slacks']))==res
            energy=dot(x,mv(H,x));assert energy>=0 and energy==F(w['minimum_energy'])
            # Also check the dual objective identity; any feasible x0 has
            # q(x0)-q(x)=q(x0-x)+2*x0.res >=0.
            assert energy==-sum(lag)
            height=hbound(F(len(U)**2),F(len(V)**2));lower=height*energy+cp+cq
            assert height==F(raw['minimum_height_square'])and lower==F(raw['unit_charge_lower'])and lower>=T
            records.append(dict(left=U,right=V,energy=str(energy),cost=str(cp+cq),lower=str(lower),margin=str(lower-T)))
        assert not expected;total+=len(records)
        rows.append(dict(cycle_size=m,complete_domain_sizes=list(map(len,domains)),cases=records))
    assert total==580
    return rows

# Affine polynomials c+t*a+u*b+v*d, sufficient for all prefix-coordinate
# lengths, distances and unit-charge margins. No symbolic algebra dependency.
def lin(c=0,t=0,u=0,v=0):return tuple(map(F,(c,t,u,v)))
def la(a,b):return tuple(x+y for x,y in zip(a,b))
def ls(a,b):return tuple(x-y for x,y in zip(a,b))
def scale(a,c):return tuple(c*x for x in a)
def nonnegative(a):return all(x>=0 for x in a)
def lmax(a,b):
    if nonnegative(ls(a,b)):return a
    assert nonnegative(ls(b,a));return b
def lmin(a,b):
    if nonnegative(ls(a,b)):return b
    assert nonnegative(ls(b,a));return a

def classical_BCD():
    rows=[]
    def check(typ,modified,omit,other,node,ds,r,k):
        def idx(z):return r if z=='r'else ls(r,lin(1))if z=='r-1'else k if z=='k'else lin(z)
        def potential(z):
            if z==0:return[]
            if typ=='C':return[(lin(1),idx(z),F(1))]
            if z==1:return[(lin(1),lin(1),F(1))]
            if typ=='D'and z=='r-1':return[(lin(1),ls(r,lin(1)),F(1,2)),(r,r,F(-1,2))]
            return[(lin(1),idx(z),F(1,2))]
        def distance(a,b):
            segments=potential(a)+[(lo,hi,-c)for lo,hi,c in potential(b)];out=lin()
            for lo,hi,c in segments:
                for low,high,d in segments:
                    length=la(ls(lmin(hi,high),lmax(lo,low)),lin(1))
                    if nonnegative(length):out=la(out,scale(length,c*d))
                    else:assert nonnegative(scale(length,-1))
            return scale(out,F(1,2)if typ=='C'else F(1))
        Rij=distance(omit,other);Rik=distance(omit,node);Rkj=distance(node,other)
        cost=scale(Rik,F(ds));loss=ls(la(Rij,Rik),Rkj);T=scale(ls(r,lin(int(modified))),F(2,9))
        assert nonnegative(ls(Rij,T))
        direct=ls(cost,loss);far=ls(cost,T)
        if nonnegative(direct):kind='loss<=cost';margin=direct
        else:assert nonnegative(far)and far[0]>0;kind='cost>T';margin=far
        rows.append(dict(type=typ,modified=modified,omissions=[omit,other],node=node,
            rank_coefficients=list(map(str,r)),node_coefficients=list(map(str,k)),
            point_coefficients=list(map(str,Rij)),cost_coefficients=list(map(str,cost)),
            loss_coefficients=list(map(str,loss)),reason=kind,margin_coefficients=list(map(str,margin))))
    for i,j in((0,'r'),('r',0)):
        for node,ds,r,k in[(0,1,lin(2,t=1),lin()),(1,1,lin(2,t=1),lin()),('k',4,lin(3,u=1,v=1),lin(2,u=1)),('r',2,lin(2,t=1),lin())]:
            if node!=i:check('B',False,i,j,node,ds,r,k)
    for i,j in(('r-1',0),(0,'r-1')):
        items=[(0,1,lin(3,t=1),lin()),(1,1,lin(3,t=1),lin()),('r',2,lin(3,t=1),lin())]
        items+=[('k',4,lin(4 if i=='r-1'else 3,u=1,v=1),lin(2,u=1))]
        for node,ds,r,k in items:
            if node!=i:check('B',True,i,j,node,ds,r,k)
    check('C',False,0,'r','k',2,lin(2,u=1,v=1),lin(1,u=1))
    check('C',False,0,'r','r',1,lin(2,t=1),lin())
    for i,j in((1,'r'),('r',1)):
        items=[(0,1,lin(2,t=1),lin())]
        if i==1:items += [('k',2,lin(3,u=1,v=1),lin(2,u=1)),('r',1,lin(2,t=1),lin())]
        else:items += [('k',2,lin(2,u=1,v=1),lin(1,u=1))]
        for node,ds,r,k in items:check('C',True,i,j,node,ds,r,k)
    for node,ds,r,k in[(1,1,lin(4,t=1),lin()),('k',4,lin(4,u=1,v=1),lin(2,u=1)),('r-1',1,lin(4,t=1),lin()),('r',1,lin(4,t=1),lin())]:
        check('D',False,0,'r',node,ds,r,k)
    assert len(rows)==24
    A=[[2,0,-2],[0,2,-2],[-2,-2,4]];ds,diags,R,_=model(A,[1,1,2])
    assert R(0,1)==1 and min(c for i in(0,1)for c in diags[i].values())>=F(1,2)>F(2,9)
    return rows

def pmul(a,b):
    out=[F()]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out
def padd(a,b):return[(a[i]if i<len(a)else 0)+(b[i]if i<len(b)else 0)for i in range(max(len(a),len(b)))]
def pscale(a,k):return[k*x for x in a]
def tails():
    out=[]
    for m,eps in[([F(18),F(2)],0),([F(19),F(2)],1)]:
        square=pmul(m,m)
        if not eps:
            num=padd(padd(square,pscale(m,6)),[-304]);den=pscale(padd(m,[-2]),36)
        else:
            num=padd(padd(padd(pmul(square,m),pscale(square,6)),pscale(m,-457)),[18]);den=pscale(pmul(m,padd(m,[-2])),36)
        assert min(num)>=0 and num[0]>0 and min(den)>=0 and den[0]>0
        # Reconstruct numerator by cross multiplication from R(h)-T-4(2+eps)/(m-2).
        if not eps:
            point_minus_target=padd(m,[8]);derived=padd(pmul(point_minus_target,padd(m,[-2])),[-288])
        else:
            point_num=pmul(padd(m,[-1]),padd(m,[9]));derived=padd(pmul(point_num,padd(m,[-2])),pscale(m,-432))
        assert derived==num
        out.append(dict(cycle_size_coefficients=list(map(str,m)),parity=eps,numerator=list(map(str,num)),denominator=list(map(str,den))))
    return out

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    data=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    ex,defects=exceptions(data);cy=cycles(data);bc=classical_BCD();tail=tails()
    assert 3*F(3,2)>4
    loss=4+2*F(2,9);assert loss==F(40,9)
    # N=n-1; one avoid-node target costs2/9. Then divide by2 and subtract2.
    near=(F(2*165,9)-F(2,9)-loss)/2-2;far=B*F(165-43,9)
    assert near==14 and far==14+F(1,135)
    assert F(1,9)>0 and B/F(9)>0
    paths=[f'scripts/{TAG}.py',f'scripts/{PRIMARY}.py',f'results/{PRIMARY}.json']
    for stem in('fresh_nikulin_unit_charge_20261015_review','fresh_nikulin_unit_charge_20261015_dimension169',
                 'fresh_nikulin_coupled_charge_20261015_review','fresh_nikulin_coupled_charge_20261015_dimension166'):
        paths += [f'scripts/{stem}.py',f'results/{stem}.json',f'results/{stem}.md']
    for p in paths:
        if p.endswith('.json'):
            for z,h in json.loads((ROOT/p).read_text()).get('sources',{}).items():assert sha(z)==h,(p,z)
    return dict(status='PASS_INDEPENDENT_EXACT_HEIGHT_CHARGE_DEFECT_BUDGET_AND_165_REVIEW',
        sources={p:sha(p)for p in paths},arithmetic='Standalone Fraction grounded matrix inverses, exact KKT verification and polynomial coefficient arithmetic; no project/producer/SymPy/optimizer imports.',
        height_ratio_limit=str(B),exceptional_cases=ex,possible_defects=defects,complete_cycle_support_cases=cy,
        all_rank_BCD=bc,all_rank_A_tails=tail,counts=dict(exceptional_records=40,singletons=71,multiple_supports=16,possible_defects=15,cycle_duals=580,BCD_classes=24),
        aggregate=dict(total_cost=4,minimum_defective_factor_cost='3/2',maximum_defective_factors=2,maximum_each_defect='2/9',
            combined_cost_and_defect_loss=str(loss),near_at_N165=str(near),far_at_N165=str(far),near_slope='1/9',far_slope='31/270',
            ideal_branch_dimension_bound=165,compact_branch_inherited_bound=143,combined_candidate_bound=165),
        scope='Independent local support completeness, exact certificates, classical tails and common-height accounting review. Published narrowness, actual framing/endpoint existence, common ratio, Schur budget and compact143 branch remain explicitly inherited from the frozen independently reviewed geometric chain. Internal research proof review; no publication-priority or formal-verification claim.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
