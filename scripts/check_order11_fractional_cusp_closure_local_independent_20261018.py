#!/usr/bin/env python3
"""Independent physical-corner checks and fresh order11 ideal convolutions."""
import fcntl
from fractions import Fraction as F
from math import prod
import signal
import sympy as s
import check_a16_rank17_spectator_extension_independent_20261018 as I

ROOT=I.ROOT;W=I.W;OLD=I.OLD;TAG='fresh_order11_fractional_cusp_closure_local_20261018'
DELTA=F(11,512)


def matrix(rows):return s.Matrix([[s.sympify(x)for x in row]for row in rows])


def inertia(g):return I.G.inertia([[F(x)for x in row]for row in g.tolist()])


def adj(g,ids):return{i:{j for j in ids if i!=j and g[i,j]}for i in ids}


def components(g,ids):
    graph=adj(g,ids);todo=set(ids);out=[]
    while todo:
        block={min(todo)};queue=list(block)
        for i in queue:
            for j in sorted(graph[i]-block):block.add(j);queue.append(j)
        todo-=block;out.append(sorted(block))
    return out


def determinant(R):
    h,k,l,m,c,d,p,q,r=s.symbols('h k l m c d p q r',nonzero=True);a=s.Rational(1,R)
    g=s.Matrix([[a,-a,-h+l*a,-m*a],[-a,a,-l*a,-k+m*a],
        [-h+l*a,-l*a,c+l*l*a-p,-l*m*a-q],[-m*a,-k+m*a,-l*m*a-q,d+m*m*a-r]])
    expected=1-(2*l/h+c/h**2+2*m/k+d/k**2)/R+(p/h**2-2*q/(h*k)+r/k**2)/R
    assert s.factor(g.det()/(h*h*k*k)-expected)==0
    return expected


def profiles(doc):
    b=doc['B_center'];g=matrix(b['gram']);nu=s.Matrix(list(map(s.sympify,b['kernel'])))
    assert g.rows==11 and g*nu==s.zeros(11,1)and g.rank()==10 and all(x>0 for x in nu)
    assert b['marks']==[0,8]and b['resistance']==8
    white=b['white'];J=g.extract(white,white).inv();va,vb=[-nu[i]*g.extract([i],white).T for i in(0,8)]
    M=s.diag(nu[0],nu[8])*g.extract([0,8],[0,8])*s.diag(nu[0],nu[8])-s.Matrix.hstack(va,vb).T*J*s.Matrix.hstack(va,vb)
    assert M.applyfunc(s.simplify)==matrix(b['old_mark_schur'])==s.Matrix([[1,-1],[-1,1]])/8
    for row in b['finite_inverse_diagonals']:
        iv=g.extract(row['nodes'],row['nodes']).inv();di=dict(zip(row['nodes'],iv.diagonal()))
        assert{str(i):str(x)for i,x in di.items()}==row['diagonal']
        assert di[row['endpoint']]==1 and di[row['next_node']]==2
        assert all(x>2 for i,x in di.items()if i not in(row['endpoint'],row['next_node']))
    t,u,v,w,z=s.symbols('t u v w z',nonnegative=True);symbols={str(x):x for x in(t,u,v,w,z)}
    parsed={};checks=0
    for row in b['profile_identities']:
        contact={int(i):s.sympify(x,locals=symbols)for i,x in row['contacts'].items()}
        h,l,c=[s.sympify(row[k],locals=symbols)for k in('height','moment','c0')]
        q=s.Matrix([contact.get(i,0)for i in white]);far=vb if row['end']=='center'else va
        assert s.simplify(sum(nu[i]*x for i,x in contact.items())-h)==0
        assert s.simplify(8*(far.T*J*q)[0]-l)==0
        assert s.simplify(2-(q.T*J*q)[0]-c-l*l/8)==0
        assert s.simplify(s.sympify(row['cost'],locals=symbols)-(2*l/h+c/h**2))==0
        parsed[row['end'],row['profile']]=(contact,h,l,c,q);checks+=1
    assert len(parsed)==8
    for row in b['cross_identities']:
        p=parsed['center',row['center']];q=parsed['terminal',row['terminal']]
        assert s.simplify((p[4].T*J*q[4])[0]-p[2]*q[2]/8)==0
        assert s.simplify(s.sympify(row['product'],locals=symbols)-(p[4].T*J*q[4])[0])==0;checks+=1
    assert len(b['cross_identities'])==16
    for row in b['universal_profile_bounds']:
        _,h,l,c,_=parsed[row['end'],row['profile']];pol=row['maximum']*h*h-2*l*h-c
        if row['end']=='terminal':pol=pol.subs(t,w+(1 if row['profile']=='isolated'else 0))
        elif row['profile']=='outer':
            branch=row['branch'];sub={'opposite_nonzero':{t:1+w},'outer_native_gap_left':{u:s.sqrt(2)+w},
                'outer_native_gap_right':{v:s.sqrt(2)+w},'two_outer_contacts':{u:1+w,v:1+z}}[branch]
            pol=pol.subs(sub)
        num=s.sympify(row['polynomial']['numerator'],locals=symbols);den=s.sympify(row['polynomial']['denominator'],locals=symbols)
        assert s.simplify(pol-num/den)==0
        variables=sorted(num.free_symbols|den.free_symbols,key=str)
        assert all(x.is_nonnegative for x in s.Poly(num,*variables).coeffs())
        assert all(x.is_nonnegative for x in s.Poly(den,*variables).coeffs())and den.subs({x:0 for x in variables})>0
        checks+=1
    assert len(b['universal_profile_bounds'])==11
    # Compare rational functions using the exact nonzero symbols of the
    # determinant construction, avoiding an unproved sampled identity.
    for R,stored in((8,b['determinant']),(10,doc['D_next']['determinant'])):
        sy={str(x):x for x in s.symbols('h k l m c d p q r',nonzero=True)}
        assert s.simplify(s.sympify(stored,locals=sy)-determinant(R))==0
    return checks


def corners(doc):
    count=0
    for row in doc['B_corner_templates']:
        g=matrix(row['gram']);ids=row['corner_nodes'];marks=row['marks'];gauge=row['gauge']
        if row['branch']=='nonextreme':
            ns=row['whole_affine_nodes'];a=g.extract(ns,ns);z=s.Matrix(list(map(s.sympify,row['positive_kernel'])))
            assert a*z==s.zeros(len(ns),1)and all(x>0 for x in z)and inertia(a)==[8,0,1]
            graph=adj(g,ns);assert len(ns)==9 and len(graph[marks[0]])==len(graph[marks[1]])==1
            assert sum(len(graph[i])==3 for i in ns)==1
            assert F(row['payment'])==1
        else:
            assert inertia(g.extract(ids,ids))==[len(ids),0,0]
            blocks=components(g,ids);a,b=[next(c for c in blocks if m in c)for m in marks]
            assert a is not b and len(b)==8;graph=adj(g,ids)
            assert len(graph[marks[1]])==1
            j=next(iter(graph[marks[1]]));assert g[marks[1],j]**2==2*gauge[marks[1]]*gauge[j]
            if row['branch']=='extreme_empty':assert len(a)==2 and len(graph[marks[0]])==1
            else:assert row['branch']=='extreme_one_A1'and len(a)==3 and len(graph[marks[0]])==2
            assert F(row['new_coupon'])==4*DELTA and F(row['selected_ideal'])==1-4*DELTA and F(row['total'])==1
        count+=1
    for row in doc['D_next']['templates']:
        g=matrix(row['gram']);white=row['white'];n=g.rows
        assert inertia(g)==[n-2,1,1]and inertia(g.extract(white,white))==[n-4,0,0]
        for root in(7,9,11,12):assert inertia(g.extract(sorted(white+[root]),sorted(white+[root])))==[n-3,0,0]
        for c in row['corners']:
            ns=sorted(white+c['marks']);assert inertia(g.extract(ns,ns))==c['signature']and c['signature'][1]==0
        ids=sorted(white+row['paying_marks']);assert inertia(g.extract(ids,ids))==[n-2,0,0]
        blocks=components(g,ids);a,b=[next(c for c in blocks if m in c)for m in row['paying_marks']]
        if row['extreme_attachment']==7:
            assert a is not b and len(a)==(2 if row['whole_A1_spectator']else 1)and len(b)==10
            assert F(row['payment'])==F(1,2)+DELTA
        else:assert a is b and len(a)==(12 if row['whole_A1_spectator']else 11)and F(row['payment'])==1
        tn=row['next_opposite_product'];assert tn==(0 if row['whole_A1_spectator']else 2)
        assert F(row['next_height'])==1+F(tn,2)and F(row['next_cost'])==4/(1+F(tn,2))
        assert F(row['total'])==F(row['payment'])+F(1,2)-DELTA>=1;count+=1
    t=s.symbols('t',nonnegative=True);dn=doc['D_next']
    assert s.simplify(s.sympify(dn['next_cost'],locals={'t':t})-4/(1+t/2))==0
    assert dn['both_next_cost_sum_upper']==8 and F(dn['both_next_determinant_gap'])==F(1,5)
    for row in doc['B_zero_repair']:
        g=matrix(row['gram']);white=row['white'];n=g.rows
        assert inertia(g.extract(white,white))==[n-4,0,0]
        for root in(8,9,11,12):assert inertia(g.extract(sorted(white+[root]),sorted(white+[root])))==[n-3,0,0]
        assert row['first_corner']!=row['second_corner']
        first=F(1)if row['fork_attachment']==10 else F(1,2)+DELTA
        second=F(1,2)+DELTA if row['heavy_profile']=='isolated'else F(1,2)-DELTA
        assert first==F(row['first_payment'])and second==F(row['second_payment'])and first+second==F(row['total'])>=1;count+=1
    assert count==13
    return count


def profiles_ideal(doc):
    X=I.I.I.X
    env={int(k):F(v)for k,v in X.object_at(ROOT/'results/fresh_order10_low8_ambient58_frontier_20261018.json',b'"current_envelopes": {').items()}
    env[9]=F(22,5);env[15]=F(51,8);env[16]=F(895,136);assert env[10]==F(53,11)
    old=dict(env);env[10]-=F(1,128)
    domain=sorted((lo,)+p if lo else p for lo in range(8)for p in X.parts(56-lo));ii=doc['ideal']
    assert len(domain)==2270==ii['complete_profiles']and domain==[tuple(r['ranks'])for r in ii['rows']]
    checks=0;margins=[]
    for rs,row in zip(domain,ii['rows']):
        h=[1]
        for r in rs:
            out=[0]*(len(h)+r)
            for i,a in enumerate(h):
                for j in range(r+1):out[i+j]+=a;checks+=1
            h=out
        v=prod(r+1 for r in rs);assert len(h)==57 and h==h[::-1]and sum(h)==v
        g=sum(min(i,56-i)*a for i,a in enumerate(h));curv=F(5*g-56*v,2*v)
        m=curv-sum(env[r]for r in rs);before=curv-sum(old[r]for r in rs)
        assert m==F(row['normalized_margin'])>0 and before==F(row['old_margin'])
        assert m-before==F(rs.count(10),128)==F(row['normalized_saving']);margins.append(m)
    assert min(margins)==F(ii['normalized_minimum'])and ii['new_negative_count']==ii['zero_count']==0
    totals={k:F(v)for k,v in doc['order11_totals'].items()}
    assert totals==dict(A=F(99,2)+11*DELTA,B=53-4*DELTA,C=F(52),D=53-4*DELTA)
    assert max(totals.values())/11==env[10]
    return checks,min(margins)


def main():
    d=W.load(TAG);hashes=0
    for p,h in d['sources'].items():assert W.digest(ROOT/p)==h;hashes+=1
    assert F(d['delta'])==DELTA; a=profiles(d);b=corners(d);c=profiles_ideal(d);OLD.U.memory_guard()
    print('PASS_INDEPENDENT_ORDER11_FRACTIONAL_CUSP_CLOSURE_LOCAL',hashes,a,b,*c,flush=True)


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);OLD.R.resources();main();OLD.R.resources(start=False)
