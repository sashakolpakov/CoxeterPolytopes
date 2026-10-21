#!/usr/bin/env python3
"""Complete bounded Vinberg enumeration using both primitive norm classes."""
import argparse
import json
import sys
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from coxeter_polytopes.exact.golden import Golden as K, PHI, ZERO
from coxeter_polytopes.search.bugaenko_vinberg import (
    _fraction_inverse,_prefix_schur_complements,_integer_upper_square_root,_mat_vec,_pair)
from coxeter_polytopes.search.h11_vinberg import (
    _coordinate_alphabet,_integral_form,_ring_integer,_from_ring_integer,
    _ri_add,_ri_subtract,_ri_scale,_ri_multiply,_ri_sign,h11_basis_gram)

TAG='fresh_h9_beta11_vinberg_20261014'
INPUTS=['results/fresh_h9_h8_cap_height_20261014_host_basis.json',
        'results/fresh_h9_h8_cap_height_20261014_norm_pool.json',
        'results/fresh_h9_h8_cap_height_20261014_norm_pool_review.json',
        'results/fresh_h9_beta11_zero_sector_20261014.json']
DELTA=1+4*PHI; BETA=3+PHI
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def old(x):
    assert all(F(t)==0 for t in x[1])
    a,b=map(F,x[0]);return K(a-b,2*b)
def enc(x):
    if isinstance(x,K):return x.to_json()
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [enc(v)for v in x]
    return x
def read(p):
    d=json.loads((ROOT/p).read_text())
    for f,h in d['sources'].items():assert sha(f)==h,f
    return d
def controls(q,T):
    U=_integer_upper_square_root(2*q*T/DELTA,0)
    V=_integer_upper_square_root(-2*q/DELTA,1)
    # |b|=|x-x'|/sqrt5<(U+V)/2, |a|<V+|b|.
    bb=(U+V)//2+1;aa=V+bb
    xs=[]
    for a in range(-aa,aa+1):
        for b in range(-bb,bb+1):
            x=K(a,b)
            if x.sign()>0 and (T-DELTA*x*x/(2*q)).sign()>=0 and (2*q+DELTA*x*x).sign(1)>=0:
                xs.append(x)
    return sorted(xs),dict(identity_upper=U,conjugate_upper=V,a_bound=aa,b_bound=bb)
def beta_divisible(pair):return (pair[0]+8*pair[1])%11==0
def relation(g,qa,qb):
    assert g.sign()<=0
    if not g:return 'right'
    a=4*g*g/(qa*qb)
    if a==4:return 'parallel'
    if a>4:return 'dotted'
    for value,label in [(K(1),3),(K(2),4),(PHI*PHI,5),(K(3),6),(PHI+2,10)]:
        if a==value:return 'pi'+str(label)
    raise AssertionError(('unclassified finite pair',g,qa,qb,a))
def old_host_isometry(G):
    oldG=h11_basis_gram()
    C=[[G[i+1][j+1]for j in range(8)]for i in range(8)]
    D=[[oldG[i+2][j+2]for j in range(8)]for i in range(8)]
    mappings=[]
    def rec(p):
        i=len(p)
        if i==8:mappings.append(p);return
        degree=sum(bool(z)for j,z in enumerate(C[i])if j!=i)
        for k in range(8):
            if k in p or sum(bool(z)for j,z in enumerate(D[k])if j!=k)!=degree:continue
            if all(C[i][j]==D[k][p[j]]for j in range(i)):rec(p+[k])
    rec([]);assert len(mappings)==1;p=mappings[0]
    M=[[K(0)for _ in range(10)]for _ in range(10)];M[0][0]=K(1);M[1][9]=K(-1)
    for i,j in enumerate(p):M[j+2][i+1]=K(1)
    cols=list(zip(*M))
    assert all(_pair(cols[i],oldG,cols[j])==G[i][j]for i in range(10)for j in range(10))
    return dict(new_E8_node_to_old=p,old_coordinates_from_new_matrix=M,
                proof='M transpose times literal oldH11 Gram times M equals literal newhost Gram; old second coordinate is minus new t.')

def enumerate_level(x,q,previous,J,schur,G):
    target=2*q+DELTA*x*x
    assert target.sign()>=0 and target.sign(1)>=0
    alphabet=_coordinate_alphabet(target);values=tuple(_ring_integer(z)for z in alphabet.values)
    forms=tuple(_integral_form(z)for z in schur)
    xi,ti=_ring_integer(x),_ring_integer(target)
    visited=[0]*10;sp=[0]*10;wp=[0]*10
    parity=0;residue=0;mismatch=0;solutions=[]
    coeff=[];bounds=[]
    for r in previous:
        L=r[0]-2*r[9]
        assert L.sign()>=0 and all(z.sign()>=0 for z in r[1:9])
        coeff.append(tuple(_ring_integer(z)for z in (L,*(2*z for z in r[1:9]))))
        bounds.append(_ring_integer(DELTA*x*r[0]))
    def extend(prefix,minimum,partials):
        nonlocal parity,residue,mismatch
        depth=len(prefix);visited[depth]+=1
        if depth:
            num,den=minimum
            rem=_ri_subtract(_ri_scale(ti,den),num)
            if _ri_sign(rem)<0 or _ri_sign(rem,1)<0:sp[depth]+=1;return
            entry=prefix[-1]
            partials=tuple(_ri_add(partials[i],_ri_multiply(entry,coeff[i][depth-1]))for i in range(len(coeff)))
            if any(_ri_sign(_ri_subtract(b,v))<0 for b,v in zip(bounds,partials)):
                wp[depth]+=1;return
        if depth==9:
            num,den=minimum
            if num!=_ri_scale(ti,den):mismatch+=1;return
            l,*lam=map(_from_ring_integer,prefix)
            y=_mat_vec(J,lam);t=(x-l)/2
            assert t.is_integral and all(z.is_integral for z in y)
            r=(x,*y,t)
            assert _pair(r,G,r)==q
            assert all((2*_pair(r,G,tuple(K(i==j)for i in range(10)))/q).is_integral for j in range(10))
            solutions.append(r);return
        den,form=forms[depth]
        base=(0,0);linear=(0,0)
        for i in range(depth):
            linear=_ri_add(linear,_ri_scale(prefix[i],2*form[i][depth]))
            for j in range(depth):base=_ri_add(base,_ri_scale(_ri_multiply(prefix[i],prefix[j]),form[i][j]))
        for entry in values:
            if q==2*BETA and not beta_divisible(entry):residue+=1;continue
            if depth==0 and any((entry[i]-xi[i])%2 for i in range(2)):parity+=1;continue
            child=_ri_add(_ri_add(base,_ri_multiply(linear,entry)),_ri_scale(_ri_multiply(entry,entry),form[depth][depth]))
            extend((*prefix,entry),(child,den),partials)
    extend((),None,tuple((0,0)for _ in coeff))
    assert all(_pair(a,G,b).sign()<=0 for i,a in enumerate(solutions)for b in solutions[i+1:])
    return dict(control=x,norm=q,height=DELTA*x*x/(2*q),target=target,
        alphabet=dict(values=alphabet.values,identity_upper=alphabet.identity_upper,
                      conjugate_upper=alphabet.conjugate_upper,a_bound=alphabet.a_bound,b_bound=alphabet.b_bound),
        visited_by_depth=visited,schur_pruned_by_depth=sp,wall_pruned_by_depth=wp,
        parity_pruned=parity,beta_residue_pruned=residue,leaf_norm_mismatches=mismatch,
        roots=solutions)

def run(T):
    basis,pool,review,zero=map(read,INPUTS)
    G=tuple(tuple(old(x)for x in row)for row in basis['gram'])
    isometry=old_host_isometry(G)
    H=tuple(map(old,zero['basepoint']));cone=tuple(tuple(map(old,r))for r in zero['simple_zero_walls'])
    assert _pair(H,G,H)==-2*DELTA and len(cone)==9
    C=tuple(tuple(G[i+1][j+1].a for j in range(8))for i in range(8))
    assert all(not G[i+1][j+1].b for i in range(8)for j in range(8))
    J=_fraction_inverse(C)
    assert all(z.denominator==1 and z>=0 for row in J for z in row)
    form=((F(1),)+(F(0),)*8,)+tuple((F(0),)+tuple(2*z for z in row)for row in J)
    schur=_prefix_schur_complements(form)
    assert all(schur[i][i][i]>=1 for i in range(9))
    # Each coordinate bound holds at both embeddings: the diagonal of the
    # inverse positive form is1 for l and C_ii/2=1 for every weight.
    Fi=_fraction_inverse(form);assert all(Fi[i][i]==1 for i in range(9))
    controlrows=[];ordered=[]
    for q in (K(2),2*BETA):
        xs,box=controls(q,T);controlrows.append(dict(norm=q,box=box,values=xs))
        ordered.extend((DELTA*x*x/(2*q),q,x)for x in xs)
    ordered.sort(key=lambda z:z[0])
    assert all(a[0]!=b[0]for a,b in zip(ordered,ordered[1:]))
    roots=[];levels=[]
    for height,q,x in ordered:
        level=enumerate_level(x,q,roots,J,schur,G)
        roots.extend(level['roots']);levels.append(level)
        print('height',height,'q',q,'x',x,'roots',len(level['roots']),'nodes',sum(level['visited_by_depth']),flush=True)
    allroots=(*cone,*roots)
    assert len(set(allroots))==len(allroots)
    assert all(_pair(r,G,H).sign()<0 for r in roots)
    assert all(_pair(a,G,b).sign()<=0 for i,a in enumerate(allroots)for b in allroots[i+1:])
    gram=[[_pair(a,G,b)for b in allroots]for a in allroots]
    cartan=[[2*gram[i][j]/gram[j][j]for j in range(len(allroots))]for i in range(len(allroots))]
    assert all(z.is_integral for row in cartan for z in row)
    relations=[[relation(gram[i][j],gram[i][i],gram[j][j])if i!=j else 'diagonal'
                for j in range(len(allroots))]for i in range(len(allroots))]
    paircounts=Counter(relations[i][j]for i in range(len(allroots))for j in range(i))
    # Move strictly into the zero sector, with a verified dyadic displacement.
    zeroG=tuple(tuple(_pair(a,G,b).a for b in cone)for a in cone)
    zeroJ=_fraction_inverse(zeroG);weights=[sum(row)for row in zeroJ]
    Z=tuple(-sum((weights[i]*cone[i][j]for i in range(9)),ZERO)for j in range(10))
    assert all(_pair(Z,G,n)==-1 for n in cone) and not _pair(H,G,Z)
    epsilon=F(1)
    while True:
        center=tuple(a+epsilon*b for a,b in zip(H,Z))
        if _pair(center,G,center)<0 and all(_pair(center,G,n)<0 for n in allroots):break
        epsilon/=2
    assert _pair(center,G,H)<0
    facets=[]
    for i,n in enumerate(allroots):
        a=_pair(center,G,n)/gram[i][i]
        point=tuple(c-a*r for c,r in zip(center,n))
        assert _pair(point,G,point)<0 and _pair(point,G,H)<0
        products=[_pair(point,G,r)for r in allroots]
        assert not products[i] and all(z<0 for j,z in enumerate(products)if j!=i)
        facets.append(dict(index=i,point=point,norm=_pair(point,G,point)))
    paths=INPUTS+[f'scripts/{TAG}.py','src/coxeter_polytopes/exact/golden.py',
        'src/coxeter_polytopes/search/bugaenko_vinberg.py','src/coxeter_polytopes/search/h11_vinberg.py',
        'src/coxeter_polytopes/lattices.py']
    return enc(dict(status='PASS_COMPLETE_BOTH_NORM_VINBERG_PREFIX',sources={p:sha(p)for p in paths},
        height_bound=T,field='Q(phi),phi²=phi+1; dictionaries a,b mean a+bphi.',
        form=G,basepoint=H,zero_walls=cone,norm_classes=[K(2),2*BETA],
        control_boxes=controlrows,positive_form=form,prefix_schur_complements=schur,
        inverse_E8=J,levels=levels,roots=allroots,gram=gram,cartan=cartan,
        pair_relations=relations,pair_counts=dict(paircounts),strict_center=center,
        center_displacement=dict(direction=Z,epsilon=epsilon),actual_facet_points=facets,
        old_host_signed_permutation=isometry,
        counts=dict(zero_walls=len(cone),positive_levels=len(levels),positive_walls=len(roots),
                    total_walls=len(allroots),nodes=sum(sum(z['visited_by_depth'])for z in levels)),
        duplicate_host='The legacy h11_vinberg module is the same discriminant11 host after a binary sign change and E8 relabeling. Its old norm2-only accepted wall list is not used as a full-group prefix.',
        completeness='All primitive integral reflection directions have unit-normalized norm2 or2beta. Exact control-coordinate boxes exhaust both classes through the stated geometric height; their levels are merged by exact height. The zero sector is independently certified. Every positive-definite coordinate box is exhaustive at both embeddings. Schur and earlier-wall pruning use exact integer-pair arithmetic only.',
        scope='Complete bounded prefix for the full integral reflection group of this lattice, not a claim of termination, reflectivity or finite volume. Geometry/cusp/vertex closure is separate. No previous15-wall chamber or norm2-only prefix is presumed complete.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--height')
    args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json'
    height=args.height if args.height is not None else (json.loads(path.read_text())['height_bound'] if args.verify else '10')
    T=F(height);assert T>0;d=run(T)
    if args.verify:assert d==json.loads(path.read_text())
    else:path.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'])
