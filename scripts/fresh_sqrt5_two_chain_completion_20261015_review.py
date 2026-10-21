#!/usr/bin/env python3
"""Independent exact two-chain review; no producer arithmetic imports.

Uses the frozen independent tuple-Fraction radical Gram implementation,
full congruence inertia, and a new exact positive-integer equation solver.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,combinations_with_replacement,product
from math import gcd,isqrt,lcm
from pathlib import Path
import fresh_sqrt5_quad_filters_20261015_review as R
import fresh_sqrt5_quad_contacts_20261015_review as Q

ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_sqrt5_two_chain_completion_20261015';TAG=BASE+'_review'
K=Q.K
load=lambda p:json.loads((ROOT/'results'/(p+'.json')).read_text())
dot=lambda a,b:sum((K(*x)*K(*y)for x,y in zip(a,b)),K())

def two_paths(core,contacts,lengths):
    n=core['nodes'];N=n+sum(lengths);g=Q.edges(core);s=n
    for contact,length in zip(contacts,lengths):
        assert length>0
        for i,m in enumerate(contact):g[i,s]=m
        for i in range(s,s+length-1):g[i,i+1]=3
        s+=length
    return dict(nodes=N,white=list(range(4,N)),labels=[g.get(p,2)for p in combinations(range(N),2)])

def gram(row):
    colors=R.coloring(row);assert colors is not None
    return R.matrix(row,colors)

def pd_corners(row,M):
    counts=[]
    for a,b in product((0,2),(1,3)):
        ids=row['white']+[a,b];G=[[M[i][j]for j in ids]for i in ids]
        counts.append(R.inertia(G))
    return counts

def exact_roots(p):
    a,b,c=p
    if not a:
        values=[]if not b else[-c/b]
    else:
        d=b*b-4*a*c
        if d<0:return []
        s=R.rroot(d)
        if s is None:return []
        values=[(-b-s)/(2*a),(-b+s)/(2*a)]
    return sorted({int(x)for x in values if x.denominator==1 and x>0})

def integer_solutions(coeff):
    """Solve both field coordinates, or factor the rational bilinear form.

    The actual catalogue uses no unbounded affine-parameter branch. The
    generic producer's linear branch is source-reviewed, not called here.
    """
    a,b,c=coeff;points=set();families=[]
    if any(x.b for x in coeff):
        elim=(a.a*b.b-a.b*b.a,c.a*b.b-c.b*b.a-a.b,-c.b)
        if any(elim):
            for x in exact_roots(elim):
                den=a*x+c;num=b*x+1
                if den:
                    y=-num/den
                    if not y.b and y.a.denominator==1 and y.a>0:points.add((x,int(y.a)))
                elif not num:families.append(('x',x))
        else:
            assert not c.b and b.b and a==b*c.a
            if c.a:
                y=-1/c.a
                if y.denominator==1 and y>0:families.append(('y',int(y)))
    else:
        den=lcm(*(x.denominator for x in (a.a,b.a,c.a,F(1))))
        A,B,C,D=[int(x*den)for x in(a.a,b.a,c.a,F(1))]
        g=gcd(gcd(abs(A),abs(B)),gcd(abs(C),abs(D)))
        A,B,C,D=[x//g for x in(A,B,C,D)]
        if A:
            H=B*C-A*D
            if H:
                divisors=set()
                for d in range(1,isqrt(abs(H))+1):
                    if H%d==0:divisors.update((d,-d,abs(H)//d,-abs(H)//d))
                for u in divisors:
                    v=H//u
                    if (u-C)%A==0 and (v-B)%A==0:
                        x,y=(u-C)//A,(v-B)//A
                        if min(x,y)>0:points.add((x,y))
            else:
                for tag,z in(('x',C),('y',B)):
                    if (-z)%A==0 and -z//A>0:families.append((tag,-z//A))
        elif not B or not C:
            if B:
                if (-D)%B==0 and -D//B>0:families.append(('x',-D//B))
            elif C:
                if (-D)%C==0 and -D//C>0:families.append(('y',-D//C))
            else:assert D
        else:
            # All actual linear cases have coefficients of the same sign;
            # positivity gives this exact finite bound, not a search cutoff.
            assert B*C>0
            if B<0:B,C,D=-B,-C,-D
            for x in range(1,max(0,(-D-C)//B)+1):
                if (-D-B*x)%C==0:
                    y=(-D-B*x)//C
                    if y>0:points.add((x,y))
    assert all(not(a*x*y+b*x+c*y+1)for x,y in points)
    for tag,z in families:
        assert (not(a*z+c)and not(b*z+1))if tag=='x'else(not(a*z+b)and not(c*z+1))
    return sorted(points),sorted(set(families))

def run():
    data=load(BASE);old=load('fresh_sqrt5_quad_contacts_20261015');unions=load('fresh_sqrt5_quad_exhaustion_20261015')['labelled_union_supports']
    sources=[r for r in old['records']if r['status']=='UNPAID_FINITE_CONTACT_DOMAIN'and not r['singular']and r['accepted']]
    expected={}
    for r in sources:
        ps=[(i,t['contacts'])for i,t in enumerate(r['trials'])if t['status']=='PASS_EXTRA_WHITE_NECESSARY_TESTS']
        assert [p for _,p in ps]==r['accepted']
        for a,b in combinations_with_replacement(ps,2):expected[r['class_index'],a[0],b[0]]=(r,[a[1],b[1]])
    assert len(sources)==130 and sum(len(r['accepted'])for r in sources)==458 and len(expected)==1455
    assert len(data['records'])==len(expected)
    assert {(r['source_class_index'],*r['trial_indices'])for r in data['records']}==set(expected)
    records=[];counts=Counter();totals=Counter();cache={};families=[];length_counts=Counter()
    for index,saved in enumerate(data['records']):
        key=(saved['source_class_index'],*saved['trial_indices']);src,contacts=expected[key]
        assert contacts==saved['contacts']and src['union_index']==saved['union_index']
        core=unions[src['union_index']];n=core['nodes'];head=two_paths(core,contacts,(1,1));N=head['nodes'];g=Q.edges(head)
        finite=all(Q.finite(head['white']+[a,b],g)for a,b in product((0,2),(1,3)))
        out=dict(key=list(key))
        if not finite:status='HEAD_REJECT_NONFINITE_CORNER'
        else:
            assert Q.no_short(N,g)and not Q.paid(head)
            M=gram(head);physical=R.inertia(M);other=R.inertia([[R.conj(x)for x in row]for row in M])
            for ii in pd_corners(head,M):assert ii==[N-2,0,0]
            totals['head_corners']+=4
            out.update(head_physical=physical,head_conjugate=other)
            if physical[1]>=2:status='HEAD_REJECT_TWO_NEGATIVE'
            elif other[1]:status='HEAD_REJECT_CONJUGATE_PSD'
            else:
                H=[row[:n]for row in M[:n]];J=R.inverse(H)
                assert R.inertia(H)==[n-1,1,0]and R.inertia([[R.conj(x)for x in row]for row in H])==[n,0,0]
                for i,j in product(range(n),repeat=2):assert dot(H[i],[J[k][j]for k in range(n)])==K(i==j)
                v=[row[n]for row in M[:n]];w=[row[n+1]for row in M[:n]]
                hv=[dot(row,v)for row in J];hw=[dot(row,w)for row in J]
                c11=sum((K(*x)*y for x,y in zip(v,hv)),K());c22=sum((K(*x)*y for x,y in zip(w,hw)),K());c12=sum((K(*x)*y for x,y in zip(v,hw)),K())
                d,e=K(*M[n][n]),K(*M[n+1][n+1]);aa=2*c11/d;bb=2*c22/e;ff=4*(c11*c22-c12*c12)/(d*e)
                co=[1-aa-bb+ff,1-aa,1-bb]
                assert co==[K(*R.unpack(x))for x in saved['polynomial_coefficients']]
                assert [d,e,c11,c22,c12]==[K(*R.unpack(x))for x in saved['schur_data']['first_norms']+saved['schur_data']['costs']]
                points,fams=integer_solutions(co);sol=saved['length_solutions']
                assert points==[tuple(x)for x in sol['finite_pairs']]
                theirs=sorted(('x',x['first'])if x['kind']=='FIXED_FIRST'else('y',x['second'])for x in sol['families'])
                assert fams==theirs
                totals['equations']+=1;totals['isolated_solutions']+=len(points)
                final=set(points);bounds=saved['family_corner_bounds'];assert len(bounds)==len(fams)
                for bound in bounds:
                    bad=two_paths(core,contacts,bound['failed_lengths']);B=gram(bad);ids=bound['corner_nodes'];BG=[[B[i][j]for j in ids]for i in ids]
                    assert BG==[[R.unpack(x)for x in row]for row in bound['corner_gram']]
                    ii=R.inertia(BG);assert ii==bound['corner_inertia']and ii!=[len(ids),0,0]
                    assert set(ids)==set(bad['white']+bound['marked_corner'])
                    for vker in bound.get('corner_kernel',[]):
                        z=[R.unpack(x)for x in vker];assert all(dot(row,z)==0 for row in BG)
                    assert bound['status']=='ALL_FREE_LENGTHS_AT_LEAST_TWO_REJECTED'
                    fam=bound['family'];assert fam['kind']in('FIXED_FIRST','FIXED_SECOND')
                    assert (fam.get('first')or fam.get('second'))==1
                    remaining=bound['remaining_positive_pairs'];assert remaining==[[1,1]];final.add((1,1))
                    families.append(dict(key=list(key),failed_lengths=bound['failed_lengths'],corner_inertia=ii))
                stars=saved['completed_stars'];assert {tuple(s['lengths'])for s in stars}==final and len(stars)==len(final)
                for star in stars:
                    ell,emm=star['lengths'];row=two_paths(core,contacts,(ell,emm));assert row==star['row']
                    G=gram(row);k=len(G);phy=R.inertia(G);con=R.inertia([[R.conj(x)for x in rr]for rr in G])
                    assert phy==star['physical_inertia']==[k-2,1,1]and con==star['conjugate_inertia']==[k-1,0,1]
                    for ii in pd_corners(row,G):assert ii==[k-2,0,0]
                    assert Q.no_short(k,Q.edges(row))and not Q.paid(row)
                    # Direct inverses of both literal path blocks independently
                    # check the all-length determinant-lemma normalization.
                    for start,length in((n,ell),(n+ell,emm)):
                        ids=list(range(start,start+length));P=[[G[i][j]for j in ids]for i in ids];PI=R.inverse(P)
                        expected_first=K(2*F(length,length+1))/K(*G[start][start])
                        assert K(*PI[0][0])==expected_first
                        totals['path_inverse_checks']+=1
                    totals['completed_stars']+=1;totals['completed_corners']+=4;length_counts[f'{ell},{emm}']+=1
                status='FINITE_UNPAID_RANK_COMPLETIONS'if stars else'NO_POSITIVE_INTEGER_LENGTH_PAIR'
                out.update(coefficients=[x.enc()for x in co],finite_pairs=[list(p)for p in points],families=[list(x)for x in fams],completed_lengths=[s['lengths']for s in stars])
        assert status==saved['status'],(key,status,saved['status'])
        counts[status]+=1;out['status']=status;records.append(out)
        if (index+1)%200==0:print('reviewed',index+1,'equations',totals['equations'],flush=True)
    assert dict(counts)==data['counts']and dict(length_counts)==data['surviving_length_pairs']
    assert totals['equations']==239 and totals['isolated_solutions']==107 and len(families)==5 and totals['completed_stars']==112
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',
        'scripts/fresh_sqrt5_quad_filters_20261015_review.py','results/fresh_sqrt5_quad_filters_20261015_review.json',
        'scripts/fresh_sqrt5_quad_contacts_20261015_review.py','results/fresh_sqrt5_quad_contacts_20261015_review.json',
        'results/fresh_sqrt5_quad_contacts_20261015.json','results/fresh_sqrt5_quad_exhaustion_20261015.json',
        'results/fresh_sqrt5_single_chain_completion_20261015_review.json']
    manifests=Q.manifest_check(paths)
    return dict(status='PASS_INDEPENDENT_ALL_LENGTH_TWO_CHAIN_REVIEW',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},recursive_manifests=manifests,profile_pairs=1455,counts=dict(counts),totals=dict(totals),surviving_length_pairs=dict(length_counts),infinite_family_corner_failures=families,records=records,
      scope='All1455 specified profile pairs, all239 exact rank equations, all107 isolated solutions and5 corner-bounded free-length families, all112 principal stars. Independent Fraction Gram/inertia and length arithmetic; no producer arithmetic imports. The prior458-profile completeness and compact/arithmetic premises are inherited. No independent exhaustive no-affine test, geometric quadrilateral realization, arbitrary branching/multiple touching nodes, dotted-core or global dimension claim.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();out=run();path=ROOT/'results'/(TAG+'.json')
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],out['totals'])
