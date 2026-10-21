#!/usr/bin/env python3
"""Exact finite-vertex/edge queue and complete spherical endpoint profiles."""
import argparse,json,sys
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_makarov_face_extension_20261013_independent as A
R=A.R
OUT='results/fresh_h9_makarov_frontier_20261013.json'
def inv(x):
    assert set(x)<=set((1,5,11,55))
    conj={n:(-a if n in (11,55) else a)for n,a in x.items()}
    norm=R.mul(x,conj);return R.mul(conj,A.invk(norm))
def div(x,y):return R.mul(x,inv(y))
def matrix_inverse(M):
    n=len(M);B=[list(row)+[R.C(int(i==j))for j in range(n)]for i,row in enumerate(M)]
    for j in range(n):
        k=next(i for i in range(j,n)if B[i][j]);B[k],B[j]=B[j],B[k]
        z=inv(B[j][j]);B[j]=[R.mul(a,z)for a in B[j]]
        for i in range(n):
            if i==j:continue
            z=B[i][j];B[i]=[A.sub(a,R.mul(z,b))for a,b in zip(B[i],B[j],strict=True)]
    return [row[n:]for row in B]
def combine(co,roots):return [R.psum(R.mul(a,r[j])for a,r in zip(co,roots,strict=True))for j in range(len(roots[0]))]
def pd_inverse(M):
    B=[[dict(z)for z in row]for row in M]
    for j in range(len(B)):
        d=B[j][j]
        if R.sign(d)<=0:return None
        for i in range(j+1,len(B)):
            for k in range(j+1,len(B)):B[i][k]=A.sub(B[i][k],div(R.mul(B[i][j],B[j][k]),d))
    return matrix_inverse(M)
def encvec(v):return list(map(R.enc,v))
def null_ray(M):
    B=[[dict(z)for z in row]for row in M];rank=0;piv=[]
    for j in range(10):
        k=next((i for i in range(rank,len(B))if B[i][j]),None)
        if k is None:continue
        B[k],B[rank]=B[rank],B[k];z=inv(B[rank][j]);B[rank]=[R.mul(a,z)for a in B[rank]]
        for i in range(len(B)):
            if i==rank:continue
            z=B[i][j];B[i]=[A.sub(a,R.mul(z,b))for a,b in zip(B[i],B[rank],strict=True)]
        piv.append(j);rank+=1
        if rank==len(B):break
    if rank!=9:return None
    free=next(j for j in range(10)if j not in piv);v=[{}for _ in range(10)];v[free]=R.C(1)
    for i,j in enumerate(piv):v[j]=R.scale(B[i][free],-1)
    assert all(not R.pair(v,row)for row in M);return v
def components(M):
    left=set(range(len(M)));out=[]
    while left:
        seen={min(left)};todo=list(seen)
        while todo:
            i=todo.pop()
            for j in sorted(left-seen):
                if M[i][j]:seen.add(j);todo.append(j)
        left-=seen;out.append(sorted(seen))
    return out
def spherical_profiles(M,ids):
    # The polar wall (index 6) is isolated in every finite support, and can be
    # normalized to norm two without changing any off-diagonal entry.
    N=[[dict(z)for z in row]for row in M]
    for i in range(len(N)):
        if N[i][i]!=R.C(2):
            assert ids[i]==6 and all(not N[i][j]for j in range(len(N))if i!=j);N[i][i]=R.C(2)
    I=pd_inverse(N);assert I is not None;cc=components(N);which={i:j for j,c in enumerate(cc)for i in c}
    rows=[{'contacts':[],'projection_square':R.enc({}),'residual_norm':R.enc(R.C(2))}];tails=[]
    squares={3:R.C(1),4:R.C(2),5:{1:F(3,2),5:F(1,2)}}
    for i in range(len(N)):
        if len(cc[which[i]])==1:
            assert I[i][i]==R.C(F(1,2));tails.append({'wall':ids[i],'labels':'every integer m >= 3','residual':'2 sin(pi/m)^2'});continue
        for m,sq in squares.items():
            p=R.mul(sq,I[i][i]);q=A.sub(R.C(2),p)
            if R.sign(q)>0:rows.append({'contacts':[[ids[i],m]],'projection_square':R.enc(p),'residual_norm':R.enc(q)})
    for n in (2,3):
        for inds in combinations(range(len(N)),n):
            # A spherical diagram is a forest: two contacts into the same
            # existing connected component would make a cycle.
            if len({which[i]for i in inds})!=n:continue
            for labels in (product(squares,repeat=n)if n==2 else [(3,3,3)]):
                p=R.psum(R.mul(squares[m],I[i][i])for i,m in zip(inds,labels,strict=True));q=A.sub(R.C(2),p)
                if R.sign(q)>0:rows.append({'contacts':[[ids[i],m]for i,m in zip(inds,labels,strict=True)],'projection_square':R.enc(p),'residual_norm':R.enc(q)})
    return rows,tails
def run(verify=False):
    d=R.read(A.SOURCE);Q=[[A.k(z)for z in row]for row in d['form']]
    roots=[[A.k(z)for z in row]for row in d['seed_normals']]+[[A.l(z)for z in d['changed_center_cap']['normal']]]
    h=list(map(A.k,d['strict_center']));gram=[[A.pair(a,Q,b)for b in roots]for a in roots]
    norms=[gram[i][i]for i in range(12)]
    dotted={(i,j)for i in range(12)for j in range(i)if R.sign(A.sub(R.mul(gram[i][j],gram[i][j]),R.mul(norms[i],norms[j])))>0}
    assert dotted=={(6,5),(11,4),(11,6)}
    stat=Counter();vertices=[];byactive={};vdata=[]
    for ss in combinations(range(12),9):
        stat['supports']+=1
        if any((max(i,j),min(i,j))in dotted for i,j in combinations(ss,2)):stat['dotted_pair_excluded']+=1;continue
        M=[[gram[i][j]for j in ss]for i in ss];I=pd_inverse(M)
        if I is None:stat['not_positive_definite']+=1;continue
        stat['positive_definite']+=1;co=A.mv(I,[A.pair(h,Q,roots[i])for i in ss]);v=R.vadd(h,R.vs(combine(co,[roots[i]for i in ss]),-1))
        q=A.pair(v,Q,v);assert R.sign(q)<0 and R.sign(A.pair(v,Q,h))<0
        ps=[A.pair(v,Q,r)for r in roots]
        if any(R.sign(z)>0 for z in ps):stat['infeasible']+=1;continue
        active=[i for i,z in enumerate(ps)if not z];assert active==list(ss)
        vi=len(vertices);byactive[tuple(active)]=vi;vdata.append((v,q,ss,I))
        vertices.append({'active':active,'point':encvec(v),'norm':R.enc(q),'all_wall_pairings':encvec(ps)})
    edges=[];profile_dict={};targetstats=Counter()
    for vi,(v,q,ss,I)in enumerate(vdata):
        for ji,released in enumerate(ss):
            support=[i for i in ss if i!=released];dual=combine(I[ji],[roots[i]for i in ss]);dualq=A.pair(dual,Q,dual)
            assert dualq==I[ji][ji]and R.sign(dualq)>0 and not A.pair(v,Q,dual)
            pp=[A.pair(v,Q,w)for w in roots];bb=[A.pair(dual,Q,w)for w in roots]
            assert all(bb[i]==R.C(int(i==released))for i in ss)
            collisions=[]
            for i in range(12):
                if i in ss:continue
                assert R.sign(pp[i])<0
                if R.sign(bb[i])<0:collisions.append((i,div(pp[i],bb[i])))
            first=None
            for i,t in collisions:
                if first is None or R.sign(A.sub(t,first[1]))<0:first=(i,t)
            row={'vertex':vi,'released':released,'support':support,'dual':encvec(dual),'dual_norm':R.enc(dualq),'start_norm':R.enc(q),'collision_parameters':[{'wall':i,'t':R.enc(t)}for i,t in collisions]}
            if first is None:kind='open'
            else:
                endq=R.add(q,R.mul(R.mul(first[1],first[1]),dualq));s=R.sign(endq)
                if s>0:kind='open'
                else:
                    endpoint=R.vadd(v,R.vm(dual,R.scale(first[1],-1)));ep=[A.pair(endpoint,Q,w)for w in roots]
                    assert all(R.sign(z)<=0 for z in ep)and R.sign(A.pair(endpoint,Q,h))<0
                    active=[i for i,z in enumerate(ep)if not z]
                    kind='finite'if s<0 else'ideal'
                    row.update({'endpoint':encvec(endpoint),'endpoint_norm':R.enc(endq),'endpoint_active':active})
                    if kind=='finite':row['endpoint_vertex']=byactive[tuple(active)]
            row['kind']=kind
            # A concrete spacelike continuation within the whole 12-wall cone
            # is supplied for each open direction by a rational bracket on the
            # null parameter; no new radical extension is needed.
            if kind=='open':
                for bits in range(8,161,8):
                    D=2**bits;lo=0;hi=1
                    while R.sign(R.add(q,R.scale(dualq,F(hi*hi,D*D))))<=0:hi*=2
                    while hi-lo>1:
                        mid=(hi+lo)//2
                        if R.sign(R.add(q,R.scale(dualq,F(mid*mid,D*D))))>0:hi=mid
                        else:lo=mid
                    t=F(hi,D);epoint=R.vadd(v,R.vs(dual,-t));ep=[A.pair(epoint,Q,w)for w in roots]
                    if all(R.sign(ep[i])<0 for i in range(12)if i not in support):break
                else:raise AssertionError('no certified open continuation')
                eq=A.pair(epoint,Q,epoint);assert R.sign(eq)>0 and R.sign(A.pair(epoint,Q,h))<0
                assert [i for i,z in enumerate(ep)if not z]==support
                row['spacelike_continuation']={'rational_parameter':str(t),'point':encvec(epoint),'norm':R.enc(eq),'all_wall_pairings':encvec(ep)}
            if released in range(7,11):
                assert kind=='open';key=tuple(support)
                M=[[gram[i][j]for j in support]for i in support];p,t=spherical_profiles(M,support)
                profile_dict[str(len(edges))]={'finite_profiles':p,'dihedral_tail_profiles':t}
                targetstats['finite_profiles']+=len(p);targetstats['dihedral_tails']+=len(t)
            edges.append(row)
    assert len(vertices)==10 and not any(11 in v['active']for v in vertices)
    assert stat['supports']==220 and len(edges)==90
    kinds=Counter(r['kind']for r in edges)
    closed={tuple(sorted((r['vertex'],r['endpoint_vertex'])))for r in edges if r['kind']=='finite'}
    assert kinds['finite']==2*len(closed)
    # A future timelike vector in the interior of the positive normal cone
    # makes its negative-one affine slice bounded. Thus all projective extreme
    # rays are exhausted by the same 220 rank-nine supports, even when their
    # incident Gram is indefinite and the ray lies outside hyperbolic space.
    weights=[A.k(z)for z in d['positive_dual_center_coefficients']]+[R.scale(h[i],-1)for i in range(6,10)]
    assert len(weights)==11 and all(R.sign(z)>0 for z in weights)
    for power in range(1,31):
        allweights=weights+[R.C(F(1,2**power))];dualcenter=combine(allweights,roots);dn=A.pair(dualcenter,Q,dualcenter)
        if R.sign(dn)<0 and R.sign(A.pair(dualcenter,Q,h))<0:break
    else:raise AssertionError('no bounded future affine slice')
    allcov=[A.mv(Q,r)for r in roots];rays={};raystats=Counter()
    for ss in combinations(range(12),9):
        raystats['supports']+=1;v=null_ray([allcov[i]for i in ss])
        if v is None:raystats['rank_deficient']+=1;continue
        ps=[A.pair(v,Q,r)for r in roots]
        if all(R.sign(z)>=0 for z in ps):v=R.vs(v,-1);ps=[R.scale(z,-1)for z in ps]
        if any(R.sign(z)>0 for z in ps):raystats['infeasible']+=1;continue
        scale=R.scale(A.pair(v,Q,dualcenter),-1);assert R.sign(scale)>0
        v=R.vm(v,inv(scale));assert A.pair(v,Q,dualcenter)==R.C(-1)
        active=tuple(i for i,z in enumerate(ps)if not z);qn=A.pair(v,Q,v);s=R.sign(qn);kind={-1:'finite',0:'ideal',1:'spacelike'}[s]
        row={'active':list(active),'point':encvec(v),'norm':R.enc(qn),'kind':kind,'all_wall_pairings':encvec([A.pair(v,Q,r)for r in roots])}
        if active in rays:assert rays[active]==row
        else:rays[active]=row
    raycounts=Counter(r['kind']for r in rays.values());assert raycounts['finite']==10
    assert {tuple(r['active'])for r in rays.values()if r['kind']=='finite'}==set(byactive)
    projective={'support_census':dict(raystats),'counts':dict(raycounts),'normal_cone_positive_weights':encvec(allweights),'timelike_slice_dual':encvec(dualcenter),'slice_dual_norm':R.enc(dn),'rays':[rays[k]for k in sorted(rays)]}
    out={'status':'PASS_EXACT_COMPLETE_FINITE_VERTEX_STAR_AND_SPHERICAL_PROFILE_QUEUE','sources':{p:R.sha(p)for p in [A.SOURCE,A.SCRIPT,'scripts/'+Path(A.__file__).name,'scripts/fresh_h22_multiquadratic_20261013_review.py','scripts/'+Path(__file__).name]},'finite_vertex_support_census':dict(stat),'finite_vertices':vertices,'directed_edge_counts':dict(kinds),'finite_edges':len(closed),'outgoing_edges':edges,'H4_released_edge_profiles':profile_dict,'H4_profile_counts':dict(targetstats),'complete_projective_ray_census':projective,'profile_completeness':'Every finite rank-nine Coxeter stabilizer is spherical. Its diagram is a forest, has maximum vertex degree three, and a degree-three vertex has only label-three edges. Degree-two attachments use labels at most five. Unbounded labels occur only when the new vertex joins an isolated old A1, giving a separate I2(m) tail. These cases are exhaustively tested by the exact Schur complement on each positive-definite support. Profiles are necessary local possibilities, not globally compatible caps.','scope':'Complete projective extreme-ray census for these twelve actual walls, together with all outgoing directions at finite vertices and exact feasible spacelike continuations for open directions. These data describe the present chamber, not an all-parameter completion obstruction or a finite-volume example.'}
    if verify:assert out==R.read(OUT)
    else:(R.P/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(vertices),dict(kinds),dict(targetstats),dict(raycounts))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
