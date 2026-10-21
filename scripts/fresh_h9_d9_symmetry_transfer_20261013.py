#!/usr/bin/env python3
"""Exact D9 chamber transport; depends explicitly on complete first-edge audits.

No producer field arithmetic, numerical search, or PARI invocation is imported.
The existing independent Fraction field engine supplies only exact operations.
"""
import argparse
import hashlib
import json
from pathlib import Path
import fresh_h9_glued_host_20261013_independent as a

ROOT=Path(__file__).resolve().parents[1]
PREFIX='fresh_h9_d9_symmetry_transfer_20261013'
OUT=ROOT/'results'/f'{PREFIX}.json'
OLD='results/fresh_h9_changed_spatial_20261013'
NEW='results/fresh_h9_changed_spatial_next_20261013'
G=tuple(tuple(a.neg(a.ALPHA)if i==j==0 else a.O if i==j else a.Z for j in range(10))for i in range(10))
IDENT=tuple(tuple(a.O if i==j else a.Z for j in range(10))for i in range(10))

def total(xs):
    z=a.Z
    for x in xs:z=a.add(z,x)
    return z

def pair(x,y):return a.dot(x,a.mv(G,y))
def matrix(x):return tuple(a.vec(r)for r in x)
def encvec(x):return list(map(a.enc,x))
def encmat(x):return list(map(encvec,x))
def tr(A):return total(A[i][i]for i in range(len(A)))
def transpose(A):return tuple(zip(*A))
def mm(A,B):return tuple(tuple(a.dot(r,c)for c in zip(*B))for r in A)
def vmul(s,v):return tuple(a.mul(s,x)for x in v)
def lattice(v):return all(map(a.integral,v))and a.integral(a.scale(total(v[1:]),a.Q(1,2)))
def basis():
    out=[(a.O,*(a.Z for _ in range(9)))]
    for i in range(1,9):
        r=[a.Z]*10;r[i]=a.O;r[i+1]=a.neg(a.O);out.append(tuple(r))
    r=[a.Z]*10;r[8]=r[9]=a.O;out.append(tuple(r))
    return tuple(out)
BASIS=basis();B=transpose(BASIS);BI=a.minverse(B)
INITIAL=tuple(tuple(a.neg(x)for x in r)for r in BASIS[1:9])+((a.Z,*(a.Z for _ in range(8)),a.scale(a.O,-2)),)


def primitive_reflective(v):
    assert lattice(v) and not lattice(tuple(a.divide(x,a.T)for x in v))
    q=pair(v,v);assert all(a.sign(q,e)>0 for e in range(4))
    assert all(a.integral(a.divide(a.scale(pair(v,b),2),q))for b in BASIS)


def positive_definite(C):
    M=[list(r)for r in C];pivots=[]
    for k in range(len(M)):
        p=M[k][k];assert a.sign(p)>0;pivots.append(p)
        for i in range(k+1,len(M)):
            for j in range(k+1,len(M)):
                M[i][j]=a.sub(M[i][j],a.divide(a.mul(M[i][k],M[k][j]),p))
    return pivots


def simple_system(roots,point):
    assert len(roots)==9 and len(set(roots))==9
    assert a.sign(pair(point,point))<0 and a.sign(point[0])>0
    C=tuple(tuple(pair(x,y)for y in roots)for x in roots)
    for i,r in enumerate(roots):
        primitive_reflective(r);assert pair(point,r)==a.Z
        for j in range(i):
            assert a.sign(C[i][j])<=0
            ratio=a.divide(a.scale(a.sq(C[i][j]),4),a.mul(C[i][i],C[j][j]))
            assert ratio in a.ALLOWED
    return C,positive_definite(C)


def mirrors(roots):
    # Literal roots are stable under the finite reflection group; opposite
    # vectors represent the same mirror. Norms are exactly preserved.
    found=set(roots);queue=list(roots)
    for r in queue:
        for s in roots:
            c=a.divide(a.scale(pair(r,s),2),pair(s,s))
            w=tuple(a.sub(x,a.mul(c,y))for x,y in zip(r,s,strict=True))
            if w not in found:found.add(w);queue.append(w)
        assert len(found)<=200
    assert all(tuple(map(a.neg,r))in found for r in found)
    return len(found)//2


def checked_sources(paths):
    sources={};seen=set()
    def visit(rel):
        if rel in seen:return
        seen.add(rel);p=ROOT/rel
        sources[rel]=hashlib.sha256(p.read_bytes()).hexdigest()
        if p.suffix!='.json':return
        d=json.loads(p.read_text())
        for key in ('sources','sources_sha256'):
            manifest=d.get(key)
            if isinstance(manifest,dict):
                for q,h in manifest.items():
                    assert isinstance(h,str)and hashlib.sha256((ROOT/q).read_bytes()).hexdigest()==h,(rel,q)
                    visit(q)
    for p in paths:visit(p)
    return dict(sorted(sources.items()))


def path_step(start_roots,start,end_roots,end,newroot,label):
    edge=tuple(r for r in start_roots if r in end_roots)
    released=tuple(r for r in start_roots if r not in end_roots)
    arriving=tuple(r for r in end_roots if r not in start_roots)
    assert len(edge)==8 and len(released)==len(arriving)==1 and arriving[0]==newroot
    assert all(pair(start,r)==pair(end,r)==a.Z for r in edge)
    assert a.sign(pair(start,newroot))<0 and a.sign(pair(end,released[0]))<0
    assert a.sign(pair(start,end))<0
    simple_system(start_roots,start);simple_system(end_roots,end)
    return {'label':label,'edge_normals':encmat(edge),'released_normal':encvec(released[0]),
            'arriving_normal':encvec(newroot),'start':encvec(start),'end':encvec(end),
            'start_on_arriving_wall':a.enc(pair(start,newroot)),
            'end_on_released_wall':a.enc(pair(end,released[0])),
            'same_oriented_edge_roots':True,
            'open_segment_mirror_free_dependency':'The hash-bound complete independent first-edge census, including simultaneous-root cone choices.'}


def supplementary(ck):
    H=matrix(ck['opposite_B9_transport']['matrix'])
    J,K=[matrix(row['matrix'])for row in ck['B2_diagram_swap_tests']]
    assert mm(H,H)==mm(J,J)==mm(K,K)==IDENT
    assert mm(J,K)==mm(K,J) and mm(mm(H,J),H)==K
    group=[IDENT,H,J,K,mm(J,K),mm(H,J),mm(H,K),mm(mm(H,J),K)]
    assert len(set(group))==8
    assert all(mm(x,y)in group for x in group for y in (H,J,K))
    preserving=[i for i,x in enumerate(group)if all(lattice(a.mv(x,b))for b in BASIS)]
    assert preserving==[0]
    # Canonical odd unimodular overclosure, not needed for the global proof.
    half=a.scale(a.O,a.Q(1,2));s=(a.Z,*(half for _ in range(9)))
    assert all(a.integral(pair(s,b))for b in BASIS)
    assert not lattice(s)and not lattice(vmul(a.scale(a.O,2),s))and lattice(vmul(a.scale(a.O,4),s))
    D=tuple(tuple(pair(x,y)for y in BASIS)for x in BASIS)
    di=a.minverse(D)
    # Each dual basis vector differs from a multiple of s by an L-vector.
    for col in transpose(di):
        v=a.mv(B,col);coef=a.scale(v[1],2)
        assert a.integral(coef)
        assert any(lattice(tuple(a.sub(x,a.mul(a.add(coef,a.scale(a.O,j)),y))for x,y in zip(v,s)))for j in (0,2))
    return {'previous_three_involution_group_order':8,'only_D9_preserving_element_is_identity':True,
            'discriminant_generator':encvec(s),'discriminant_module':'O/(4)',
            'canonical_overclosure':'<-alpha> + O^9',
            'overclosure_argument':'The discriminant module is cyclic O/(4), generated by s=(0,1/2,...,1/2). The unique intermediate submodule of order Norm(2) is generated by 2s because 2 is totally ramified. Its inverse image is exactly the ambient integral lattice. Hence every D9 isometry preserves this overclosure.',
            'scope':'These supplementary calculations do not identify the odd overclosure with an earlier F4 lattice. They are not needed for the direct D9 chamber transport.'}


def run(verify=False):
    paths=[OLD+'.json',OLD+'_step2.json',OLD+'_step3.json',OLD+'_independent.json',
           OLD+'_checkpoint.json',NEW+'.json',NEW+'_independent.json']
    sources=checked_sources(paths+['scripts/fresh_h9_glued_host_20261013_independent.py'])
    source=[json.loads((ROOT/p).read_text())for p in paths]
    first,second,third,chain_audit,ck,new,new_audit=source
    assert chain_audit['status']=='PASS_INDEPENDENT_D9_COMPLETE_EDGE_CHAIN'
    for i,d in enumerate((first,second,third),1):
        audit=next(r for r in chain_audit['stages']if r['stage']==i)
        assert audit['profiles']==(len(d['profiles'])if i==1 else d['profiles'])
        assert audit['interval_points']==len(d['visited_points'])
        assert chain_audit['sources'][paths[i-1]]==sources[paths[i-1]]
    assert new_audit['status']=='PASS_INDEPENDENT_COMPLETE_D9_A2_B6_FIRST_EDGE'
    assert new_audit['source_sha256']==sources[NEW+'.json']
    assert (new_audit['profiles'],new_audit['interval_points'],new_audit['first_mirrors'])==(28,20,2)
    assert matrix(first['initial_walls'])==INITIAL
    origin=(a.O,*(a.Z for _ in range(9)))
    oldroots=INITIAL[1:]+(a.vec(first['first_mirror']),)
    oldpoint=a.vec(first['first_endpoint'])
    path=[path_step(INITIAL,origin,oldroots,oldpoint,oldroots[-1],'initial B9 to chain corner 1')]
    for stage,d in ((2,second),(3,third)):
        roots=matrix(d['endpoint_incident']);point=a.vec(d['endpoint'])
        path.append(path_step(oldroots,oldpoint,roots,point,a.vec(d['chosen_root']),f'chain corner {stage-1} to {stage}'))
        oldroots,oldpoint=roots,point
    newroots=matrix(new['endpoint_incident']);newpoint=a.vec(new['endpoint'])
    path.append(path_step(INITIAL,origin,newroots,newpoint,a.vec(new['chosen_root']),'initial B9 to new B3+B6 corner'))
    # The new prefix contains both endpoints and all four entire edge segments:
    # homogeneous halfspaces are convex along the normalized timelike segment.
    walls=matrix(new['all_walls'])
    for row in path:
        for key in('start','end'):
            p=a.vec(row[key]);assert all(a.sign(pair(p,r))<=0 for r in walls)
    perm=(8,1,2,3,4,5,6,7,0)
    pointscale=(a.Q(-23),a.Q(-4),a.Q(30),a.Q(14))
    assert a.sign(pointscale)>0
    assert pair(oldpoint,oldpoint)==a.mul(a.sq(pointscale),pair(newpoint,newpoint))
    C,pivots=simple_system(oldroots,oldpoint);CN,pivotsn=simple_system(newroots,newpoint)
    assert all(C[i][j]==CN[perm[i]][perm[j]]for i in range(9)for j in range(9))
    oldframe=transpose((*oldroots,oldpoint))
    newframe=transpose((*(newroots[i]for i in perm),vmul(pointscale,newpoint)))
    F=mm(newframe,a.minverse(oldframe));FI=a.minverse(F)
    assert mm(F,FI)==mm(FI,F)==IDENT and mm(mm(transpose(F),G),F)==G
    assert a.mv(F,oldpoint)==vmul(pointscale,newpoint)
    assert all(a.mv(F,r)==newroots[perm[i]]for i,r in enumerate(oldroots))
    module=mm(mm(BI,F),B);module_inverse=mm(mm(BI,FI),B)
    assert all(a.integral(x)for A in (module,module_inverse)for row in A for x in row)
    assert all(lattice(a.mv(A,b))for A in(F,FI)for b in BASIS)
    trace=tr(F);expected=tuple(map(a.Q,(-7378736,-3993344,12596300,6817068)))
    assert trace==expected and a.sign(a.sub(trace,a.scale(a.O,10)))>0
    old_count=mirrors(oldroots);new_count=mirrors(newroots)
    assert old_count==new_count==45
    forged=[list(r)for r in F];forged[0][0]=a.add(forged[0][0],a.O)
    assert mm(mm(transpose(forged),G),forged)!=G
    result={'status':'PASS_EXACT_INFINITE_FULL_D9_CHAMBER_STABILIZER','sources':sources,
            'arithmetic':'Exact Fraction arithmetic from the pre-existing independent quartic engine; no producer imports. Complete edge censuses and their certified PARI norm/unit interface are explicit hash-bound dependencies.',
            'host':'<-(T-1)> + D9 over O=Z[T], T^4-4T^2+2=0, T=2cos(pi/8)',
            'gram':encmat(G),'D9_basis':encmat(BASIS),'initial_zero_chamber':'B9, 81 mirrors',
            'oriented_path':path,'required_complete_edges':4,
            'old_complete_corner_type':'B3+B6','new_complete_corner_type':'B3+B6',
            'old_finite_mirrors':old_count,'new_finite_mirrors':new_count,
            'old_corner':encvec(oldpoint),'new_corner':encvec(newpoint),
            'old_incident':encmat(oldroots),'new_incident':encmat(newroots),
            'old_Gram':encmat(C),'new_Gram':encmat(CN),'old_LDL':encvec(pivots),'new_LDL':encvec(pivotsn),
            'permutation':list(perm),'positive_point_scale':a.enc(pointscale),
            'transport_matrix':encmat(F),'inverse_matrix':encmat(FI),
            'module_matrix':encmat(module),'module_inverse':encmat(module_inverse),
            'trace':a.enc(trace),'trace_above_dimension':True,
            'exact_oriented_primitive_root_mapping':True,'future_preserving':True,
            'preserves_D9_and_inverse':True,'forged_matrix_rejected':True,
            'supplementary':supplementary(ck),
            'conclusion':'The four complete nearest-edge steps connect two full finite corners of one D9 reflection chamber. The D9 lattice isometry maps their oriented complete simple systems, so preserves that full chamber. Its trace exceeds10, hence it has infinite order. Discreteness then forces infinite chamber volume. No finite-volume reflection subgroup using D9-lattice reflections exists.',
            'scope':'This exact D9 lattice is nonreflective. No conclusion about all lattices in its rational quadratic space, other arithmetic hosts, nonarithmetic groups, or compact H9 in general.'}
    if verify:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'trace',result['trace'],'required edges',len(path),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
