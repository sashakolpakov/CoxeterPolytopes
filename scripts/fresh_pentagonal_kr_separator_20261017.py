#!/usr/bin/env python3
"""Numerical S-cap search replacing the dotted K--R gap by a quadrilateral.

All saved dotted-KR pentagonal prefixes are tested. Four prescribed contacts
determine S linearly; its required norm is a separate gate. This is bounded
discovery, not an exact or all-label nonexistence result.
"""
import argparse
import json
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import numpy as np
import fresh_cu_release47_pentagonal_cap_20261017 as N

ROOT=N.ROOT;TAG='fresh_pentagonal_kr_separator_20261017'
ALPHABET=(2,3,4,5,6,8,12,24)
FRAME=(68,72,70,71)  # A,K,P,R in cyclic order around S.


def relation(x):
    old=N.B.relation(float(x))
    if not isinstance(old,int):
        return old
    if old==2:
        return old
    gap=min(abs(N.B.theta(old)-N.B.theta(old-1)),
            abs(N.B.theta(old+1)-N.B.theta(old)))
    # Unprescribed finite contacts must resolve the local label spacing.
    if abs(x+N.B.theta(old)) > min(2e-7,gap/32):
        return None
    return old


def disk():
    prior=N.disk();removed=({70,71,72},{71,68,72})
    fs=[t for t in prior['triangles'] if set(t) not in removed]
    fs += [[70,71,73],[70,73,72],[68,71,73],[68,73,72]]
    es=Counter(tuple(sorted(e)) for t in fs for e in combinations(t,2))
    vs=set(sum(fs,[]));bc=prior['boundary']
    assert (len(vs),len(es),len(fs))==(10,22,13)
    assert len(vs)-len(es)+len(fs)==1 and (71,72) not in es
    assert {e for e,c in es.items() if c==1}=={tuple(sorted(e)) for e in zip(bc,bc[1:]+bc[:1])}
    assert set(es.values())=={1,2}
    links=[]
    for v in sorted(vs):
        ee=[tuple(i for i in t if i!=v) for t in fs if v in t]
        ds=Counter(i for e in ee for i in e);seen={next(iter(ds))}
        while True:
            new=seen|{b for a,b in ee if a in seen}|{a for a,b in ee if b in seen}
            if new==seen:break
            seen=new
        assert seen==set(ds)
        assert sorted(ds.values())==([1,1]+[2]*(len(ds)-2) if v in bc else [2]*len(ds))
        links.append(dict(vertex=v,edges=ee))
    return dict(boundary=bc,triangles=fs,vertices=10,edges=22,faces=13,
                removed_edge=[71,72],S_link=[70,71,68,72,70],links=links)


@lru_cache(None)
def recipes(types,bonds):
    """All finite-corner label tuples in the explicitly bounded alphabet."""
    out=[]
    for sp in range(7):
        pairs=[]
        for i in range(4):
            j=(i+1)%4
            pairs.append({(m,n) for m in ALPHABET for n in ALPHABET
                          if N.P.corner3(types[i],types[j],sp,bonds[i],m,n) is not None})
        for a,k in sorted(pairs[0]):
            for kk,p in sorted(pairs[1]):
                if kk!=k:continue
                for pp,r in sorted(pairs[2]):
                    if pp==p and (r,a) in pairs[3]:
                        out.append((sp,(a,k,p,r)))
    return tuple(out)


def run():
    path=ROOT/f'results/{N.TAG}.json';old=json.loads(path.read_text())
    topology=disk()
    cur,oldlabels,base,h,profiles,*_=N.U.setup()
    base.pop(47);base.pop(67)
    oi={i:j for j,i in enumerate(oldlabels)}
    actual=[i for i in cur['actual_sides'] if i!=47]
    off=lambda a,b:float(N.H.A.native(np.zeros(4),a,np.zeros(4),b,profiles))
    counts=Counter();fail=Counter();unresolved=[];candidates=[];near=[];domains={}
    domain_arrays={};minres=float('inf');worst_cond=0.;tested=[]
    def native(v,sp,nd,ty):
        vals={i:float(v@N.M@w)+(profiles[sp]['offsets'][oi[i]] if i in oi else off(sp,ty.get(i,0))) for i,w in nd.items()}
        bad=[i for i,x in vals.items() if relation(x) in (None,'parallel_or_unresolved')]
        return vals,bad,float(v@N.M@h)+profiles[sp]['center_offset']
    for r in old['R_caps']:
        if r['KR_relation']!='dotted':continue
        rid=r['index'];tested.append(rid);e=old['events'][r['event']]
        nd=dict(base);nd.update({int(i):np.array(v) for i,v in e['normals'].items()});nd[71]=np.array(r['normal'])
        ty={int(i):p for i,p in e['profiles'].items()};ty[71]=r['profile']
        types=tuple(ty[i] for i in FRAME)
        k=old['K_seeds'][e['K']]
        bonds=(e['AK'],k['labels'][2],r['labels'][2],r['labels'][1])
        key=(types,bonds)
        if key not in domain_arrays:
            dom=recipes(types,bonds)
            rhs=np.array([[-N.B.theta(m)-off(sp,tp) for m,tp in zip(ms,types)] for sp,ms in dom])
            norms=np.array([profiles[sp]['square'] for sp,ms in dom])
            domain_arrays[key]=(dom,rhs,norms)
            domains[key]=dict(profiles=types,cycle_bonds=bonds,recipes=[dict(profile=sp,labels=ms) for sp,ms in dom])
        dom,rhs,norms=domain_arrays[key]
        counts['prefixes']+=1;counts['four_contact_recipes']+=len(dom)
        if not dom:continue
        frame=np.array([nd[i] for i in FRAME]);mat=frame@N.M
        cond=float(np.linalg.cond(mat));worst_cond=max(worst_cond,cond)
        if not np.isfinite(cond) or cond>1e10:
            unresolved.append(dict(R=rid,reason='ill_conditioned_four_frame',condition=cond));continue
        inv=np.linalg.inv(mat)
        normals=rhs@inv.T
        squares=np.einsum('ij,jk,ik->i',normals,N.M,normals)
        residuals=np.abs(squares-norms)
        best=int(np.argmin(residuals));res=float(residuals[best]);minres=min(minres,res)
        if len(near)<24 or res<near[-1]['absolute_norm_residual']:
            sp,ms=dom[best]
            near.append(dict(R=rid,event=e['index'],profile=sp,labels=ms,
                absolute_norm_residual=res,computed_square=float(squares[best]),required_square=float(norms[best]),condition=cond))
            near=sorted(near,key=lambda z:z['absolute_norm_residual'])[:24]
        for j in np.flatnonzero(residuals<=2e-7):
            counts['near_norm_candidates']+=1;sp,ms=dom[int(j)];v=normals[j]
            vals,bad,cp=native(v,sp,nd,ty)
            if bad or cp>=-1e-8:
                fail[('S_native',tuple(bad),'center' if cp>=-1e-8 else '')]+=1;continue
            # Recheck every new prefix wall with the spacing-aware native gate.
            gatebad=[]
            for i in [67,68,69,70,71,72]:
                _,bb,cc=native(nd[i],ty.get(i,0),{q:w for q,w in nd.items() if q!=i},ty)
                if bb or cc>=-1e-8:gatebad.append((i,bb,cc))
            if gatebad:
                fail[('prefix_native',tuple(i for i,bb,cc in gatebad),'')]+=1;continue
            nn=dict(nd);nn[73]=v;tt=dict(ty);tt[73]=sp
            ids=actual+[i for i in sorted(nn) if i not in oi]
            bd=N.S.C.detailed_cone(nn,ids)
            sq=[profiles[tt.get(i,0)]['square'] if i not in oi else float(nn[i]@N.M@nn[i]) for i in ids]
            ic=N.H.incidence_score([nn[i] for i in ids],ids,sq)
            good=not bd['warnings'] and ic['status']=='numerically_separated_acute_incidence'
            good=good and bd['open_ends']==ic['open_ends'] and bd['counts'].get('finite',0)==ic['finite']
            good=good and {tuple(sorted(z['active'])) for z in bd['rays'] if z['kind']=='finite'}=={tuple(sorted(z)) for z in ic['finite_supports']}
            if not good:
                unresolved.append(dict(R=rid,profile=sp,labels=ms,reason='cone_PD_disagreement',boundary=bd,PD=ic));continue
            sc=N.D.compact_score(dict(boundary=bd,PD_incidence=ic))
            rayset={tuple(sorted(z['active'])) for z in sc['rays'] if z['kind']=='finite'}
            whole={}
            for i in [69,70,71,72,73]:
                expected={tuple(sorted(t)) for t in topology['triangles'] if i in t}
                rays=[z for z in sc['rays'] if i in z['active']]
                whole[i]=all(z['kind']=='finite' for z in rays) and {tuple(sorted(z['active'])) for z in rays}==expected
            complete=all(tuple(sorted(t)) in rayset for t in topology['triangles']) and all(whole.values())
            candidates.append(dict(R=rid,event=e['index'],profile=sp,labels=ms,normal=v.tolist(),
                native_products=sorted(vals.items()),center_product=cp,norm_residual=float(residuals[j]),
                contact_residual=float(np.max(np.abs(mat@v-rhs[j]))),whole_faces=whole,completed_disk=complete,score=sc))
            counts['native_S_candidates']+=1
            if sc['open_ends']<10:counts['numerical_open_end_improvements']+=1
        if counts['prefixes']%2000==0:
            print('KR_SEPARATOR',counts['prefixes'],counts['four_contact_recipes'],counts['near_norm_candidates'],len(candidates),minres,flush=True)
    assert len(tested)==32019
    paths=[Path(__file__),Path(N.__file__),path,Path(N.P.__file__),Path(N.B.__file__),Path(N.U.__file__),Path(N.H.__file__),Path(N.S.__file__)]
    return dict(status='NUMERICAL_PENTAGONAL_KR_SEPARATOR_COMPLETE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        topology=topology,alphabet=ALPHABET,tested_R_records=tested,
        finite_corner_domains=list(domains.values()),counts=dict(counts),
        nearest_norm_rows=near,minimum_absolute_norm_residual=minres,
        maximum_four_frame_condition=worst_cond,candidates=candidates,unresolved=unresolved,
        blockers=[dict(stage=k[0],walls=k[1],other=k[2],count=n) for k,n in sorted(fail.items())],
        scope='All32019 saved dotted-KR prefixes; every seven-profile/four-label tuple in the explicit bounded alphabet satisfying all four finite corner classifications. Four-frame linear solution and numerical norm/native gates; full cone+PD comparison only after all native gates. No exact rejection, all-label exhaustion, ambient finite-volume or new construction dimension claim. Frozen predecessors unchanged.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:
        assert out==json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'unresolved',len(out['unresolved']),flush=True)
