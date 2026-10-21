#!/usr/bin/env python3
"""Compare both endpoint caps for all five native wall66 alternatives.

Complete finite-profile shapes, including explicitly isolated label tails;
bounded numerical root/pair discovery, not an exact geometry certificate.
"""
import argparse,json
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_nonarithmetic_release30_66_hinge_20261017 as H
import fresh_nonarithmetic_release49_endpoint_gate_20261017 as G

ROOT,B,M=H.ROOT,H.B,H.M
TAG='fresh_nonarithmetic_alternative66_endcaps_20261017'
CHOICES=[[2,2,4],[3,2,3],[3,2,4],[4,2,3],[4,2,4]]
FRAMES=[[33,30,66],[58,60,66]]

@lru_cache(None)
def components(key):
    n=int(len(key)**.5);g=[key[i*n:(i+1)*n]for i in range(n)];left=set(range(n));out=[]
    while left:
        cc=set();todo=[min(left)]
        while todo:
            a=todo.pop()
            if a in cc:continue
            cc.add(a);todo.extend(b for b in left if g[a][b]>2 and b not in cc)
        left-=cc;ids=sorted(cc)
        if len(ids)==2:out.append('I2('+str(g[ids[0]][ids[1]])+')');continue
        code=G.canonical([[g[i][j]for j in ids]for i in ids]);name=G.CAT.get((len(ids),code))
        if name is None:return None
        out.append(name)
    return sorted(out)

def corner(frame,choice,pi,ms,pair):
    a,b=pair;anchors=[frame[a],frame[b]];g=[[2]*5 for _ in range(5)]
    for i in range(5):g[i][i]=1
    g[0][1]=g[1][0]=3
    for j,i in enumerate(anchors):
        for k,m in enumerate(G.profile_labels(6 if i==60 else 0)):g[k][j+2]=g[j+2][k]=m
    for k,m in enumerate(G.profile_labels(pi)):g[k][4]=g[4][k]=m
    old=choice[2]if set(anchors)=={33,66}else choice[0]if set(anchors)=={58,66}else 2
    g[2][3]=g[3][2]=old
    g[2][4]=g[4][2]=ms[a];g[3][4]=g[4][3]=ms[b]
    return components(tuple(v for row in g for v in row))

def domain(choice,frame):
    low=[];tails=[];counts=Counter()
    for pi in range(7):
        for ms in product(range(2,7),repeat=3):
            ty=[corner(frame,choice,pi,ms,p)for p in[(0,1),(0,2)]]
            if any(v is None for v in ty):counts['nonfinite']+=1;continue
            high=[j for j,m in enumerate(ms)if m==6]
            assert len(high)<=1
            row=dict(profile=pi,labels=list(ms),corner_types=ty)
            if high:row['unbounded_label_position']=high[0];tails.append(row)
            else:low.append(row)
    return dict(low=low,tails=tails,counts=dict(counts))

def frame_solver(frame,labels,ns,h,profiles,recipe,tail_c=None):
    f=np.array([ns[i]for i in frame]);gram=f@M@f.T
    if abs(np.linalg.det(gram))<1e-10:return []
    inv=np.linalg.inv(gram);k=np.linalg.svd(f@M)[2][-1];kq=k@M@k
    if abs(kq)<1e-10:return []
    p=profiles[recipe['profile']];off=np.array(p['offsets'])
    native=-np.array([B.theta(m)for m in recipe['labels']])
    if tail_c is not None:native[recipe['unbounded_label_position']]=-tail_c
    rhs=native-off[[labels.index(i)for i in frame]];part=rhs@inv@f
    rem=(p['square']-part@M@part)/kq
    if rem<=1e-10:return []
    out=[]
    for sign in[-1,1]:
        y=part+sign*np.sqrt(rem)*k;vv=np.array([y@M@ns[i]+off[labels.index(i)]for i in labels])
        for j,i in enumerate(frame):vv[labels.index(i)]=native[j]
        rel={i:B.relation(x)for i,x in zip(labels,vv)}
        for j,i in enumerate(frame):rel[i]=recipe['labels'][j]if tail_c is None or j!=recipe.get('unbounded_label_position')else 'assigned_continuous_tail'
        bad=[i for i in labels if rel[i]is None or rel[i]=='parallel_or_unresolved']
        out.append(dict(branch=sign,projection=y.tolist(),center_product=float(y@M@h+p['center_offset']),
            blockers=[dict(wall=i,product=float(vv[labels.index(i)]))for i in bad],
            finite_contacts=[[i,r]for i,r in rel.items()if isinstance(r,int)],
            nonframe_dotted=all(rel[i]=='dotted'for i in labels if i not in frame),
            max_nonframe_native=float(max(vv[j]for j,i in enumerate(labels)if i not in frame))))
    return out

def shared_corner(a,b,profiles,mutual):
    g=[[2]*5 for _ in range(5)]
    for i in range(5):g[i][i]=1
    g[0][1]=g[1][0]=3
    for at,row in[(3,a),(4,b)]:
        for j,m in enumerate(G.profile_labels(row['profile'])):g[j][at]=g[at][j]=m
        g[2][at]=g[at][2]=dict(row['finite_contacts'])[66]
    g[3][4]=g[4][3]=mutual
    return components(tuple(v for row in g for v in row))

def run():
    z,labels,ns,h,profiles,primary=H.setup();old=H.A.setup()[0]
    gram=np.array([[B.mid(x)for x in row]for row in old['residual_gram']]);vals,vecs=np.linalg.eigh((gram+gram.T)/2);tr=vecs@np.diag(1/np.sqrt(abs(vals)))
    saved=json.loads((ROOT/'results'/f'{G.TAG}.json').read_text());proofs=saved['records'][1]['exact_recipes']
    squares={r['wall']:B.mid(r['square'])for r in z['projections']};stages=[];all_improvements=[]
    for choice in CHOICES:
        pr=next(r for r in proofs if r['profile']==0 and r['labels']==choice)
        branch=next(r for r in pr['branches']if r['sign']==-1)
        ns[66]=np.linalg.solve(tr,np.array([B.mid(v)for v in branch['coordinates']]))
        actual=z['actual_sides'];control=B.ledger(np.array([ns[i]for i in actual]),actual)
        assert(control['finite'],control['exterior'],control['open_ends'])==(26,6,10)
        groups=[];domains=[];tail_scans=[];single=[];rejected=[];counts=Counter()
        for group,frame in enumerate(FRAMES):
            d=domain(choice,frame);domains.append(d);pool=[]
            for recipe in d['low']:
                roots=frame_solver(frame,labels,ns,h,profiles,recipe);counts['real_low_branches']+=len(roots)
                for r in roots:
                    rec=dict(**r,profile=recipe['profile'],assigned_labels=recipe['labels'],frame=frame,group=group,corner_types=recipe['corner_types'])
                    if r['blockers']or r['center_product']>=-1e-8:rejected.append(rec);continue
                    rec['index']=len(single);y=np.array(r['projection']);p=profiles[recipe['profile']]
                    rec['acute_incidence']=H.incidence_score([ns[i]for i in actual]+[y],actual+[67],[squares[i]for i in actual]+[p['square']])
                    rec['projective_boundary']=B.ledger(np.array([ns[i]for i in actual]+[y]),actual+[67])
                    single.append(rec);pool.append(rec['index'])
            groups.append(pool)
            # A fixed compact parameter scout, not a label cutoff or proof.
            # Endpoints c=sqrt3 and c=2 include the limiting parallel value.
            for recipe in d['tails']:
                samples=[]
                for c in np.linspace(np.sqrt(3.),2.,33):
                    for r in frame_solver(frame,labels,ns,h,profiles,recipe,float(c)):
                        samples.append(dict(c=float(c),branch=r['branch'],center_product=r['center_product'],
                            blockers=r['blockers'],all_nonframe_dotted=r['nonframe_dotted']))
                tail_scans.append(dict(group=group,frame=frame,recipe=recipe,samples=samples,
                    scope='33 evenly spaced continuous cosine parameters, including limiting c=2; not discrete Coxeter normals except when separately exactified, not an exhaustive interval gate.'))
        pairs=[];relation_counts=Counter()
        for ia,ib in product(*groups):
            a,b=single[ia],single[ib];ya,yb=np.array(a['projection']),np.array(b['projection'])
            native=H.A.native(ya,a['profile'],yb,b['profile'],profiles);rel=B.relation(native);relation_counts[str(rel)]+=1
            if not isinstance(rel,int):continue
            ty=shared_corner(a,b,profiles,rel)
            if ty is None:continue
            rr=[ns[i]for i in actual]+[ya,yb];labs=actual+[67,68];ss=[squares[i]for i in actual]+[profiles[a['profile']]['square'],profiles[b['profile']]['square']]
            score=H.incidence_score(rr,labs,ss);boundary=B.ledger(np.array(rr),labs)
            rec=dict(indices=[ia,ib],mutual_label=rel,native_product=native,shared66_corner_types=ty,acute_incidence=score,projective_boundary=boundary)
            pairs.append(rec)
            if score.get('open_ends',99)<10:all_improvements.append(dict(choice=choice,**rec,roots=[a,b]));print('IMPROVEMENT',choice,[ia,ib],score['open_ends'],flush=True)
        stages.append(dict(choice=choice,control=control,finite_domains=domains,groups=groups,counts=dict(counts),
            accepted_low_caps=single,rejected_low_caps=rejected,tail_continuous_scans=tail_scans,pair_relation_counts=dict(relation_counts),finite_shared_corner_pairs=pairs))
        print('ALTERNATIVE66',choice,'domains',[(len(d['low']),len(d['tails']))for d in domains],'caps',list(map(len,groups)),'pairs',dict(relation_counts),'valid',len(pairs),flush=True)
    paths=[Path(__file__),Path(H.__file__),Path(G.__file__),H.PARENT,H.BOUNDARY,H.A.INPUT,H.A.C.PARENT,ROOT/'results'/f'{G.TAG}.json']
    return dict(status='BOUNDED_NUMERICAL_FIVE_ALTERNATIVE66_ENDPOINT_COMPARISON',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        choices=CHOICES,frames=FRAMES,white_profiles=profiles,stages=stages,improvements=all_improvements,
        finite_profile_scope='Complete rank<=5 endpoint shapes: enumerate labels2..5 plus symbolic6 representing every isolated I2(m),m>=6. A high label belongs only to a two-node connected component; the common cap cannot have two high contacts in one finite corner, and the fixed33/66 or60/A2 bonds rule out the remaining possibilities. Every high family remains explicit.',
        scope='Native midpoint root/gate tests and numerical complete-boundary diagnostics only. No exact new geometry or all-parameter exclusion. Continuous high-label samples are not finite Coxeter root certificates. A shared66 finite corner is required before a low cap pair is nominated; exact affine screening and all physical gates remain prerequisites to promotion.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/(TAG+'.json')
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'improvements',len(out['improvements']),flush=True)
