#!/usr/bin/env python3
"""Bounded cap proposals onremaining complete geodesic faces without finite vertices."""
import argparse,gzip,json,math
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import mpmath as mp
import numpy as np
import fresh_h9_polar_endpoint_profiles_20261015 as P
import fresh_h9_quartet_polar_second_20261014 as B
A,D,R,ROOT=B.A,B.D,B.P,B.ROOT
TAG='fresh_h9_after_third_vertex_free_scout_20261015'
FILES=['results/fresh_h9_replacement_escape_extension_20261015.json','results/fresh_h9_after_third_vertex_free_edges_20261015.json']
ANCHORS=[2,3,4,5,6,7,8,9,10,11,12,16,20,24,32,64];TAILS=[6,7,8,10,12,16,24,32,64]
mp.mp.dps=100;S5=mp.sqrt(5);RAD=6+16*S5/5

def scalar(x):return mp.mpf(x[0][0])+mp.mpf(x[0][1])*S5+(mp.mpf(x[1][0])+mp.mpf(x[1][1])*S5)*mp.sqrt(RAD)
def qval(x):return float(mp.mpf(str(x.a))+mp.mpf(str(x.b))*S5)
def positive_root(a,b,c):
    assert a<0 and c>0
    disc=b*b-4*a*c
    # One positive root; use cancellation-safe expression.
    return(-b-math.sqrt(disc))/(2*a)if b>=0 else 2*c/(-b+math.sqrt(disc))
def escape_samples(lo,a,b,c):
    cuts=[lo]
    if abs(a)<1e-18:
        if abs(b)>1e-18 and -c/b>lo:cuts.append(-c/b)
    else:
        disc=b*b-4*a*c
        if disc>=0:
            for z in((-b-math.sqrt(disc))/(2*a),(-b+math.sqrt(disc))/(2*a)):
                if z>lo:cuts.append(z)
    cuts=sorted(set(cuts));out=[]
    for u,v in zip(cuts,cuts[1:]+[math.inf]):
        z=(u+v)/2 if math.isfinite(v)else max(2*u,u+1)
        if a*z*z+b*z+c>1e-8:out.append(z)
    return out

def run():
    d,f=[json.loads((ROOT/p).read_text())for p in FILES]
    for doc in(d,f):
        for p,h in doc['sources'].items():assert R.sha(p)==h,p
    A.RAD=A.A.Q(*d['field']['new_positive_radical_square']);exactG=D.dm(d['gram'])
    labels={}
    for i,j in combinations(range(len(d['normals'])),2):
        typ,m,_=R.relation(exactG[i][j],exactG[i][i],exactG[j][j])
        if typ in('right','finite'):labels[i,j]=m
    Q=mp.matrix([[scalar(x)for x in r]for r in d['form']]);N=[mp.matrix([scalar(x)for x in r])for r in d['normals']]
    E,H=[mp.matrix([scalar(x)for x in d[key]['point']])if key=='strict_escape'else mp.matrix([scalar(x)for x in d[key]])for key in('strict_escape','center')]
    scales=[mp.sqrt(2/scalar(x))for x in d['norms']]
    G=np.array([[float(scalar(x)*scales[i]*scales[j])for j,x in enumerate(row)]for i,row in enumerate(d['gram'])])
    NE=np.array([float((n.T*Q*E)[0]*scales[i])for i,n in enumerate(N)])
    NH=np.array([float((n.T*Q*H)[0]*scales[i])for i,n in enumerate(N)])
    fvals=-2*np.cos(np.pi/np.arange(2,257));counts=Counter();found={};stats=[]
    for row in f['faces']:
        S=row['support'];outside=[i for i in range(len(d['normals']))if i not in S]
        ll={(i,j):labels[tuple(sorted((S[i],S[j])))]for i,j in combinations(range(8),2)}
        pd=P.profiles(list(range(8)),ll);J=np.array([[qval(z)for z in r]for r in pd['effective_projection_inverse']])
        P0=mp.matrix([scalar(x)for x in row['start_point']]);V0=mp.matrix([scalar(x)for x in row['direction']])
        pq,dq=scalar(row['norm_polynomial'][0]),scalar(row['norm_polynomial'][2]);U=P0/mp.sqrt(-pq);V=V0/mp.sqrt(dq)
        un=np.array([float((U.T*Q*n)[0]*scales[i])for i,n in enumerate(N)]);vn=np.array([float((V.T*Q*n)[0]*scales[i])for i,n in enumerate(N)])
        ue,ve,uh,vh=[float((v.T*Q*w)[0])for v,w in((U,E),(V,E),(U,H),(V,H))]
        profiles=[r['contacts']for r in pd['finite_profiles']]+[[[r['wall'],m]]for r in pd['dihedral_tail_profiles']for m in TAILS]
        before=len(found);accepted=cutters=0
        for contacts in profiles:
            c=np.zeros(8)
            for i,m in contacts:c[i]=-2*math.cos(math.pi/m)
            rho=2-c@J@c;assert rho>0;sr=math.sqrt(rho)/2;b=c@J@G[S];pe=c@J@NE[S];ph=c@J@NH[S]
            for orientation in(1,-1):
                a=sr*(un+orientation*vn);cc=sr*(orientation*vn-un)
                assert all(a[i]<0 and cc[i]>0 for i in outside)
                ah,ch=sr*(uh+orientation*vh),sr*(orientation*vh-uh);assert ah<0 and ch>0
                ae,ce=sr*(ue+orientation*ve),sr*(orientation*ve-ue)
                def record(z,recipe):
                    nonlocal accepted,cutters
                    if not math.isfinite(z)or z<=1e-10:return
                    hp=ah*z+ph+ch/z
                    if hp>=-1e-8:return
                    vals=a*z+b+cc/z;err=np.min(abs(vals[:,None]-fvals[None,:]),axis=1)
                    if not np.all((vals< -2-1e-7)|(err<1e-8)):return
                    t=orientation*float(mp.sqrt(-pq/dq))*(z*z-1)/(z*z+1)
                    en=float(pq)+t*t*float(dq)
                    if en>=-1e-8:counts['near_null_skips']+=1;return
                    ev=pe+ae*z+ce/z;accepted+=1;cutters+=ev>1e-8
                    rec=dict(face_id=row['id'],support=S,contacts=[[S[i],m]for i,m in contacts],orientation=orientation,
                             exponential_parameter=z,**recipe,normalized_wall_pairings=vals.tolist(),
                             center_pairing=float(hp),escape_pairing=float(ev),cuts_escape=bool(ev>1e-8),
                             endpoint_parameter=t,endpoint_norm=en,
                             finite_neighbors=[i for i in range(len(d['normals']))if vals[i]>-2+1e-7])
                    key=tuple(np.round(vals,8))
                    if key not in found:found[key]=rec
                threshold=max(positive_root(a[i],b[i]+2,cc[i])for i in outside)
                threshold=max(threshold,positive_root(ah,ph,ch));threshold+=1e-7*max(1,threshold)
                record(max(2*threshold,threshold+1),dict(method='all_dotted_sample'))
                for z in escape_samples(threshold,ae,pe,ce):record(z,dict(method='all_dotted_escape_interval'))
                for anchor in outside:
                    for m in ANCHORS:
                        z=positive_root(a[anchor],b[anchor]+2*math.cos(math.pi/m),cc[anchor]);counts['anchor_events']+=1
                        record(z,dict(method='finite_anchor',anchor=anchor,anchor_label=m))
        counts['sampled_profiles']+=len(profiles)
        stats.append(dict(face_id=row['id'],profiles=len(profiles),accepted_events=accepted,new_caps=len(found)-before,cutter_events=int(cutters)))
        print(stats[-1],'total',len(found),'cutters',sum(r['cuts_escape']for r in found.values()),flush=True)
    bank=sorted(found.values(),key=lambda r:(not r['cuts_escape'],-r['escape_pairing'],r['face_id']))
    paths=FILES+[f'scripts/{TAG}.py','scripts/fresh_h9_polar_endpoint_profiles_20261015.py',
                 'scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014.py','scripts/fresh_h9_quartet_polar_second_20261014.py']
    return dict(status='NUMERICAL_COMPLETE_LINE_ENDPOINT_PROPOSALS',sources={p:R.sha(p)for p in paths},
                counts=dict(counts),anchor_labels=ANCHORS,tail_labels=TAILS,faces=stats,
                candidate_count=len(bank),escape_cutter_count=sum(r['cuts_escape']for r in bank),candidates=bank,
                scope='Bounded numerical Coxeter-event proposals on allremaining complete geodesic one-faces, both cap orientations and unrestricted real line parameter. Local profiles are complete except sampled I2 tails; outside-anchor labels and finite-angle gate remain bounded. Includesall-dotted interval escape proposals. No exact cap or all-label exclusion is claimed.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();bank=d.pop('candidates')
    blob=gzip.compress(json.dumps(bank,separators=(',',':')).encode(),mtime=0);bp=ROOT/f'data/enumerations/{TAG}_bank.json.gz'
    d['candidate_bank']=dict(path=str(bp.relative_to(ROOT)),sha256=sha256(blob).hexdigest(),count=len(bank));d['best_escape_pairing_samples']=bank[:24]
    out=ROOT/f'results/{TAG}.json'
    if args.verify:assert blob==bp.read_bytes()and d==json.loads(out.read_text())
    else:bp.write_bytes(blob);out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['candidate_count'],d['escape_cutter_count'])
