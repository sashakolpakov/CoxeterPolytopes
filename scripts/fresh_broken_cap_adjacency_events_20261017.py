#!/usr/bin/env python3
"""Directed continuous R44-contact events for the broken-adjacency scout.

Solve X's two observed blocking old contacts first. R44 is then an output
of a norm equation, rather than a prescribed sampled finite label.
"""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_broken_cap_adjacency_20261017 as T

S,H,B,G,M,ROOT=T.S,T.H,T.B,T.G,T.M,T.ROOT
TAG='fresh_broken_cap_adjacency_events_20261017'

def run():
    z,labels,ns,h,profiles,primary=H.setup();actual=[i for i in z['actual_sides']if i not in S.DELETE]
    oldsq={r['wall']:B.mid(r['square'])for r in z['projections']}
    k0,k1=T.a_curve(ns,h);xbank=[];rejected_x=Counter();counts=Counter();unresolved=[]
    for ms in product(range(2,6),repeat=3):
        roots,norm=T.N.norm_roots([ns[i]for i in[8,30,63]],[-B.theta(m)for m in ms])
        assert not norm['status'].startswith('unresolved')
        for sg,xy in roots:
            vals=T.native_values(xy,0,labels,ns,profiles);bad=T.gates(vals,[8,30,63]);center=float(xy@M@h)
            if center>=-1e-8:counts['X_wrong_center']+=1;continue
            if bad:rejected_x[tuple(i for i,v in bad)]+=1;continue
            xbank.append(dict(index=len(xbank),old_frame=[8,30,63],old_labels=ms,branch=sg,projection=xy.tolist(),center_product=center))
    assert len(xbank)==57
    # Each entry represents a complete finite shape; a symbolic old8 or63
    # high label is sampled explicitly. The old44 contact remains unknown.
    rdomains=defaultdict(set);rd=S.domain([8,44,63])
    for rec in rd['low']+rd['tails']:
        for m8,m63 in product(T.HIGH_LABELS if rec['labels'][0]==6 else[rec['labels'][0]],
                             T.HIGH_LABELS if rec['labels'][2]==6 else[rec['labels'][2]]):
            rdomains[rec['profile'],m8,m63].add(rec['labels'][1])
    blocks=Counter();near=[];rbank=[];bridges=[]
    for xx in xbank:
        xy=np.array(xx['projection']);x8=xx['old_labels'][0]
        for(rpi,r8,r63),possible44 in sorted(rdomains.items()):
            for rx in [2,3,4,5]+T.HIGH_LABELS:
                if T.corner(rpi,r8,rx,0,x8)is None:continue
                roots,norm=T.N.norm_roots([ns[8],ns[63],xy],[-B.theta(r8),-B.theta(r63),-B.theta(rx)],profiles[rpi]['square'])
                counts['R_event_norm_equations']+=1
                if norm['status'].startswith('unresolved'):
                    unresolved.append(dict(X=xx['index'],R_profile=rpi,R8=r8,R63=r63,RX=rx,norm=norm));continue
                for sg,ry in roots:
                    counts['R_event_branches']+=1
                    assert abs(ry@M@ry-profiles[rpi]['square'])<2e-6
                    center=float(ry@M@h+profiles[rpi]['center_offset'])
                    if center>=-1e-8:counts['R_wrong_center']+=1;continue
                    vals=T.native_values(ry,rpi,labels,ns,profiles);bad=T.gates(vals,[8,63]);rel44=B.relation(vals[44]);counts['future_R_events']+=1
                    blocks[tuple(i for i,v in bad)]+=1
                    record=dict(X=xx['index'],profile=rpi,old8=r8,old63=r63,RX=rx,branch=sg,
                                forced_old44_product=vals[44],forced_old44_relation=rel44,possible44_shapes=sorted(possible44),
                                projection=ry.tolist(),center_product=center,blockers=bad)
                    if len(bad)<=1:near.append(record)
                    if bad:continue
                    if not isinstance(rel44,int):counts['R44_not_finite']+=1;continue
                    if (rel44 if rel44<6 else 6)not in possible44:counts['R44_wrong_finite_shape']+=1;continue
                    assert S.finite_corner(8,44,rpi,r8,rel44)is not None and S.finite_corner(44,63,rpi,rel44,r63)is not None
                    record['index']=len(rbank);rbank.append(record);counts['coherent_R_events']+=1
                    for ax in [2,3,4,5]+T.HIGH_LABELS:
                        if T.corner(0,x8,ax,0,2)is None:continue
                        aa=float(np.sqrt(.5)*(k0+k1)@M@xy);bb=B.theta(ax);cc=float(np.sqrt(.5)*(-k0+k1)@M@xy)
                        if abs(aa)<1e-10:
                            unresolved.append(dict(X=xx['index'],R=record['index'],AX=ax,reason='degenerate_A_parameter_equation'));continue
                        for root in np.roots([aa,bb,cc]):
                            if abs(root.imag)>1e-8 or root.real<=0:continue
                            t=float(root.real);ay=np.sqrt(.5)*(t*(k0+k1)+(-k0+k1)/t)
                            aval=T.native_values(ay,0,labels,ns,profiles);abad=T.gates(aval,[7,8]);ac=float(ay@M@h)
                            if abad or ac>=-1e-8:continue
                            ar=float(ay@M@ry);arl=B.relation(ar)
                            if arl is None or arl=='parallel_or_unresolved' or ar>=-np.sqrt(2*profiles[rpi]['square'])-1e-7:continue
                            nd=dict(ns);nd.update({67:ay,68:ry,69:xy})
                            bd=S.C.detailed_cone(nd,actual+[67,68,69])
                            inc=H.incidence_score([nd[i]for i in actual+[67,68,69]],actual+[67,68,69],
                                [oldsq[i]for i in actual]+[2,profiles[rpi]['square'],2])
                            bridges.append(dict(X=xx['index'],R=record['index'],A_t=t,AX=ax,A_projection=ay.tolist(),AR_product=ar,
                                boundary=bd,PD_incidence=inc))
                            counts['coherent_ARX_event_triples']+=1
        if(xx['index']+1)%10==0:print('BROKEN_EVENT_X',xx['index']+1,'of57',dict(counts),'bridges',len(bridges),flush=True)
    old=json.loads((ROOT/'results'/f'{T.TAG}.json').read_text())
    paths=[Path(__file__),Path(T.__file__),Path(S.__file__),Path(H.__file__),ROOT/'results'/f'{T.TAG}.json',H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_NUMERICAL_BROKEN_ADJACENCY_CONTINUOUS_R44_EVENTS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        deleted_walls=S.DELETE,boundary_cycle=S.CYCLE,exact_old_deletion_cost=old['exact_old_deletion_cost'],
        counts=dict(counts,X_low_old_contact_recipes=64,X_native_coherent_roots=57),
        X_bank=xbank,X_rejections=[[list(k),v]for k,v in sorted(rejected_x.items())],
        complete_R_finite_shapes=rd,R_event_blockers=[[list(k),v]for k,v in sorted(blocks.items())],
        near_R_events=near,coherent_R_events=rbank,coherent_ARX_triples=bridges,unresolved_events=unresolved,
        equation='X is solved from old8/30/63 contacts and norm2. For each required finite R endpoint shape and R-X contact, solve R in frame[old8,old63,X]; its old44 product is determined continuously by the norm branches, then checked against the full native Coxeter gates. Any surviving A parameter solves a*t^2+theta(AX)*t+c=0 exactly at the numerical equation level.',
        scope='A bounded numerical contact-event continuation driven by the two old blockers30/63 of the first scout. The R44 contact is not prescribed from a finite label grid; its derived value must pass the native gate. X old labels are2..5 and other symbolic high labels are sampled6,10,24, so no all-label or all-parameter obstruction is asserted. Any accepted triples still need B/C completion, permanent-cusp screening, exact physical promotion, and complete ambient volume analysis. No new exact chamber.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=json.loads(json.dumps(run()));dest=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert out==json.loads(dest.read_text())
    else:dest.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],flush=True)
