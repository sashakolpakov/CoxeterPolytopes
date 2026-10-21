#!/usr/bin/env python3
"""Move a cavity duty from old30 to an actual cap X; bounded numerical scout.

The new frame(X,67,60) is absent from the frozen299-equation search.
Every accepted addition is scored on the complete S6 cone and by independent
PD-support incidence. No saved coordinate is promoted to an exact root.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import numpy as np
import fresh_release58_66_seven_boundary_20261017 as N
ROOT,M,G=N.ROOT,N.M,N.G
TAG='fresh_release58_66_shifted_anchor_patch_20261017'


def run():
    src=ROOT/'results'/f'{N.TAG}.json';parent=json.loads(src.read_text());s=N.Scout()
    anchor=parent['anchor_domain']['caps'][parent['chosen_anchor_index']]
    assert anchor['profile']==0 and anchor['labels']==[3,2,2]
    nd=dict(s.base);nd[67]=np.array(anchor['normal']);types={67:0};actual=s.actual+[67]
    v,bad,cp=s.native(nd[67],0,s.base,{})
    assert not bad and cp<-1e-8
    assert max(abs(v[i]+N.B.theta(m))for i,m in zip(anchor['frame'],anchor['labels']))<1e-8
    sid=s.score(nd,types,actual,'anchor');assert s.ledgers[sid]['open_ends']==11
    xs=parent['wedge_banks'][1]['caps'];ys=parent['wedge_banks'][4]['caps'];assert len(xs)==39 and len(ys)==8
    events=[];joins=[];joinhist=Counter();compact_anchor=0
    for x in xs:
        nx,tx=s.add(nd,types,[x]);vx,badx,cx=s.native(nx[68],tx[68],nd,types);assert not badx and cx<-1e-8
        assert max(abs(vx[i]+N.B.theta(m))for i,m in zip(x['frame'],x['labels']))<2e-7
        assert s.corner([33,30,68],nx,tx)and s.corner([30,67,68],nx,tx)
        xid=s.score(nx,tx,actual+[68],'X_prefix');assert xid is not None
        edges={tuple(sorted(e['sides'])):e['kind']for e in s.ledgers[xid]['edges']}
        assert edges[(67,68)]=='one_open'and edges[(60,67)]=='one_open'
        bank=s.bank([68,67,60],[(0,1),(1,2)],nx,tx,actual+[68],69,'shifted_Z',score=True)
        for z in bank['caps']:
            nz,tz=s.add(nx,tx,[z],69)
            assert s.corner([68,67,69],nz,tz)and s.corner([67,60,69],nz,tz)
            ledger=s.ledgers[z['score']];ar=[r for r in ledger['rays']if 67 in r['active']]
            z['anchor_face_rays']=ar;z['whole_anchor_compact']=bool(ar)and all(r['kind']=='finite'for r in ar)
            compact_anchor+=int(z['whole_anchor_compact'])
            for y in ys:
                nn,tt=s.add(nz,tz,[y],70)
                vals,bad,cy=s.native(nn[70],tt[70],nz,tz)
                relation=G.relation(vals[69]);joinhist[str(relation)]+=1
                if bad or cy>=-1e-8:continue
                sid4=s.score(nn,tt,actual+[68,69,70],'X_Z_Y')
                if sid4 is None:continue
                rec=dict(X=x['index'],Z=z['index'],Y=y['index'],YZ_product=vals[69],YZ_relation=relation,score=sid4)
                joins.append(rec)
                if s.ledgers[sid4]['open_ends']<10:s.improvements.append(dict(stage='X_Z_Y',record=rec,roots=[x,z,y]))
        events.append(dict(X=x,new_boundary=[7,33,68,67,60,59,55],X_score=xid,bank=bank))
        print('SHIFTED_ANCHOR',x['index'],'Xprofile',x['profile'],'recipes',len(bank['symbolic_recipes']),'Z',len(bank['caps']),'joins',len(joins),'best',min(l['open_ends']for l in s.ledgers),'unresolved',len(s.unresolved),flush=True)
    paths=[Path(__file__),Path(N.__file__),src,Path(N.H.__file__),Path(N.E.__file__),Path(N.N.__file__),Path(N.S.__file__),Path(N.P.__file__),Path(N.G.__file__),N.H.PARENT,N.H.BOUNDARY]
    return dict(status='BOUNDED_SHIFTED_ANCHOR67_PATCH',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      exact_deletion_ledger=s.deletion,retained_old_walls=61,anchor=anchor,anchor_score=sid,old_control_score=s.control,deletion_score=s.deleted,
      original_boundary=N.CYCLE,new_boundary_template=[7,33,68,67,60,59,55],new_frame=[68,67,60],required_corner_pairs=[[0,1],[1,2]],
      X_prefixes=len(xs),Y_companions=len(ys),high_label_samples=N.HIGH,events=events,coherent_Y_joins=joins,YZ_relation_histogram=dict(joinhist),
      compact_anchor_faces=compact_anchor,counts=dict(s.counts),boundary_ledgers=s.ledgers,improvements=s.improvements,unresolved=s.unresolved,
      histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n)for(a,b,c,d),n in sorted(s.hist.items())],
      blockers=[dict(stage=a,walls=list(b),count=n)for(a,b),n in sorted(s.blockers.items())],near_native_failures=s.near,
      scope='Changed incidence: every one of39 saved native X caps at33/30/67 replaces old30 on the actual seven-edge cavity boundary, then Z is solved in X/67/60 rather than old30/67/60. All seven A2 profiles and every finite low2..5 recipe; symbolic high slots sampled only6/8/12/24. All61 retained old native contacts, all prior/new products and strict common center are checked. Every accepted cap/compatible Y join has complete numerical S6cone/PD-support agreement. Exact old58/66 deletion debit is inherited; no exact new normal, ambient finite census, all-label exclusion or finite-volume claim. The exact63-wall199-finite10-open-end control remains separate.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'improvements',len(out['improvements']),'unresolved',len(out['unresolved']),flush=True)
