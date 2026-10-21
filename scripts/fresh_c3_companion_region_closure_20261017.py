#!/usr/bin/env python3
"""Whole five-edge B-region cap at the coherent C3 companion triples.

The necessary finite corner domain has no unbounded label: each contact in
the three-frame shares a required finite corner with a retained bond >=4.
Prefix coordinates and root/cone acceptance remain numerical.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_c3_companion_continuation_20261017 as E
C,J,H,B,M,N,S,P,ROOT=E.C,E.J,E.H,E.B,E.M,E.N,E.S,E.P,E.ROOT
TAG='fresh_c3_companion_region_closure_20261017'

def run():
    pp=ROOT/'results'/f'{C.TAG}.json';cp=ROOT/'results'/f'{E.TAG}.json'
    parents=json.loads(pp.read_text());prior=json.loads(cp.read_text())
    cur,labels,base,h,profiles,primary=C.setup();counts=Counter();blockers=Counter();rows=[];hits=[];unresolved=[];scores=Counter();improvements=[]
    cycle=[30,52,70,68,69]
    def off(a,b):return float(H.A.native(np.zeros(4),a,np.zeros(4),b,profiles))
    for record in prior['records']:
        pr=parents['companions'][record['parent']]
        for pair in record['coherent_pairs']:
            a,b=[record['continuations'][i]for i in pair['indices']]
            assert a['frame']==[30,67,68]and b['frame']==[47,52,68]
            nd=dict(base);nd.update({68:np.array(pr['normal']),69:np.array(a['normal']),70:np.array(b['normal'])})
            ps={68:pr['profile'],69:a['profile'],70:b['profile']};ps=lambda i,ps=ps:ps.get(i,0)
            desired_edges={tuple(sorted(t))for t in zip(cycle,cycle[1:]+cycle[:1])}
            assert desired_edges<={tuple(e['sides'])for e in pair['boundary']['edges']if e['kind']=='one_open'}
            m3052=B.relation(float(nd[30]@M@nd[52]));m5270=b['labels'][1]
            assert m3052==4 and m5270>=4
            domain=[];native=[];local=Counter()
            for pi,ms in product(range(7),product(range(2,7),repeat=3)):
                if P.corner3(0,0,pi,m3052,ms[0],ms[1])is None:continue
                if P.corner3(0,ps(70),pi,m5270,ms[1],ms[2])is None:continue
                assert 6 not in ms
                domain.append([pi,*ms]);rhs=[-B.theta(m)-off(pi,ps(i))for i,m in zip([30,52,70],ms)]
                roots,info=N.norm_roots([nd[i]for i in[30,52,70]],rhs,profiles[pi]['square']);counts['norm_equations']+=1
                if info['status'].startswith('unresolved'):
                    unresolved.append(dict(parent=pr['index'],pair=pair['indices'],profile=pi,labels=ms,norm=info));continue
                for sg,v in roots:
                    counts['norm_branches']+=1
                    err=max(abs(float(v@M@v)-profiles[pi]['square']),max(abs(float(v@M@nd[i])-r)for i,r in zip([30,52,70],rhs)))
                    if err>1e-5:unresolved.append(dict(parent=pr['index'],pair=pair['indices'],residual=err));continue
                    vals={i:float(v@M@nd[i])+(profiles[pi]['offsets'][labels.index(i)]if i in labels else off(pi,ps(i)))for i in nd}
                    bad=[i for i,x in vals.items()if B.relation(x)in[None,'parallel_or_unresolved']]
                    if bad:blockers[tuple(bad)]+=1;local['native_blocked']+=1;continue
                    center=float(v@M@h)+profiles[pi]['center_offset']
                    if center>=-1e-8:local['center_rejected']+=1;continue
                    rel={i:B.relation(vals[i])for i in cycle}
                    finite=[]
                    for i,j in zip(cycle,cycle[1:]+cycle[:1]):
                        mij=B.relation(float(nd[i]@M@nd[j])+off(ps(i),ps(j)))
                        assert isinstance(mij,int)
                        finite.append(isinstance(rel[i],int)and isinstance(rel[j],int)and P.corner3(ps(i),ps(j),pi,mij,rel[i],rel[j])is not None)
                    rec=dict(profile=pi,labels=list(ms),branch=sg,normal=v.tolist(),center_product=center,native_products=sorted(vals.items()),five_finite_corners=finite)
                    native.append(rec);counts['native_candidates']+=1
                    if not all(finite):local['missing_required_corner']+=1
                    nn=dict(nd);nn[71]=v;sel=cur['actual_sides']+[67,68,69,70,71]
                    bd=S.C.detailed_cone(nn,sel);ic=H.incidence_score([nn[i]for i in sel],sel,[float(nn[i]@M@nn[i])for i in sel])
                    good=not bd['warnings']and ic['status']=='numerically_separated_acute_incidence'and bd['open_ends']==ic['open_ends']and bd['counts'].get('finite',0)==ic['finite']
                    good=good and {tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}=={tuple(sorted(r))for r in ic['finite_supports']}
                    if not good:unresolved.append(dict(parent=pr['index'],pair=pair['indices'],stage='boundary',boundary=bd,PD_incidence=ic));continue
                    fr=[r for r in bd['rays']if 71 in r['active']]
                    rec.update(boundary=bd,PD_incidence=ic,complete_new_facet=fr)
                    key=(bd['counts']['finite'],bd['counts']['exterior'],bd['open_ends']);scores[key]+=1
                    if bd['open_ends']<10:improvements.append(dict(parent=pr['index'],pair=pair['indices'],cap=rec))
                    if all(finite):
                        assert len(fr)==5 and all(r['kind']=='finite'for r in fr)
                        hits.append(dict(parent=pr['index'],pair=pair['indices'],cap=rec));counts['whole_region_caps']+=1
            rows.append(dict(parent=pr['index'],pair=pair['indices'],first_companion_labels=pr['labels'],L_labels=a['labels'],R_labels=b['labels'],R_profile=b['profile'],complete_finite_domain=domain,native_candidates=native,counts=dict(local)))
    paths=[Path(__file__),Path(C.__file__),Path(E.__file__),pp,cp,H.PARENT,H.BOUNDARY]
    return dict(status='NUMERICAL_C3_WHOLE_REGION_CAP_DOMAIN',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      deleted_walls=[],retained_control_walls=63,boundary_cycle=cycle,prefixes=len(rows),counts=dict(counts),records=rows,hits=hits,improvements=improvements,unresolved=unresolved,
      complete_boundary_histogram=[dict(finite=f,exterior=e,open_ends=i,count=n)for(f,e,i),n in sorted(scores.items())],
      blockers=[dict(walls=list(k),count=n)for k,n in sorted(blockers.items())],
      label_domain_proof='A norm2 acute finite triple with one bond>=6 and another nonright bond has determinant<=0. Every assigned cap contact in30/52/R shares a required finite corner with the fixed30-52 label4 or52-R label>=4. Thus no assigned cap label>=6 is possible, including every nonzeroA2profile; the complete finite forest classifier then tests labels2..5. Numerical root realization at the100 fixed prefix triples is not an exact exclusion.',
      scope='A fresh whole-B-region closure duty with no deletions. The100 coherent numerical prefixes have five open boundaryedges30-52-R-A-L-30. A successful cap must give allfive finite corners and an entire compact pentagonal facet, followed by complete target cone/PD agreement. This source exhausts the finite necessary frame-domain at these literal numerical prefixes, but neither exactifies their coordinates nor excludes deformed prefixes. No finite-volume claim without an exact full ambient certificate.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['prefixes'],out['counts'],'hits',len(out['hits']),'unresolved',len(out['unresolved']),flush=True)
