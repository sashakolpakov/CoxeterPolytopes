#!/usr/bin/env python3
"""Freeze a norm-coupled reflection pair and its exactly checked escape."""
import json
from fractions import Fraction as Q
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_crossplane_repair_20261013_pair_verify as P
V=P.V;K=P.K

def run():
    source=P.B.PREFIX+'_caps.json'
    d={'status':'EXACT_ORTHOGONAL_PAIR_REPAIR','source':source,'source_sha256':V.sha(source),
       'constraints':P.CONSTRAINTS,'difference_pairing':2,
       'targets':[V.read(source)['escape'],V.read(P.T.PREFIX+'_caps.json')['escape']],
       'definition':{'D':'seed6-seed7, norm4','A':'The eight listed old normals followed by D.',
                     'u':'sum_i (Gram(A)^-1 p)_i A_i, p ending in2',
                     'v':'common center minus its projection onto span(A)',
                     'r':'u + branch*sqrt((2-u²)/v²)*v','s':'r-D; norm2 and r.s=0 exactly'},
       'scope':'Bounded adaptive metric-equation search, with a prescribed orthogonal reflection coupling. No exhaustive family claim.'}
    for branch in(-1,1):
        d['branch']=branch
        try:out=P.check(d)
        except AssertionError:continue
        break
    else:raise AssertionError('Neither norm branch passes the exact pair gates')
    G,I,h,g,old,names,caps,center,_=P.reconstruct(d)
    scout=V.read(P.T.PREFIX+'_geometry_scout.json');case=next(c for c in scout['cases']if c['proposal']==0 and c['caps']==2)
    f=lambda p:sum(p)/(2*V.SCALE)
    hn=[f(V.golden(x))for x in h];scale=sum(hn);hn=[x/scale for x in hn]
    for eps in(1e-4,1e-6,1e-8):
        y=tuple((Q(round((x+eps*c)*10**12),10**12),Q(0))for x,c in zip(case['new_escape'],hn,strict=True))
        if K.sign(K.pair(y,I,y))<=0 or K.sign(K.pair(y,I,h))>=0:continue
        pp=[V.dot(r,list(map(V.golden,y)))for r in old+caps]
        if min(p[0]for p in pp)>0:break
    else:raise AssertionError('Scout escape did not exactify')
    d['escape']=list(map(K.enc,y));d['numerical_geometry']={k:v for k,v in case.items()if k not in('normals','old_escape','new_escape')}
    d['numerical_geometry']['source']=P.T.PREFIX+'_geometry_scout.json'
    d['numerical_geometry']['source_sha256']=V.sha(P.T.PREFIX+'_geometry_scout.json')
    out=P.check(d);path=P.PREFIX+'_caps.json';(V.ROOT/path).write_text(json.dumps(d,indent=2)+'\n');out['source_sha256']=V.sha(path)
    (V.ROOT/(P.PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'], 'EXACT_ESCAPE',out['affine_and_vertex_census'])

if __name__=='__main__':run()
