#!/usr/bin/env python3
"""Numerical intersections of two cap-orbit norm curves in seed parameters."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import brentq
import networkx as nx

import fresh_h9_multiblock_seed_20261013 as base
from fresh_h9_asymmetric_blocks_20261013 import reduced_geometry

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/enumerations/fresh_h9_asymmetric_intersections_20261013_scout.json'


def coeffs(ns,b,P):
    _,_,tv,dv,sq,_=P;n=np.array(ns);den=4-b*n
    C=np.sum(tv*tv,axis=1)/(2*b)+sq/8+b*np.sum(dv*dv/den,axis=1)/8-2
    A=np.sum(tv,axis=1)**2/(2*b)
    B=np.sum(dv/den,axis=1)**2/2;L=np.sum(n/den)
    return C*2*L-A*L-B*2,C*(2-b*L)-A+B*b,-C*b


def curves(ns,b,P,indices):
    aa,bb,cc=(z[indices]for z in coeffs(ns,b,P))
    d=bb*bb-4*aa*cc
    with np.errstate(invalid='ignore',divide='ignore'):
        roots=np.stack([(-bb-np.sqrt(d))/(2*aa),(-bb+np.sqrt(d))/(2*aa)],axis=1)
    roots.sort(axis=1)
    roots[(roots<=b+2+1e-6)|(roots>=20)]=np.nan
    return roots.ravel()


def run(output,grid,groups):
    started=time.monotonic();records=[];counts=[]
    for ns in groups:
        P=base.profiles(ns)
        # Deduplicate identical curve functions using samples. Discovery
        # heuristic only; exact equality is not claimed from these samples.
        signatures=np.column_stack([z for b in (.041,.113,.217)for z in coeffs(ns,b,P)])
        signatures/=np.maximum(1,np.linalg.norm(signatures,axis=1))[:,None]
        _,indices=np.unique(np.round(signatures,9),axis=0,return_index=True)
        bs=np.linspace(.005,4/max(ns)-.005,grid)
        old=None;previous=None;points={};brackets=0
        ii,jj=np.triu_indices(2*len(indices),1)
        ii,jj=ii[ii//2!=jj//2],jj[ii//2!=jj//2]
        for b in bs:
            roots=curves(ns,b,P,indices);diff=roots[ii]-roots[jj]
            if old is not None:
                crossed=np.flatnonzero(np.isfinite(diff)&np.isfinite(old)&(diff*old<0)&(abs(diff)>1e-9)&(abs(old)>1e-9))
                for t in crossed:
                    u,v=int(ii[t]),int(jj[t]);brackets+=1
                    def fun(x):
                        z=curves(ns,x,P,indices);return z[u]-z[v]
                    try:b0=brentq(fun,previous,b,xtol=1e-13)
                    except ValueError:continue
                    r=curves(ns,b0,P,indices);a0=float((r[u]+r[v])/2)
                    if not np.isfinite(a0)or abs(r[u]-r[v])>1e-6:continue
                    points.setdefault((round(a0,8),round(b0,8)),
                                      (a0,float(b0),[int(indices[u//2]),int(indices[v//2])]))
            old=diff;previous=b
        ranked=[]
        for a,b,identities in points.values():
            qa,qb,qc=coeffs(ns,b,P);res=qa*a*a+qb*a+qc
            scale=abs(qa*a*a)+abs(qb*a)+abs(qc)+1
            selected=np.flatnonzero(abs(res)/scale<1e-10).tolist()
            ranked.append((sum(P[-1][selected]),a,b,selected,identities))
        ranked.sort(reverse=True)
        count={'groups':ns,'distinct_sampled_curves':len(indices),'grid':grid,
               'crossing_brackets':brackets,'parameter_intersections':len(points),
               'top_cap_count':int(ranked[0][0])if ranked else 0}
        counts.append(count);print(count,flush=True)
        for _,a,b,profiles,identities in ranked:
            G,R,refs=base.bank(ns,a,b,P,profiles)
            geo=reduced_geometry(ns,G,R)
            geo.pop('N')
            record={'groups':ns,'a':a,'b':b,'profile_indices':profiles,
                    'crossing_profiles':identities,'cap_count':len(R),'geometry_all_caps':geo}
            if geo['outside']==0:
                C=base.compat(R@G@R.T);record['compatible_pair_count']=int(C.sum()//2)
                record['cap_covectors']=refs
            records.append(record)
        output.write_text(json.dumps({'status':'NUMERICAL_CAP_CURVE_INTERSECTIONS',
            'finite_pairing_magnitudes':base.VALS.tolist(),
            'census':counts,'records':records,'seconds':round(time.monotonic()-started,3),
            'scope':'Numerical sign-change sampling of intersections of two cap-orbit norm curves for listed seed partitions, b in (0,4/max(n)) and a in (b+2,20), with the recorded finite pairing alphabet. Tangencies and unsampled intersections may be missed. Complete prescribed finite pairing orbit banks at each proposed parameter pair; numerical ordered-sector geometry of their intersection, without mutual compatibility requirement. No exact exclusion or finite-volume claim.'},indent=2)+'\n')
        print({'completed_groups':ns,'records':len(records),
               'best_max_klein_square':min(r['geometry_all_caps']['max_klein_square']for r in records),
               'seconds':round(time.monotonic()-started,2)},flush=True)
    print({'done':len(records),'no_outside':[i for i,r in enumerate(records)if r['geometry_all_caps']['outside']==0]},flush=True)


def select(source):
    data=json.loads(source.read_text());cases=[]
    base.VALS=np.array(data.get('finite_pairing_magnitudes',[0,1,base.PHI]))
    for index,row in enumerate(data['records']):
        if row['geometry_all_caps']['outside']:continue
        P=base.profiles(row['groups'])
        G,R,_=base.bank(row['groups'],row['a'],row['b'],P,row['profile_indices'])
        C=base.compat(R@G@R.T)
        cliques=list(nx.find_cliques(nx.from_numpy_array(C)));checks=[]
        for clique in cliques:
            geometry=base.geometry(G,R,clique)
            checks.append({'clique':clique,'geometry':geometry})
        record={'source_record':index,'cap_count':len(R),
                'compatibility_edges':int(C.sum()//2),
                'maximal_cliques':len(cliques),'maximum_clique_size':max(map(len,cliques)),
                'best_max_klein_square':min(c['geometry']['max_klein_square']for c in checks),
                'all_maximal_cliques_have_outside_witness':all(c['geometry']['outside']for c in checks),
                'checks':checks}
        cases.append(record)
        print({k:v for k,v in record.items()if k!='checks'},flush=True)
    source.with_name(source.stem.replace('_scout','_selection_scout')+'.json').write_text(json.dumps({
        'status':'NUMERICAL_MAXIMAL_CLIQUE_GEOMETRY_SCREEN','cases':cases,
        'scope':'All maximal cliques of the numerical compatibility graphs using finite labels2..120 or dotted pairs for every sampled all-bank bounded candidate. Full Klein geometry is recomputed without symmetry restriction for each clique. Numerical witnesses only; no exact algebraic signs or exhaustive parameter claim.'},indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=OUT)
    p.add_argument('--grid',type=int,default=500)
    p.add_argument('--select',action='store_true')
    p.add_argument('--labels',type=int,nargs='+',default=[2,3,5])
    p.add_argument('--groups',type=int,nargs='+',default=[4,5,6,7],help='First block size; second is8 minus it')
    args=p.parse_args()
    if args.select:select(args.output)
    else:
        base.VALS=np.array([0. if m==2 else 2*np.cos(np.pi/m)for m in args.labels])
        run(args.output,args.grid,[[n,8-n]for n in args.groups])
