#!/usr/bin/env python3
"""Complete local first-edge probes using shifted integral pairing cosets."""
import argparse,hashlib,json,time
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path
from fresh_h9_new_host_20261013 import *
from fresh_h9_quartic_mixed_20261013_endpoint import square_root
from fresh_h9_quartic_first_endpoint_20261013 import line_key,finite_orbit,root_test
from fresh_h9_quartic_edge_walk_20261013 import finite_vertex

SOURCE=ROOT/'results/fresh_h9_new_host_20261013_corrected.json'
OUT=ROOT/'results/fresh_h9_new_host_next_20261013_census.json'


def probe(start,b,g,C):
    N=pair(start,g,start);B=pair(b,g,b);target=(-N/B)*(1-K(1)/C)
    assert target>0
    low=K(1)
    while low*low>target:low/=2
    while 4*low*low<=target:low*=2
    high=2*low
    for _ in range(28):
        mid=(low+high)/2
        if mid*mid<=target:low=mid
        else:high=mid
    p=tuple(x-low*y for x,y in zip(start,b));assert pair(p,g,p)<0
    p=tuple(x/p[0]for x in p)
    actual=pair(start,g,p)**2/(N*pair(p,g,p));assert 1<actual<=C
    return p,low,actual


def run(edge_index,Cvalue):
    started=time.monotonic();data=json.loads(SOURCE.read_text());fd=data['forced_H5_face'];g=Host(1).g;g9=Host(2).g
    walls=[tuple(map(K.from_json,r))for r in fd['all_walls']];walls9=[tuple(map(K.from_json,r))for r in data['all_walls']]
    er=data['next_existing_wall_frontier'][edge_index];ids=er['edge_indices'];released=er['released_wall'];edge=[walls[i]for i in ids];w=walls[released]
    k0=tuple(map(K.from_json,er['start']));N=pair(k0,g,k0);assert N<0
    eg=tuple(tuple(pair(r,g,s)for s in edge)for r in edge);ei=inverse(eg)
    cf=mv(ei,tuple(pair(w,g,e)for e in edge));b=tuple(w[j]-sum((c*e[j]for c,e in zip(cf,edge)),Z)for j in range(6));B=pair(b,g,b)
    assert all(B.sign(e)>0 for e in range(4))and pair(k0,g,b)==0
    k1,step,C=probe(k0,b,g,K(Cvalue));assert all(pair(k1,g,r)<=0 for r in walls)and pair(k1,g,w)<0
    profiles=[];bounds=[];span={};folded={};records=[];stats=Counter()
    output=OUT.with_name(OUT.stem+f'_e{edge_index}_c{Cvalue}.json')
    for q in NORMS:
        alph=[balanced_alphabet((q*pair(e,g,e),)*4)for e in edge];before=len(profiles)
        for fixed in product(*alph):
            cf=mv(ei,tuple(-z for z in fixed));r0=tuple(sum((c*e[j]for c,e in zip(cf,edge)),Z)for j in range(6));R=q-pair(r0,g,r0)
            if R==0:
                if root_test(r0,q,g):span[line_key(r0)]=r0
                continue
            if any(R.sign(e)<=0 for e in range(4)):continue
            offset=pair(r0,g,w);caps=(R*B*C,*(R*B for _ in range(3)))
            center=K(*(round(Fraction(x))for x in offset.to_json()));shift=center-offset
            loose=tuple(2*x+2*shift*shift for x in caps);positive=balanced_alphabet(loose)
            values=sorted({center+z for x in positive for z in(x,-x)if all((caps[e]-(center+z-offset)**2).sign(e)>=0 for e in range(4))})
            idx=len(profiles);profiles.append({'norm':q.to_json(),'fixed_magnitudes':enc(fixed),'fixed_component':enc(r0),'residual':R.to_json(),
                'offset':offset.to_json(),'caps':enc(caps),'integral_center':center.to_json(),'loose_caps':enc(loose),'pairing_values':enc(values)})
            for value in values:
                anchor=tuple(z+(value-offset)*bb/B for z,bb in zip(r0,b));assert pair(anchor,g,k0)==0 and pair(anchor,g,w)==value
                square=(q-pair(anchor,g,anchor))/N;stats['quadratics']+=1;radical=square_root(square)
                if radical is None:continue
                stats['field_squares']+=1
                for t in {radical,-radical}:
                    r=tuple(x+t*y for x,y in zip(anchor,k0));assert pair(r,g,r)==q
                    if not root_test(r,q,g):continue
                    stats['primitive_integral']+=1;pa,pb=pair(k0,g,r),pair(k1,g,r)
                    if pa*pb>0:continue
                    assert pa or pb
                    lam=-pa/(pb-pa);assert 0<=lam<=1
                    key=line_key(r)
                    if key in folded:continue
                    folded[key]=r;records.append({'profile':idx,'vector':enc(r),'norm':q.to_json(),'parameter':lam.to_json()})
        bounds.append({'norm':q.to_json(),'pairing_alphabets':[enc(v)for v in alph]})
        print({'edge':edge_index,'C':Cvalue,'norm':str(q),'profiles':len(profiles)-before,'quadratics':stats['quadratics'],'mirrors':len(records),'seconds':round(time.monotonic()-started,3)},flush=True)
    edge_group=finite_orbit(edge,edge,g);assert set(finite_orbit(list(span.values()),edge,g))==set(edge_group)
    expanded=finite_orbit(list(folded.values()),edge,g)
    start_group={key:r for key,r in expanded.items()if pair(k0,g,r)==0};start_group.update(edge_group)
    expected_start=finite_orbit(edge+[w],edge+[w],g);assert set(start_group)==set(expected_start)
    positive=[r for r in records if K.from_json(r['parameter'])>0]
    endpoint=None
    if positive:
        first=min(K.from_json(r['parameter'])for r in positive);point=tuple(x+first*(y-x)for x,y in zip(k0,k1))
        full={key:r for key,r in expanded.items()if pair(point,g,r)==0};full.update(edge_group);options=[]
        for row in positive:
            if K.from_json(row['parameter'])!=first:continue
            r=tuple(map(K.from_json,row['vector']));q=K.from_json(row['norm'])
            if pair(k0,g,r)>0:
                if any(pair(r,g,e)!=0 for e in edge):continue
                r=tuple(-z for z in r)
            r9=(r[0],*r[1:5],*(Z,)*4,r[5])
            group=finite_orbit(edge+[r],edge+[r],g)
            if set(group)!=set(full)or not compatible(r9,q,walls9,g9):continue
            options.append((pair(k0,g,r)**2/q,q,r))
        assert options,'first crossing has no complete compatible chamber sector'
        _,q,r=min(options);assert pair(point,g,point)<0 and all(pair(point,g,v)<=0 for v in walls)
        ng=tuple(tuple(pair(a,g,b)for b in edge+[r])for a in edge+[r])
        endpoint={'parameter':first.to_json(),'point':enc(point),'norm':pair(point,g,point).to_json(),'root':enc(r),'root_norm':q.to_json(),
            'root_height':(ALPHA*r[0]*r[0]/q).to_json(),'incident_gram':menc(ng),'full_stabilizer_mirrors':len(full),'H9_stabilizer_mirrors':len(full)+24}
    out={'status':'COMPLETE_FIRST_EDGE_MIRROR'if endpoint else'COMPLETE_PROBE_NO_NEW_MIRROR',
        'source':str(SOURCE.relative_to(ROOT)),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'edge_index':edge_index,'edge_indices':ids,'released_wall':released,'gram':menc(g),'edge_roots':menc(edge),
        'start':enc(k0),'probe':enc(k1),'target_cosh_squared':Cvalue,'actual_cosh_squared':C.to_json(),'probe_axis_parameter':step.to_json(),
        'released_wall_residual':enc(b),'residual_norm':B.to_json(),'norm_bounds':bounds,'profiles':profiles,
        'span_roots':menc(list(span.values())),'folded_segment_roots':records,'expanded_segment_mirrors':len(expanded),
        'edge_mirrors':len(edge_group),'complete_start_mirrors':len(start_group),'endpoint':endpoint,
        'H9_lift_source':'results/fresh_h9_new_host_20261013_first.json',
        'H9_lift_source_sha256':hashlib.sha256((ROOT/'results/fresh_h9_new_host_20261013_first.json').read_bytes()).hexdigest(),
        'stats':dict(stats),'seconds':round(time.monotonic()-started,3),
        'scope':'Complete all-norm shifted integral-pairing census on one finite segment leaving a previously certified full corner. Both radical signs and all residual classes are included. Full finite edge unfolding checks the start and, when present, first endpoint stabilizer. The unchanged complete second-F4 census lifts the interval statement to H9. A no-crossing probe leaves the rest of the edge open.'}
    output.write_text(json.dumps(out,indent=2)+'\n');print({'status':out['status'],'edge':edge_index,'C':Cvalue,'profiles':len(profiles),'stats':dict(stats),'endpoint':endpoint,'seconds':out['seconds']},flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--edge',type=int,choices=range(4),required=True);p.add_argument('--cosh2',type=int,default=4);args=p.parse_args();run(args.edge,args.cosh2)
