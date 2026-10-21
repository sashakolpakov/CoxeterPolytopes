#!/usr/bin/env python3
"""Independent all-norm shifted-pairing edge census in the changed host."""
import argparse,hashlib,json,time
from itertools import product
from collections import Counter
from functools import cmp_to_key
from pathlib import Path
import fresh_h9_new_host_20261013_independent as a
from fresh_h9_new_host_20261013_edge_independent import balanced
from fresh_h9_quartic_first_endpoint_20261013_independent import quartic_roots,closure,line,primitive


def replay(path):
    data=json.loads(path.read_text());host=a.Host(1);h9=a.Host(2)
    source=a.ROOT/data['source'];assert hashlib.sha256(source.read_bytes()).hexdigest()==data['source_sha256']
    base=json.loads(source.read_text());face=base['forced_H5_face'];walls=[a.vec(r)for r in face['all_walls']];walls9=[a.vec(r)for r in base['all_walls']]
    edge=[walls[i]for i in data['edge_indices']];w=walls[data['released_wall']]
    assert edge==[a.vec(r)for r in data['edge_roots']]and tuple(a.vec(r)for r in data['gram'])==host.g
    k0,k1=a.vec(data['start']),a.vec(data['probe']);N=host.pair(k0,k0)
    assert a.sign(N)<0 and a.sign(host.pair(k1,k1))<0 and a.sign(host.pair(k0,k1))<0
    assert all(host.pair(k,e)==a.Z for k in(k0,k1)for e in edge)and host.pair(k0,w)==a.Z and a.sign(host.pair(k1,w))<0
    assert all(a.sign(host.pair(k,r))<=0 for k in(k0,k1)for r in walls)
    eg=tuple(tuple(host.pair(r,s)for s in edge)for r in edge);a.positive_definite(eg);ei=a.minverse(eg)
    cf=a.mv(ei,tuple(host.pair(w,e)for e in edge));b=tuple(a.sub(w[j],a.dot(cf,tuple(e[j]for e in edge)))for j in range(6));B=host.pair(b,b)
    assert b==a.vec(data['released_wall_residual'])and B==a.dec(data['residual_norm'])and all(a.sign(B,e)>0 for e in range(4))
    C=a.divide(a.sq(host.pair(k0,k1)),a.mul(N,host.pair(k1,k1)))
    assert C==a.dec(data['actual_cosh_squared'])and a.sign(a.sub(C,a.O))>0 and a.sign(a.sub(a.scale(a.O,data['target_cosh_squared']),C))>=0
    step=a.dec(data['probe_axis_parameter']);raw=tuple(a.sub(x,a.mul(step,y))for x,y in zip(k0,b));assert tuple(a.divide(x,raw[0])for x in raw)==k1
    full=edge+[w];fullgram=tuple(tuple(host.pair(r,s)for s in full)for r in full);a.positive_definite(fullgram);full_inverse=a.minverse(fullgram)
    saved={(a.dec(row['norm']),a.vec(row['fixed_magnitudes'])):row for row in data['profiles']};profiles=set();span={};found={};counts=Counter()
    for q,bound in zip(a.POOL,data['norm_bounds'],strict=True):
        assert a.dec(bound['norm'])==q
        alph=[a.alphabet((a.mul(q,host.pair(e,e)),)*4)for e in edge]
        assert [set(v)for v in alph]==[set(a.vec(v))for v in bound['pairing_alphabets']]
        for fixed in product(*alph):
            cf=a.mv(ei,tuple(a.neg(v)for v in fixed));r0=tuple(a.dot(cf,tuple(e[j]for e in edge))for j in range(6));R=a.sub(q,host.pair(r0,r0))
            if R==a.Z:
                if primitive(r0,q,host):span[line(r0)]=r0
                continue
            if any(a.sign(R,e)<=0 for e in range(4)):continue
            key=q,fixed;profiles.add(key);row=saved[key]
            assert r0==a.vec(row['fixed_component'])and R==a.dec(row['residual'])
            offset=host.pair(r0,w);caps=(a.mul(a.mul(R,B),C),*(a.mul(R,B)for _ in range(3)))
            center=a.dec(row['integral_center']);assert a.integral(center)
            shift=a.sub(center,offset);loose=tuple(a.scale(a.add(cap,a.sq(shift)),2)for cap in caps)
            assert offset==a.dec(row['offset'])and caps==a.vec(row['caps'])and loose==a.vec(row['loose_caps'])
            values={a.add(center,z)for v in balanced(loose)for z in(v,a.neg(v))if all(a.sign(a.sub(caps[e],a.sq(a.sub(a.add(center,z),offset))),e)>=0 for e in range(4))}
            assert values==set(a.vec(row['pairing_values']))
            for value in values:
                # Independent anchor: invert the entire rank-five positive Gram.
                coeff=a.mv(full_inverse,(*map(a.neg,fixed),value))
                anchor=tuple(a.dot(coeff,tuple(r[j]for r in full))for j in range(6))
                assert host.pair(anchor,k0)==a.Z and host.pair(anchor,w)==value
                square=a.divide(a.sub(q,host.pair(anchor,anchor)),N);counts['quadratics']+=1;roots=quartic_roots(square)
                if roots:counts['field_squares']+=1
                for t in roots:
                    r=tuple(a.add(x,a.mul(t,y))for x,y in zip(anchor,k0));assert host.pair(r,r)==q
                    if not primitive(r,q,host):continue
                    counts['primitive_integral']+=1;pa,pb=host.pair(k0,r),host.pair(k1,r)
                    if a.sign(a.mul(pa,pb))>0:continue
                    assert pa!=a.Z or pb!=a.Z
                    lam=a.divide(a.neg(pa),a.sub(pb,pa));assert a.sign(lam)>=0 and a.sign(a.sub(a.O,lam))>=0
                    found[line(r)]=(r,q,lam)
    assert profiles==set(saved)and dict(counts)==data['stats']
    assert set(span)=={line(a.vec(r))for r in data['span_roots']}
    expected={line(a.vec(r['vector'])):(a.dec(r['norm']),a.dec(r['parameter']))for r in data['folded_segment_roots']}
    assert {key:(r[1],r[2])for key,r in found.items()}==expected
    edgegroup=closure(edge,edge,host);assert set(closure(list(span.values()),edge,host))==set(edgegroup)
    expanded=closure([r[0]for r in found.values()],edge,host)
    start={key for key,r in expanded.items()if host.pair(k0,r)==a.Z}|set(edgegroup)
    assert start==set(closure(full,full,host))
    assert [len(edgegroup),len(expanded),len(start)]==[data[k]for k in('edge_mirrors','expanded_segment_mirrors','complete_start_mirrors')]
    positive=[row[2]for row in found.values()if a.sign(row[2])>0];result={'source':str(path.relative_to(a.ROOT)),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'status':'PASS_'+data['status'],'profiles':len(profiles),'stats':dict(counts),'start_mirrors':len(start)}
    if not positive:assert data['endpoint']is None
    else:
        endpoint=data['endpoint'];first=min(positive,key=cmp_to_key(lambda x,y:a.sign(a.sub(x,y))));assert first==a.dec(endpoint['parameter'])
        point=tuple(a.add(x,a.mul(first,a.sub(y,x)))for x,y in zip(k0,k1));assert point==a.vec(endpoint['point'])and host.pair(point,point)==a.dec(endpoint['norm'])and a.sign(host.pair(point,point))<0
        r=a.vec(endpoint['root']);q=a.dec(endpoint['root_norm']);assert primitive(r,q,host)and host.pair(r,r)==q and host.pair(point,r)==a.Z
        assert a.sign(host.pair(k0,r))<0 and host.compatible(r,walls)
        r9=(r[0],*r[1:5],*(a.Z,)*4,r[5]);assert h9.compatible(r9,walls9)
        end={key for key,v in expanded.items()if host.pair(point,v)==a.Z}|set(edgegroup)
        assert end==set(closure(edge+[r],edge+[r],host))and len(end)==endpoint['full_stabilizer_mirrors']
        gram=tuple(tuple(host.pair(u,v)for v in edge+[r])for u in edge+[r]);a.positive_definite(gram)
        assert gram==tuple(a.vec(row)for row in endpoint['incident_gram'])
        result.update(endpoint_mirrors=len(end),H9_endpoint_mirrors=len(end)+24,root=[a.enc(v)for v in r])
    liftpath=a.ROOT/data['H9_lift_source'];assert hashlib.sha256(liftpath.read_bytes()).hexdigest()==data['H9_lift_source_sha256']
    liftdata=json.loads(liftpath.read_text())
    for q,packet in zip(a.POOL,liftdata['H9_other_F4_census'],strict=True):
        assert q==a.dec(packet['norm']);bs=a.block_vectors(q);assert set(bs)=={(a.dec(row['norm']),a.vec(row['vector']))for row in packet['dominant_F4_vectors']}
        assert all(en in(a.Z,q)for en,v in bs)
    result['complete_H9_lift']=True
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('sources',type=Path,nargs='+');p.add_argument('--output',type=Path,required=True);p.add_argument('--verify',action='store_true');args=p.parse_args()
    result={'status':'PASS_INDEPENDENT_SHIFTED_PAIRING_CENSUS','arithmetic':'Separate Fraction arithmetic and complete quadratic-tower roots; independent greedy unit-balanced trace enumeration; full rank-five Gram inversion for anchors; finite reflection unfolding.', 'cases':[replay(path.resolve())for path in args.sources]}
    if args.verify:assert json.loads(args.output.read_text())==result
    else:args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],len(result['cases']),flush=True)
