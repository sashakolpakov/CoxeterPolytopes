#!/usr/bin/env python3
"""Separate Fraction enumeration of the new host's complete first edge."""
import argparse,hashlib,json
from itertools import product
from collections import Counter
import fresh_h9_new_host_20261013_independent as a
from fresh_h9_new_host_20261013_edge_independent import balanced
from fresh_h9_quartic_first_endpoint_20261013_independent import quartic_roots,closure,line,primitive

SOURCE=a.ROOT/'results/fresh_h9_new_host_20261013_first.json'
OUT=a.ROOT/'results/fresh_h9_new_host_20261013_first_independent.json'

def replay(verify=False):
    data=json.loads(SOURCE.read_text());host=a.Host(1)
    for p,h in data['sources'].items():assert hashlib.sha256((a.ROOT/p).read_bytes()).hexdigest()==h
    assert tuple(a.vec(r)for r in data['gram'])==host.g and a.dec(data['alpha'])==a.ALPHA
    fd=json.loads((a.ROOT/'results/fresh_h9_new_host_20261013_h5_h100.json').read_text())
    walls=[a.vec(r)for r in fd['all_walls']];edge=[walls[i]for i in data['edge_indices']]
    assert edge==[a.vec(r)for r in data['edge_roots']]
    k0,k1=a.vec(data['start']),a.vec(data['end']);assert k0==(a.O,*(a.Z,)*5)
    assert k1[0]==a.O and a.sign(host.pair(k1,k1))<0 and host.pair(k1,k1)==a.dec(data['end_norm'])
    assert all(host.pair(k0,e)==a.Z and host.pair(k1,e)==a.Z for e in edge)
    assert all(a.sign(host.pair(k1,w))<=0 for w in walls)
    H=a.sub(a.divide(a.ALPHA,a.neg(host.pair(k1,k1))),a.O);assert H==a.dec(data['height_cap'])
    eg=tuple(tuple(host.pair(r,s)for s in edge)for r in edge);a.positive_definite(eg);ei=a.minverse(eg)
    # Direct five-equation elimination, rather than the producer's fixed
    # fundamental-weight completed square.
    eq=[(a.O,*(a.Z,)*5)]+[a.mv(host.g,e)for e in edge]
    for free in range(5,-1,-1):
        pivots=[j for j in range(6)if j!=free]
        try:mi=a.minverse(tuple(tuple(row[j]for j in pivots)for row in eq))
        except(StopIteration,ZeroDivisionError):continue
        break
    direction=[a.Z]*6;direction[free]=a.O
    for j,z in zip(pivots,a.mv(mi,tuple(a.neg(row[free])for row in eq))):direction[j]=z
    direction=tuple(direction);qa=host.pair(direction,direction);assert all(a.sign(qa,e)>0 for e in range(4))
    saved={(a.dec(p['norm']),a.vec(p['fixed_magnitudes'])):p for p in data['profiles']}
    profiles=set();span={};found={};counts=Counter()
    for q,bound in zip(a.POOL,data['norm_bounds'],strict=True):
        assert q==a.dec(bound['norm'])
        alph=[a.alphabet((a.mul(q,host.pair(e,e)),)*4)for e in edge]
        assert [set(v)for v in alph]==[set(a.vec(v))for v in bound['pairing_alphabets']]
        for fixed in product(*alph):
            co=a.mv(ei,tuple(a.neg(v)for v in fixed));r0=tuple(a.dot(co,tuple(e[j]for e in edge))for j in range(6));R=a.sub(q,host.pair(r0,r0))
            if R==a.Z:
                if primitive(r0,q,host):span[line(r0)]=r0
                continue
            if any(a.sign(R,e)<=0 for e in range(4)):continue
            key=(q,fixed);profiles.add(key);row=saved[key]
            assert a.vec(row['fixed_component'])==r0 and a.dec(row['residual'])==R
            caps=(a.divide(a.mul(H,q),a.ALPHA),*(a.divide(a.neg(R),a.ALPHA)for _ in range(3)))
            times={z for x in balanced(caps)for z in(x,a.neg(x))}
            assert caps==a.vec(row['time_caps'])and times==set(a.vec(row['time_coordinates']))
            for x in times:
                solved=a.mv(mi,(x,*map(a.neg,fixed)));anchor=[a.Z]*6
                for j,z in zip(pivots,solved):anchor[j]=z
                anchor=tuple(anchor);qb=a.scale(host.pair(anchor,direction),2);qc=a.sub(host.pair(anchor,anchor),q)
                discr=a.sub(a.sq(qb),a.scale(a.mul(qa,qc),4));counts['quadratics']+=1
                roots=quartic_roots(discr)
                if roots:counts['field_squares']+=1
                for radical in roots:
                    t=a.divide(a.sub(radical,qb),a.scale(qa,2));r=tuple(a.add(v,a.mul(t,d))for v,d in zip(anchor,direction))
                    assert host.pair(r,r)==q and tuple(a.neg(host.pair(r,e))for e in edge)==fixed
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
    egroots=closure(edge,edge,host);assert set(closure(list(span.values()),edge,host))==set(egroots)
    expanded=closure([r[0]for r in found.values()],edge,host)
    initial={key for key,r in expanded.items()if host.pair(k0,r)==a.Z}|set(egroots)
    end={key for key,r in expanded.items()if host.pair(k1,r)==a.Z}|set(egroots)
    assert initial==set(closure(walls[:5],walls[:5],host))
    assert {r[2]for r in found.values()if a.sign(r[2])>0}=={a.O}
    assert a.dec(data['nearest_parameter'])==a.O
    selected=a.vec(data['selected_choice']['vector']);q=a.dec(data['selected_choice']['norm'])
    assert primitive(selected,q,host)and host.pair(selected,selected)==q and host.pair(k1,selected)==a.Z
    assert host.compatible(selected,walls)
    assert end==set(closure(edge+[selected],edge+[selected],host))
    gram=tuple(tuple(host.pair(r,s)for s in edge+[selected])for r in edge+[selected]);a.positive_definite(gram)
    assert gram==tuple(a.vec(r)for r in data['selected_incident_gram'])
    for i,r in enumerate(edge+[selected]):assert host.compatible(r,(edge+[selected])[:i])
    assert [len(egroots),len(initial),len(end),len(expanded)]==[data[k]for k in('edge_mirrors','complete_initial_mirrors','complete_endpoint_mirrors','expanded_segment_mirrors')]
    for q,packet in zip(a.POOL,data['H9_other_F4_census'],strict=True):
        assert q==a.dec(packet['norm']);bs=a.block_vectors(q)
        assert set(bs)=={(a.dec(r['norm']),a.vec(r['vector']))for r in packet['dominant_F4_vectors']}
        assert all(en in(a.Z,q)for en,v in bs)
    h9=a.Host(2);lift=lambda r:(r[0],*r[1:5],*(a.Z,)*4,r[5])
    other=[tuple(a.neg(a.O)if i==j else a.Z for i in range(10))for j in range(5,9)]
    othergroup=closure(other,other,h9);assert len(othergroup)==24
    for gen,small in((edge,egroots),(walls[:5],{k:k for k in initial}),(edge+[selected],{k:k for k in end})):
        allgen=[*map(lift,gen),*other];expected={line(lift(r))for r in small}|set(othergroup)
        assert set(closure(allgen,allgen,h9))==expected
    assert [len(egroots)+24,len(initial)+24,len(end)+24]==[data[k]for k in('H9_edge_mirrors','H9_initial_mirrors','H9_endpoint_mirrors')]
    result={'status':'PASS_INDEPENDENT_FULL_H9_FIRST_EDGE','source':str(SOURCE.relative_to(a.ROOT)),
        'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'profiles':len(profiles),'stats':dict(counts),
        'folded_segment_mirrors':len(found),'expanded_segment_mirrors':len(expanded),
        'H5_stabilizer_mirrors':[len(egroots),len(initial),len(end)],'H9_stabilizer_mirrors':[len(egroots)+24,len(initial)+24,len(end)+24],
        'nearest_root':[a.enc(v)for v in selected],'nearest_norm':a.enc(q),
        'arithmetic':'Separate Fraction field arithmetic; independent greedy unit-balanced trace alphabets; direct five-equation elimination and complete quadratic-tower roots; full finite orbit/stabilizer replay. No PARI norm-orbit completeness is needed for this first-edge certificate.'}
    if verify:assert json.loads(OUT.read_text())==result
    else:OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],result['stats'],result['H9_stabilizer_mirrors'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');replay(p.parse_args().verify)
