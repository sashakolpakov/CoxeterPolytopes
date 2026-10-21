#!/usr/bin/env python3
"""Independent Fraction root census of a finite quartic edge segment."""
from fresh_h9_quartic_edge_walk_20261013_independent import (
    ROOT, Q, Z, O, T, ALPHA, POOL, Host, add, sub, neg, mul, scale, sq,
    inv, divide, dot, mv, integral, sign, dec, vec, enc, matrix, pd, a, old, ALLOWED)
from fractions import Fraction
from math import isqrt
from itertools import product
from collections import Counter
import argparse
import hashlib
import json
import time
from functools import cmp_to_key

SOURCE='results/fresh_h9_quartic_first_endpoint_20261013.json'
OUT=ROOT/'results/fresh_h9_quartic_first_endpoint_20261013_independent.json'


def rational_root(x):
    if x<0:return None
    n,d=isqrt(x.numerator),isqrt(x.denominator)
    return Q(n,d)if n*n==x.numerator and d*d==x.denominator else None


def qm(x,y):return x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def qd(x,y):
    z=qm(x,(y[0],-y[1]));n=y[0]**2-2*y[1]**2
    return z[0]/n,z[1]/n


def quadratic_roots(x):
    """All square roots in Q(sqrt2), from the two coefficient equations."""
    u,v=x;out=set()
    if not v:
        r=rational_root(u)
        if r is not None:out.update(((r,Q(0)),(-r,Q(0))))
        r=rational_root(u/2)
        if r is not None:out.update(((Q(0),r),(Q(0),-r)))
    else:
        delta=rational_root(u*u-2*v*v)
        if delta is not None:
            for d in {delta,-delta}:
                p=rational_root((u+d)/2)
                if p:
                    out.update(((p,v/(2*p)),(-p,-v/(2*p))))
    assert all(qm(r,r)==x for r in out)
    return out


def quartic_roots(x):
    """All K-square roots; independent implementation, no producer square gate."""
    u=(x[0]+2*x[2],x[2]);v=(x[1]+2*x[3],x[3]);d=(Q(2),Q(1));solutions=[]
    if not any(v):
        solutions += [(p,(Q(0),Q(0)))for p in quadratic_roots(u)]
        solutions += [((Q(0),Q(0)),q)for q in quadratic_roots(qd(u,d))]
    else:
        aa=qm(u,u);bb=qm(d,qm(v,v));n=(aa[0]-bb[0],aa[1]-bb[1])
        for s in quadratic_roots(n):
            for p in quadratic_roots(((u[0]+s[0])/2,(u[1]+s[1])/2)):
                if any(p):solutions.append((p,qd(v,(2*p[0],2*p[1]))))
    out={(p[0]-2*p[1],q[0]-2*q[1],p[1],q[1])for p,q in solutions}
    assert all(sq(r)==x for r in out)
    return out


def line(r):
    first=next(z for z in r if z!=Z)
    return tuple(divide(z,first)for z in r)


def reflect(r,e,host):
    factor=divide(scale(host.pair(r,e),2),host.pair(e,e))
    return tuple(sub(x,mul(factor,y))for x,y in zip(r,e))


def closure(roots,generators,host):
    out={line(r):r for r in roots};queue=list(out.values())
    for r in queue:
        for e in generators:
            s=reflect(r,e,host);key=line(s)
            if key not in out:out[key]=s;queue.append(s)
        assert len(out)<1000
    return out


def primitive(r,q,host):
    return (all(map(integral,r)) and not all(integral(divide(z,T))for z in r)
            and all(integral(divide(scale(z,2),q))for z in mv(host.g,r)))


def normalized_corner(roots,host):
    cov=[mv(host.g,r)for r in roots]
    v=(O,*mv(a.minverse(tuple(tuple(c[j]for j in range(1,6))for c in cov)),
              tuple(neg(c[0])for c in cov)))
    assert all(host.pair(v,r)==Z for r in roots)and sign(host.pair(v,v))<0
    pd(tuple(tuple(host.pair(r,s)for s in roots)for r in roots))
    return v


def replay(stage=1):
    suffix={1:'',2:'_second',3:'_third'}[stage]
    source=SOURCE.replace('.json',suffix+'.json')
    start=time.monotonic();data=json.loads((ROOT/source).read_text());host=Host(1)
    expected={
        1:('results/fresh_h9_quartic_mixed_20261013_endpoint.json',None,
           'results/fresh_h9_quartic_edge_walk_20261013_profiles.json',[1,2,4,12],6),
        2:('results/fresh_h9_quartic_edge_walk_20261013_profiles.json',0,
           'results/fresh_h9_quartic_edge_walk_20261013_next.json',[1,2,4,70],12),
        3:('results/fresh_h9_quartic_edge_walk_20261013_next.json',0,
           'results/fresh_h9_quartic_edge_walk_20261013_f4.json',[1,2,70,72],4)}[stage]
    assert data['stage']==stage
    assert tuple(data[k]for k in('base_source','base_source_case','proposal_source','edge_indices','start_corner_wall'))==expected
    assert set(data['source_sha256'])=={expected[0],expected[2],'results/fresh_h9_glued_host_20261013_structure.json'}
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in data['source_sha256'].items())
    base=json.loads((ROOT/data['base_source']).read_text())
    if data['base_source_case']is not None:base=base['cases'][data['base_source_case']]
    walls=tuple(vec(r)for r in base['forced_H5_face']['all_walls'])
    edge=tuple(walls[i]for i in data['edge_indices']);assert edge==matrix(data['edge_roots'])
    eg=tuple(tuple(host.pair(r,s)for s in edge)for r in edge);pd(eg);ei=a.minverse(eg)
    proposal=json.loads((ROOT/data['proposal_source']).read_text())
    raw=proposal['compatible_dictionary'][proposal['cases'][0]['dictionary_index']]['vector']
    r9=vec(raw);proposed=(*r9[:5],r9[9])
    start_wall=data['start_corner_wall']
    k0=normalized_corner((*edge,walls[start_wall]),host);k1=normalized_corner((*edge,proposed),host)
    assert k0==vec(data['start_corner'])and k1==vec(data['proposed_endpoint'])
    assert all(sign(host.pair(k0,r))<=0 and sign(host.pair(k1,r))<=0 for r in walls)
    assert sign(host.pair(k1,walls[start_wall]))<0
    h0=sub(divide(ALPHA,neg(host.pair(k0,k0))),O)
    h1=sub(divide(ALPHA,neg(host.pair(k1,k1))),O);H=h0 if sign(sub(h0,h1))>=0 else h1
    assert h0==dec(data['start_height_cap'])and h1==dec(data['end_height_cap'])and H==dec(data['height_cap'])
    # Solve five linear equations (time plus four edge pairings) directly.
    equations=[(O,*(Z,)*5)]+[mv(host.g,e)for e in edge]
    # Any nonzero 5x5 minor leaves a complete free coordinate over K.
    chosen=None
    for free in range(5,-1,-1):
        pivots=[j for j in range(6)if j!=free]
        try:small_inverse=a.minverse(tuple(tuple(row[j]for j in pivots)for row in equations))
        except(StopIteration,ZeroDivisionError):continue
        chosen=free,pivots,small_inverse;break
    assert chosen is not None
    free,pivots,small_inverse=chosen
    direction=[Z]*6;direction[free]=O
    for j,z in zip(pivots,mv(small_inverse,tuple(neg(row[free])for row in equations))):direction[j]=z
    direction=tuple(direction);qa=host.pair(direction,direction)
    assert all(sign(qa,e)>0 for e in range(4))and direction[0]==Z
    profiles={};span={};found={};counts=Counter();per_norm=[]
    saved_profiles={(dec(row['norm']),vec(row['fixed_magnitudes'])):row for row in data['profiles']}
    for q,bounds in zip(POOL,data['norm_bounds'],strict=True):
        assert q==dec(bounds['norm'])
        alph=[a.alphabet((mul(q,host.pair(e,e)),)*4)for e in edge]
        assert [set(xs)for xs in alph]==[set(vec(xs))for xs in bounds['pairing_alphabets']]
        caps=(divide(mul(H,q),ALPHA),*(divide(neg(q),ALPHA)for _ in range(3)))
        assert caps==vec(bounds['time_caps'])
        positive=old.balanced(caps);times={x for z in positive for x in(z,neg(z))}
        assert times==set(vec(bounds['time_alphabet']))
        qp_count=0
        for fixed in product(*alph):
            coeff=mv(ei,tuple(neg(z)for z in fixed))
            r0=tuple(dot(coeff,tuple(e[j]for e in edge))for j in range(6))
            residual=sub(q,host.pair(r0,r0))
            if residual==Z:
                if primitive(r0,q,host):span[line(r0)]=r0
                continue
            if any(sign(residual,e)<=0 for e in range(4)):continue
            key=q,fixed;profiles[key]=(r0,residual);row=saved_profiles[key];qp_count+=1
            assert r0==vec(row['fixed_component'])and residual==dec(row['residual_norm'])
            assert bool(quartic_roots(residual))==row['residual_square']
            for x in times:
                solved=mv(small_inverse,(x,*map(neg,fixed)));anchor=[Z]*6
                for j,z in zip(pivots,solved):anchor[j]=z
                anchor=tuple(anchor);qb=scale(host.pair(anchor,direction),2);qc=sub(host.pair(anchor,anchor),q)
                D=sub(sq(qb),scale(mul(qa,qc),4));counts['quadratics']+=1
                roots=quartic_roots(D)
                if roots:counts['square_discriminants']+=1
                for radical in roots:
                    t=divide(sub(radical,qb),scale(qa,2))
                    r=tuple(add(z,mul(t,d))for z,d in zip(anchor,direction))
                    assert host.pair(r,r)==q and tuple(neg(host.pair(r,e))for e in edge)==fixed
                    if not primitive(r,q,host):continue
                    counts['primitive_integral_solutions']+=1
                    pa,pb=host.pair(k0,r),host.pair(k1,r)
                    if sign(mul(pa,pb))>0:continue
                    assert pa!=Z or pb!=Z
                    lam=divide(neg(pa),sub(pb,pa));assert sign(lam)>=0 and sign(sub(O,lam))>=0
                    found[line(r)]=(r,q,lam,pa,pb)
        per_norm.append({'norm':enc(q),'time_alphabet':len(times),'integral_profiles':qp_count})
    assert set(profiles)==set(saved_profiles)
    assert counts==data['stats']
    assert set(span)=={line(vec(r))for r in data['span_folded_roots']}
    expected={line(vec(r['vector'])):(dec(r['norm']),dec(r['segment_parameter']))for r in data['folded_segment_roots']}
    assert {key:(row[1],row[2])for key,row in found.items()}==expected
    edge_group=closure(edge,edge,host);assert set(closure(list(span.values()),edge,host))==set(edge_group)
    all_segment=closure([r[0]for r in found.values()],edge,host)
    initial={key for key,r in all_segment.items()if host.pair(k0,r)==Z}|set(edge_group)
    expected_initial=set(closure((*edge,walls[start_wall]),(*edge,walls[start_wall]),host))
    assert initial==expected_initial and data['complete_start_stabilizer']
    parameters={row[2]for row in found.values()}
    first=min((x for x in parameters if sign(x)>0),key=cmp_to_key(lambda x,y:sign(sub(x,y))))
    assert first==dec(data['nearest_parameter'])
    point=tuple(add(mul(sub(O,first),x),mul(first,y))for x,y in zip(k0,k1))
    near=vec(data['nearest_root']);q=dec(data['nearest_norm'])
    assert (line(near)==line(proposed))==data['nearest_is_proposed_mirror']
    assert host.pair(near,near)==q and host.pair(point,near)==Z and sign(host.pair(k0,near))<0
    closest=[row for row in found.values()if row[2]==first and
             (sign(row[3])<0 or all(host.pair(row[0],e)==Z for e in edge))]
    for r,n,lam,pa,pb in closest:
        chosen=divide(sq(host.pair(k0,near)),q);priority=divide(sq(pa),n)
        assert sign(sub(priority,chosen))>=0
        if chosen==priority:assert sign(sub(n,q))>=0
    nearest={key for key,r in all_segment.items()if host.pair(point,r)==Z}|set(edge_group)
    assert nearest==set(closure((*edge,near),(*edge,near),host))
    ng=tuple(tuple(host.pair(r,s)for s in(*edge,near))for r in(*edge,near));pd(ng)
    assert ng==matrix(data['nearest_incident_gram'])and point==vec(data['nearest_point'])
    assert host.pair(point,point)==dec(data['nearest_point_norm'])
    for roots in((*edge,walls[start_wall]),(*edge,near)):
        for i,r in enumerate(roots):
            for s in roots[:i]:
                p=host.pair(r,s);assert sign(p)<=0
                assert divide(scale(sq(p),4),mul(host.pair(r,r),host.pair(s,s)))in ALLOWED
    assert [len(edge_group),len(initial),len(nearest),len(all_segment)]==[
        data['edge_stabilizer_mirrors'],data['start_stabilizer_mirrors'],
        data['nearest_full_stabilizer_mirrors'],data['expanded_segment_mirror_count']]
    for nq,packet in zip(POOL,data['H9_lifting_block_census'],strict=True):
        assert nq==dec(packet['norm'])
        bs=a.block_vectors(nq)
        assert set(bs)=={(dec(row['norm']),vec(row['vector']))for row in packet['bounded_dominant_F4_vectors']}
        assert all(en in(Z,nq)for en,v in bs)
        assert all(v==(Z,)*4 for en,v in bs if all(sign(sub(nq,en),e)>0 for e in range(4)))
    h9=Host(2);lift=lambda r:(r[0],*r[1:5],*(Z,)*4,r[5])
    other=[tuple(neg(O)if i==j else Z for i in range(10))for j in range(5,9)]
    other_mirrors=closure(other,other,h9);assert len(other_mirrors)==24
    for generators,smallkeys in((edge,set(edge_group)),((*edge,walls[start_wall]),initial),((*edge,near),nearest)):
        both=[*map(lift,generators),*other]
        expected={line(lift(r))for r in smallkeys}|set(other_mirrors)
        assert set(closure(both,both,h9))==expected
    full_counts=[len(edge_group)+24,len(initial)+24,len(nearest)+24]
    assert full_counts==[data[k]for k in('H9_edge_stabilizer_mirrors','H9_start_stabilizer_mirrors','H9_end_stabilizer_mirrors')]
    return {'status':'PASS_INDEPENDENT_COMPLETE_FIRST_ENDPOINT_CENSUS','source':source,
            'source_sha256':hashlib.sha256((ROOT/source).read_bytes()).hexdigest(),
            'norm_packets':per_norm,'integral_profiles':len(profiles),
            'nonsquare_residual_profiles':sum(not row['residual_square']for row in data['profiles']),
            'quadratics':counts['quadratics'],'primitive_integral_solutions':counts['primitive_integral_solutions'],
            'folded_segment_mirrors':len(found),'expanded_segment_mirrors':len(all_segment),
            'mirrors_meeting_open_segment':sum(Z!=row[2]!=O for row in found.values()),
            'edge_mirrors':len(edge_group),'complete_start_stabilizer_mirrors':len(initial),
            'complete_end_stabilizer_mirrors':len(nearest),'nearest_is_proposed_mirror':data['nearest_is_proposed_mirror'],
            'H9_stabilizer_mirror_counts_edge_start_end':full_counts,
            'seconds':round(time.monotonic()-start,3),
            'scope':'Separate Fraction arithmetic, independently chosen affine free coordinate, full time/pairing alphabets, independent complete quadratic-tower square test, full primitive integrality and finite reflection-orbit comparisons. Independent dominant F4 census excludes mixed-block crossing roots and lifts the complete step and stabilizer checks to H9. No volume or full-chamber completion claim.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--stage',type=int,choices=(1,2,3),default=1)
    args=p.parse_args();result=replay(args.stage)
    suffix={1:'',2:'_second',3:'_third'}[args.stage]
    output=OUT.with_name(OUT.stem.replace('_independent',suffix+'_independent')+'.json')
    if args.verify:
        saved=json.loads(output.read_text());saved.pop('seconds');result.pop('seconds');assert saved==result
    else:output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
