#!/usr/bin/env python3
"""Fraction-only replay of thirty-three facets, fourteen new vertices, and the edge queue."""
import argparse,hashlib,json
from collections import Counter
from fresh_h9_new_host_next_20261013_a9_independent import a,G,P,compatible

SOURCE=a.ROOT/'results/fresh_h9_a9_outer_20261013_geometry.json'
FRONTIER=a.ROOT/'results/fresh_h9_a9_outer_20261013_frontier.json'
OUT=a.ROOT/'results/fresh_h9_a9_outer_20261013_geometry_independent.json'

def run(verify=False):
    data,fd=[json.loads(p.read_text())for p in(SOURCE,FRONTIER)]
    for d in(data,fd):
        for p,h in d['sources'].items():assert hashlib.sha256((a.ROOT/p).read_bytes()).hexdigest()==h
    assert tuple(a.vec(r)for r in data['gram'])==G
    walls=[a.vec(r)for r in data['all_walls']];inside=a.vec(data['strict_interior'])
    assert len(walls)==33 and a.sign(P(inside,inside))<0 and a.sign(inside[0])>0
    for i,r in enumerate(walls):
        assert P(r,r)==a.scale(a.O,2)and all(map(a.integral,r))
        assert not all(a.integral(a.divide(x,a.T))for x in r)and not all(a.integral(a.scale(x,a.Q(1,5)))for x in r)
        assert all(a.integral(x)for x in a.mv(G,r))and compatible(r,walls[:i])and a.sign(P(inside,r))<0
    for row in data['actual_finite_vertices']:
        p=a.vec(row['point']);c=row['certificate'];ids=c['indices'];v=a.vec(c['vector'])
        assert p[0]==a.O and tuple(a.divide(x,v[0])for x in v)==p
        gram=tuple(tuple(P(walls[i],walls[j])for j in ids)for i in ids);a.positive_definite(gram)
        assert len(ids)==9 and gram==tuple(a.vec(r)for r in c['gram'])
        assert P(v,v)==a.dec(c['norm'])and a.sign(P(v,v))<0 and P(v,inside)==a.dec(c['center_pairing'])and a.sign(P(v,inside))<0
        assert tuple(P(v,r)for r in walls)==a.vec(c['all_wall_pairings'])and all(a.sign(P(v,r))<=0 for r in walls)
        assert [i for i,r in enumerate(walls)if P(v,r)==a.Z]==c['all_active_indices']==sorted(ids)
    geo=data['geometry'];escape=a.vec(geo['strict_vector'])
    assert P(escape,escape)==a.dec(geo['strict_norm'])and a.sign(P(escape,escape))>0
    assert a.sign(escape[0])>0 and P(escape,inside)==a.dec(geo['strict_center_pairing'])and a.sign(P(escape,inside))<0
    assert all(a.sign(P(escape,r))<0 for r in walls)
    points=[a.vec(r['point'])for r in fd['vertices']];rev=lambda r:(r[0],*reversed(r[1:]))
    assert fd['wall_reversal']==[walls.index(rev(r))for r in walls]
    assert fd['vertex_reversal']==[points.index(rev(p))for p in points]
    rows=[]
    for vi,row in enumerate(fd['vertices']):
        p=points[vi];ids=row['incident_indices'];assert p[0]==a.O and P(p,p)==a.dec(row['norm'])and a.sign(P(p,p))<0
        assert a.sign(P(p,inside))<0 and all(a.sign(P(p,r))<=0 for r in walls)
        assert ids==[i for i,r in enumerate(walls)if P(p,r)==a.Z]and len(ids)==9
        a.positive_definite(tuple(tuple(P(walls[i],walls[j])for j in ids)for i in ids))
        assert sorted(r['released_wall']for r in row['outgoing_edges'])==ids
        for r in row['outgoing_edges']:
            rows.append(r);d=a.vec(r['direction']);released=r['released_wall']
            assert r['vertex']==vi and set(r['edge_indices'])==set(ids)-{released}and d[0]==a.Z
            assert all(P(d,walls[i])==(a.neg(a.O)if i==released else a.Z)for i in ids)
            hits=[(a.divide(a.neg(P(p,w)),P(d,w)),i)for i,w in enumerate(walls)if a.sign(P(d,w))>0]
            if r['status']=='NO_EXISTING_WALL_AHEAD':assert not hits;continue
            t=a.dec(r['parameter']);assert a.sign(t)>0 and all(a.sign(a.sub(s,t))>=0 for s,i in hits)
            assert r['existing_walls']==[i for s,i in hits if s==t]
            end=tuple(a.add(x,a.mul(t,y))for x,y in zip(p,d));n=P(end,end)
            assert end==a.vec(r['point'])and n==a.dec(r['norm'])and all(a.sign(P(end,w))<=0 for w in walls)
            target=points.index(end)if end in points else None;assert target==r['known_target']
            expected='KNOWN_CERTIFIED_VERTEX'if target is not None else'FINITE_EXISTING_ENDPOINT'if a.sign(n)<0 else'IDEAL_EXISTING_ENDPOINT'if n==a.Z else'LIGHT_CONE_BEFORE_EXISTING_WALL'
            assert r['status']==expected
    assert fd['outgoing_count']==len(rows)==216 and fd['status_counts']==dict(Counter(r['status']for r in rows))
    lookup={(r['vertex'],r['released_wall']):r for r in rows};orbits=set()
    for key,r in lookup.items():
        if r['status']=='KNOWN_CERTIFIED_VERTEX':continue
        other=(fd['vertex_reversal'][key[0]],fd['wall_reversal'][key[1]])
        assert lookup[other]['status']==r['status'];orbits.add(tuple(sorted({key,other})))
    assert orbits=={tuple(map(tuple,r))for r in fd['unresolved_reversal_orbits']}and len(orbits)==fd['unresolved_orbit_count']
    known={(min(r['vertex'],r['known_target']),max(r['vertex'],r['known_target']))for r in rows if r['status']=='KNOWN_CERTIFIED_VERTEX'}
    expected={(0,1)}|{(1,i)for i in range(2,10)}|{(2+i%2,10+i)for i in range(14)}
    assert known==expected==set(map(tuple,fd['known_undirected_edges']))==set(map(tuple,fd['certified_path_edges']))
    # Each new norm-two facet is orthogonal to the central cap.
    assert all(P(walls[10],r)==a.Z for r in walls[11:])
    out={'status':'PASS_INDEPENDENT_A9_THIRTY_THREE_WALL_GEOMETRY_AND_FRONTIER','sources':{str(p.relative_to(a.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,FRONTIER)},
      'walls':33,'pair_checks':528,'new_actual_finite_vertices':14,'known_actual_vertices':24,'known_undirected_edges':23,'outgoing_directions':216,
      'frontier_status_counts':fd['status_counts'],'unresolved_reversal_orbits':len(orbits),'strict_escape_norm':a.enc(P(escape,escape)),
      'scope':'Separate Fraction field operations replay all literal norms, integrality, finite-angle gates, interiors, fourteen rank-nine vertex certificates, strict escape, and every existing-wall frontier direction. The nearest-mirror completeness comes from the separate seven-profile/orbit audit, not from this finite-wall lookahead.'}
    if verify:assert json.loads(OUT.read_text())==out
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['frontier_status_counts'],len(orbits),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
