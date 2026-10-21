#!/usr/bin/env python3
"""Producer-free exact replay of the new A9 Pell cap and its edge frontier."""
import argparse,hashlib,json
from fresh_h9_new_host_next_20261013_a9_independent import a,ALPHA,G,POOL,P,compatible
SOURCE=a.ROOT/'results/fresh_h9_new_host_next_20261013_a9_endpoint.json'
OUT=a.ROOT/'results/fresh_h9_new_host_next_20261013_a9_endpoint_independent.json'

def run(verify=False):
    d=json.loads(SOURCE.read_text());base=a.ROOT/d['source'];assert hashlib.sha256(base.read_bytes()).hexdigest()==d['source_sha256']
    assert tuple(a.vec(r)for r in d['gram'])==G
    z,x=a.dec(d['z']),a.dec(d['time']);assert a.sub(a.sq(z),a.mul(a.scale(ALPHA,10),a.sq(x)))==a.scale(a.O,20)
    root=a.vec(d['root']);assert root==(x,*(a.scale(z,a.Q(i,10))for i in(1,2,3,4,5,4,3,2,1)))
    walls=[a.vec(r)for r in d['all_walls']];inside=a.vec(d['strict_interior']);assert walls[-1]==root
    for i,r in enumerate(walls):
        q=P(r,r);assert q in POOL and all(map(a.integral,r))
        assert not all(a.integral(a.divide(v,a.T))for v in r)and not all(a.integral(a.scale(v,a.Q(1,5)))for v in r)
        assert all(a.integral(a.divide(a.scale(v,2),q))for v in a.mv(G,r))and compatible(r,walls[:i])
    assert P(root,root)==a.scale(a.O,2)and a.sign(P(inside,inside))<0 and all(a.sign(P(inside,r))<0 for r in walls)
    point=a.vec(d['endpoint']);c=a.divide(a.mul(ALPHA,x),z)
    assert point==(a.O,*(a.scale(c,i)for i in(1,2,3,4,5,4,3,2,1)))
    assert P(point,point)==a.divide(a.scale(ALPHA,-20),a.sq(z))==a.dec(d['endpoint_norm'])and a.sign(P(point,point))<0
    assert all(a.sign(P(point,r))<=0 for r in walls)
    proof=d['actual_finite_vertex'];ids=proof['indices'];incident=[walls[i]for i in ids]
    gram=tuple(tuple(P(r,s)for s in incident)for r in incident);a.positive_definite(gram)
    assert gram==tuple(a.vec(row)for row in proof['gram'])and ids==[0,1,2,3,5,6,7,8,10]
    v=a.vec(proof['vector']);assert tuple(a.divide(t,v[0])for t in v)==point
    assert P(v,v)==a.dec(proof['norm'])and P(v,inside)==a.dec(proof['center_pairing'])
    assert tuple(P(v,r)for r in walls)==a.vec(proof['all_wall_pairings'])
    assert [i for i,r in enumerate(walls)if P(v,r)==a.Z]==proof['all_active_indices']==ids
    for row in d['outgoing_existing_wall_frontier']:
        released=row['released_wall'];edge=row['edge_indices'];direction=a.vec(row['direction'])
        assert set(edge)==set(ids)-{released}and direction[0]==a.Z
        assert all(P(direction,walls[i])==a.Z for i in edge)and P(direction,walls[released])==a.neg(a.O)
        hits=[(a.divide(a.neg(P(point,r)),P(direction,r)),i)for i,r in enumerate(walls)if a.sign(P(direction,r))>0]
        if row['status']=='NO_EXISTING_WALL_AHEAD':assert not hits;continue
        t=a.dec(row['parameter']);assert all(a.sign(a.sub(s,t))>=0 for s,i in hits)
        assert [i for s,i in hits if s==t]==row['existing_walls']
        endpoint=tuple(a.add(u,a.mul(t,v))for u,v in zip(point,direction))
        assert endpoint==a.vec(row['point'])and P(endpoint,endpoint)==a.dec(row['norm'])
        expected='FINITE_EXISTING_ENDPOINT'if a.sign(P(endpoint,endpoint))<0 else'IDEAL_EXISTING_ENDPOINT'if P(endpoint,endpoint)==a.Z else'LIGHT_CONE_BEFORE_EXISTING_WALL'
        assert row['status']==expected
    geom=d['geometry'];v=a.vec(geom['strict_vector']);assert P(v,v)==a.dec(geom['strict_norm'])and a.sign(P(v,v))>0
    assert a.sign(v[0])>0 and all(a.sign(P(v,r))<0 for r in walls)and a.sign(P(v,inside))<0
    result={'status':'PASS_INDEPENDENT_A9_ENDPOINT','source':str(SOURCE.relative_to(a.ROOT)),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
      'walls':len(walls),'all_pair_checks':len(walls)*(len(walls)-1)//2,'norm_identity':True,'actual_finite_vertex_rank':9,
      'strict_escape_norm':a.enc(P(v,v)),'outgoing_edges':len(d['outgoing_existing_wall_frontier']),
      'open_existing_wall_directions':sum(r['status']!='FINITE_EXISTING_ENDPOINT'for r in d['outgoing_existing_wall_frontier']),
      'scope':'Separate Fraction arithmetic replays the literal norm identity, every root/pair/interior/escape, actual finite vertex, and all existing-wall edge lookahead. No PARI call, nearest-mirror, or full-stabilizer claim.'}
    if verify:assert json.loads(OUT.read_text())==result
    else:OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],result['walls'],result['open_existing_wall_directions'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
