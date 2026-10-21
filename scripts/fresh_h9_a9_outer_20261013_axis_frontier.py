#!/usr/bin/env python3
"""Exact existing-wall lookahead after the complete outer H9 vertex star."""
import argparse,hashlib,json
from collections import Counter
from fresh_h9_a9_outer_20261013_profiles import ROOT,K,Z,O,G,inverse,mv,dot,pair,enc,menc

SOURCE=ROOT/'results/fresh_h9_a9_outer_20261013_axis_geometry.json'
OLD=ROOT/'results/fresh_h9_a9_outer_20261013_frontier.json'
OUT=ROOT/'results/fresh_h9_a9_outer_20261013_axis_frontier.json'

def run(verify=False):
    data,old=[json.loads(p.read_text())for p in(SOURCE,OLD)]
    walls=[tuple(map(K.from_json,r))for r in data['all_walls']]
    cov=[mv(G,r)for r in walls]
    inside=tuple(map(K.from_json,data['strict_interior']))
    points=[tuple(map(K.from_json,r['point']))for r in old['vertices']]
    names=[r['name']for r in old['vertices']]
    for r in data['actual_finite_vertices']:
        p=tuple(map(K.from_json,r['point']));assert p not in points
        points.append(p);names.append(f"axis{r['case']}{'_reversed'if r['reversed']else''}")
    rev=lambda r:(r[0],*reversed(r[1:]))
    wall_swap=[walls.index(rev(r))for r in walls]
    point_swap=[points.index(rev(p))for p in points]
    records=[];all_rows=[]
    for vi,p in enumerate(points):
        assert p[0]==O and pair(p,G,p)<0 and pair(p,G,inside)<0
        values=[dot(c,p)for c in cov];assert all(v<=0 for v in values)
        ids=[i for i,v in enumerate(values)if v==Z];assert len(ids)==9
        gram=tuple(tuple(pair(walls[i],G,walls[j])for j in ids)for i in ids);gi=inverse(gram)
        rows=[]
        for k,released in enumerate(ids):
            u=tuple(-sum((gi[j][k]*walls[i][z]for j,i in enumerate(ids)),Z)for z in range(10))
            d=tuple(u[j]-u[0]*p[j]for j in range(10));assert d[0]==Z
            deriv=[dot(c,d)for c in cov]
            assert all(deriv[i]==(-O if i==released else Z)for i in ids)
            hits=[(-values[i]/s,i)for i,s in enumerate(deriv)if s>0]
            row={'vertex':vi,'released_wall':released,'edge_indices':[i for i in ids if i!=released],'direction':enc(d)}
            if not hits:row['status']='NO_EXISTING_WALL_AHEAD'
            else:
                t=min(s for s,i in hits);assert t>0
                endpoint=tuple(x+t*y for x,y in zip(p,d));norm=pair(endpoint,G,endpoint)
                assert all(dot(c,endpoint)<=0 for c in cov)
                target=points.index(endpoint)if endpoint in points else None
                status='KNOWN_CERTIFIED_VERTEX'if target is not None else'FINITE_EXISTING_ENDPOINT'if norm<0 else'IDEAL_EXISTING_ENDPOINT'if norm==Z else'LIGHT_CONE_BEFORE_EXISTING_WALL'
                row.update(status=status,parameter=t.to_json(),point=enc(endpoint),norm=norm.to_json(),existing_walls=[i for s,i in hits if s==t],known_target=target)
            rows.append(row);all_rows.append(row)
        records.append({'name':names[vi],'point':enc(p),'norm':pair(p,G,p).to_json(),'incident_indices':ids,'outgoing_edges':rows})
    lookup={(r['vertex'],r['released_wall']):r for r in all_rows}
    orbits=[];seen=set()
    for key,row in lookup.items():
        if row['status']=='KNOWN_CERTIFIED_VERTEX':continue
        other=(point_swap[key[0]],wall_swap[key[1]])
        assert lookup[other]['status']==row['status']
        orbit=tuple(sorted({key,other}))
        if orbit in seen:continue
        seen.add(orbit);orbits.append([list(x)for x in orbit])
    counts=dict(Counter(r['status']for r in all_rows))
    known={(min(r['vertex'],r['known_target']),max(r['vertex'],r['known_target']))for r in all_rows if r['status']=='KNOWN_CERTIFIED_VERTEX'}
    expected=set(map(tuple,old['known_undirected_edges']))|{(2+i%2,24+i)for i in range(2)}
    assert expected<=known
    out={'status':'EXACT_TWENTY_SIX_CORNER_A9_FRONTIER','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,OLD)},
      'vertices':records,'wall_reversal':wall_swap,'vertex_reversal':point_swap,'known_undirected_edges':[list(x)for x in sorted(known)],
      'certified_path_edges':[list(x)for x in sorted(expected)],'outgoing_count':len(all_rows),'status_counts':counts,'unresolved_reversal_orbits':orbits,'unresolved_orbit_count':len(orbits),
      'scope':'Twenty-five source-anchored edges connect the twenty-six complete chamber corners. Any extra known-corner edge is also a full-chamber segment by convexity. Other rows are exact lookahead against only the thirty-five saved walls and do not exclude additional lattice mirrors before the light cone or a listed finite endpoint. Reversal is integral.'}
    if verify:assert json.loads(OUT.read_text())==out
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],counts,'unresolved orbits',len(orbits),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
