#!/usr/bin/env python3
"""Reuse the certified two-field orbit engine on seven new outer edges."""
import argparse,hashlib,json
import fresh_h9_a9_next_20261013_orbits as engine
from fresh_h9_a9_outer_20261013_profiles import ROOT,SOURCE,OUT as PROFILE,FIRST
OUT=ROOT/'results/fresh_h9_a9_outer_20261013_orbits.json'

def run(scout_powers=8):
    pd=json.loads(PROFILE.read_text());assert pd['status']=='COMPLETE_SEVEN_OUTER_EDGE_FINITE_PROFILES'and len(pd['cases'])==7
    assert pd['source_corner_complete_mirrors']==21
    for p,h in pd['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    # This engine uses only V²=-20alpha, b²=B and the saved covectors;
    # b may have nonzero time coordinate. Normalizing V-sb by its actual
    # time coordinate keeps its exact parameterization valid here.
    engine.PROFILE=PROFILE;engine.SOURCE=SOURCE;engine.OUT=OUT
    engine.run(scout_powers)
    out=json.loads(OUT.read_text());assert len(out['cases'])==7
    out['status']='COMPLETE_SEVEN_OUTER_EDGE_NORM_INTERFACE_ROUND'
    out['sources'][str(FIRST.relative_to(ROOT))]=hashlib.sha256(FIRST.read_bytes()).hexdigest()
    out['scope']='Seven actual outgoing edges at a previously completely certified outer A4+A4+A1 corner. Every finite profile and every relative-norm orbit through its selected upper endpoint is exhaustive, with the prior certified PARI interfaces reused exactly. The bounded unit scout only proposes an upper endpoint; unsuccessful cases remain bounded. This batch keeps the central H8 facet and does not examine the direction releasing it.'
    OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],[(c['released_wall'],c['status'],len(c.get('visited_points',[])))for c in out['cases']],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--scout-powers',type=int,default=8);run(p.parse_args().scout_powers)
