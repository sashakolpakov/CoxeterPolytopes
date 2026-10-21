#!/usr/bin/env python3
"""An exact cap beyond the cubic-seven height cutoff; it does not close H9."""
import json
from pathlib import Path

from fresh_cubic7_dihedral_h9_20261013 import K,T,Z,decode_matrix,enc,pair,mv,finite_angle_counts

ROOT=Path(__file__).resolve().parents[1]
source=ROOT/'results/fresh_cubic7_dihedral_h9_2theta_h30_20261013.json'
data=json.loads(source.read_text())
gram=decode_matrix(data['gram']);walls=decode_matrix(data['all_walls'])
alpha=2*T
# The sextic unit z-w*sqrt(alpha) has unit relative norm n.
# Its square divided by -n gives a relative norm-one unit.
z=(T+1)**2;w=T+2;n=z*z-alpha*w*w
assert n==2-T-T*T and abs(n.norm())==1
X=-2*z*w/n;Y=-(z*z+alpha*w*w)/n
assert X==6+18*T+8*T*T and Y==11+28*T+12*T*T
assert X.is_integral and Y.is_integral and X.sign()>0 and Y.sign()>0
assert Y*Y-alpha*X*X==1
root=(X,*([Z]*8),Y)
assert pair(root,gram,root)==1
assert all((2*a).is_integral for a in mv(gram,root))
wall_pairings=tuple(pair(root,gram,r) for r in walls)
assert all(p.sign()<=0 for p in wall_pairings)
inside=tuple(map(K.from_json,data['strict_interior']))
escape=tuple(map(K.from_json,data['geometry']['hyperideal_witness']['vector']))
assert pair(root,gram,inside).sign()<0 and pair(root,gram,escape).sign()<0
assert pair(escape,gram,escape).sign()>0
out={
    'status':'EXACT_COMPATIBLE_PELL_CAP_WITH_REMAINING_ESCAPE',
    'source':str(source.relative_to(ROOT)),
    'root':enc(root),'norm':K(1).to_json(),'height':(alpha*X*X).to_json(),
    'wall_pairings':enc(wall_pairings),
    'strict_interior_pairing':pair(root,gram,inside).to_json(),
    'retained_escape_pairing':pair(root,gram,escape).to_json(),
    'finite_angle_counts':finite_angle_counts((*walls,root),gram),
    'scope':'An exact wall beyond height30; neither a completed Vinberg enumeration nor a nonreflectivity proof.',
}
path=ROOT/'results/fresh_cubic7_dihedral_pell_cap_20261013.json'
path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps({'output':str(path),'status':out['status'],'height':str(alpha*X*X)}))
