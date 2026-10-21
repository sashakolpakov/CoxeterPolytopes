#!/usr/bin/env python3
"""Produce a rational escape near a null end of the unchanged S61 face.

Numerical scouting only; the accepted rational point needs a separate exact
all-wall replay. The squared norm of each rounded proposal is exact here.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import mpmath as mp
from fresh_h20_coupled_roles_20261014 import load

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_coupled_escape_20261014_scout'
FILES=['results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json',
       'results/fresh_h20_e7_bridge_face_20261014_exact.json',
       'results/fresh_h20_bridge_edge_closure_20261014.json',
       'results/fresh_h20_coupled_roles_20261014_exact.json']

def run():
    old,bridge,edge,new=[json.loads((ROOT/p).read_text()) for p in FILES]
    N,Q,G,H,oldE,dot=load()
    def f(x):
        a=F(x);return mp.mpf(a.numerator)/a.denominator
    rad=[mp.sqrt(x) for x in old['field']['base_square_root_basis']]
    def v0(x):return sum(f(c)*r for c,r in zip(x,rad))
    D=mp.sqrt(v0(old['field']['discriminant']))
    def v1(x):return v0(x[0])+D*v0(x[1])
    K=mp.sqrt(v1(bridge['field']['new_positive_sqrt_radical']))
    def v2(x):return v1(x[0])+K*v1(x[1])
    R=mp.matrix([v1(x) for x in new['cap_R']['normal']])
    qR=v1(new['cap_R']['norm'])
    assert abs(dot(R,R)-qR)<mp.mpf('1e-65')*qR
    N.append(R*mp.sqrt(2/qR))
    frame=[mp.matrix([v1(x) for x in new['hex_basis'][key]])/mp.sqrt(v1(q))
           for key,q in zip(('C','U','V'),new['hex_basis']['squared_lengths'])]
    for coeff in new['hexagon']['frame_coefficients']:
        N.append(sum((v0(a)*v for a,v in zip(coeff,frame)),mp.zeros(21,1)))
    assert len(N)==107 and all(abs(dot(v,v)-2)<mp.mpf('1e-65') for v in N)
    S=edge['recipe']['support']
    X=mp.matrix([v2(x) for x in edge['recipe']['X']]);A=v2(edge['recipe']['A'])
    W=mp.matrix([v2(x) for x in edge['recipe']['W']]);B=v2(edge['recipe']['W_norm'])
    time=X/mp.sqrt(A);space=N[95]/mp.sqrt(2);side=W/mp.sqrt(B)
    assert max(abs(dot(a,b)-(int(i==j)*(1 if i else -1)))
               for i,a in enumerate((time,space,side)) for j,b in enumerate((time,space,side)))<mp.mpf('1e-65')
    J=mp.matrix([[dot(N[i],N[j]) for j in S] for i in S])**-1
    Y=-sum((a*N[i] for a,i in zip(J*mp.ones(len(S),1),S)),mp.zeros(21,1))
    y2=dot(Y,Y);H=H/mp.sqrt(-dot(H,H))
    exactQ=[[F(x) for x in row] for row in old['form']]
    tries=0
    # First try the midpoint direction opposite the old cutter; then scan
    # other null directions only if that one fails a nonsupport inequality.
    angles=[mp.pi]+[mp.pi+mp.pi*j/256 for j in range(-127,128) if j]
    for ai,theta in enumerate(angles):
        nu=time+mp.cos(theta)*space+mp.sin(theta)*side
        if any(dot(nu,n)>=-mp.mpf('1e-55') for i,n in enumerate(N) if i not in S):continue
        tau=y2/(-4*dot(nu,H))
        assert tau>0
        for power in range(0,22):
            eps=mp.mpf(10)**(-power)
            proposal=nu/eps+Y+eps*tau*H
            if dot(proposal,proposal)<=0 or any(dot(proposal,n)>=0 for n in N):continue
            for digits in (6,9,12,15,18,24,30):
                tries+=1;den=10**digits
                point=[F(int(mp.nint(x*den)),den) for x in proposal]
                norm=sum(point[i]*exactQ[i][j]*point[j] for i in range(21) for j in range(21))
                test=mp.matrix([f(x) for x in point])
                pairings=[dot(test,n) for n in N]
                if norm>0 and all(x<0 for x in pairings) and dot(test,H)<0:
                    paths=FILES+[f'scripts/{TAG}.py','scripts/fresh_h20_coupled_roles_20261014.py']
                    return dict(status='RATIONAL_ESCAPE_NUMERICALLY_PROPOSED_EXACT_NORM',
                        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                        support=S,angle_index=ai,angle_radians=mp.nstr(theta,40),
                        epsilon_power=power,rounding_digits=digits,rounding_attempts=tries,
                        rational_point=list(map(str,point)),exact_norm=str(norm),
                        normalized_pairings=[mp.nstr(x,40) for x in pairings],
                        center_pairing=mp.nstr(dot(test,H),40),
                        scope='The norm is exact rational. All107 pair signs and orientation here are100-decimal-digit numerical proposals, requiring separate exact arithmetic/interval acceptance.')
    raise AssertionError('No rounded escape found')

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    result=run()
    out=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result==json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['angle_index'],result['epsilon_power'],result['rounding_digits'],result['exact_norm'])
