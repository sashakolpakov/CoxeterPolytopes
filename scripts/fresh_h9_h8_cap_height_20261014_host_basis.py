#!/usr/bin/env python3
"""A free O_K basis for the maximal-height reflection module."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sympy as s
import fresh_h9_mixed_d4_coupling_20261014_review as A
from fresh_h9_h8_cap_height_20261014_integral_host import integer_pair
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_cap_height_20261014_host_basis'
HOST='results/fresh_h9_h8_cap_height_20261014_integral_host.json'
BASE='results/fresh_h9_h8_cap_height_20261014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def mm(X,Y):return[[sum((x*y for x,y in zip(row,col)),A.Z)for col in zip(*Y)]for row in X]
def run():
    d,h=[json.loads((ROOT/p).read_text())for p in(BASE,HOST)]
    for obj in(d,h):
        for p,v in obj['sources'].items():assert sha(p)==v,p
    Q,N=A.dm(d['form']),A.dm(h['normals']);phi=A.PHI;a=(phi-1)/2
    B=[[A.E(i==j)for j in range(10)]for i in range(9)]+[[a]+[A.Z]*8+[A.E(-1)]]
    R=mm(N,A.inverse(B));Rint=[[integer_pair(x)for x in row]for row in R]
    # Explicit inverse containment: the eight E8 simple-coordinate vectors,
    # one time vector, and the thirteenth cap are O-linear combinations.
    T=[[A.Z]*15 for _ in range(10)]
    wall_for_axis={1:1,2:9,3:2,4:3,5:4,6:5,7:6,8:7}
    for axis,wall in wall_for_axis.items():T[axis][wall]=A.E(-1)
    assert N[8][0]==1 and not N[8][-1]
    T[0][8]=A.ONE
    for axis,wall in wall_for_axis.items():T[0][wall]=N[8][axis]
    T[9][13]=A.ONE
    Tint=[[integer_pair(x)for x in row]for row in T]
    assert mm(T,N)==B and mm(R,B)==N
    G=[[A.pair(x,y,Q)for y in B]for x in B]
    assert G[0][0]==-2*phi and G[0][9]==G[9][0]==-1 and G[9][9]==2
    assert all(not G[i][j]for i in(0,9)for j in range(1,9))
    assert [row[1:9]for row in G[1:9]]==[row[1:9]for row in Q[1:9]]
    def expr(x):return s.Rational(str(x.a.a))+s.Rational(str(x.a.b))*s.sqrt(5)
    spatial=s.Matrix([[expr(x)for x in row[1:9]]for row in G[1:9]])
    assert spatial.det()==1
    determinant=s.simplify(s.Matrix([[expr(x)for x in row]for row in G]).det())
    assert determinant==-3-2*s.sqrt(5)
    assert s.expand(determinant*determinant.xreplace({s.sqrt(5):-s.sqrt(5)}))==-11
    paths=[BASE,HOST,f'scripts/{TAG}.py','scripts/fresh_h9_h8_cap_height_20261014_integral_host.py',
           'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py']
    return A.enc(dict(status='PASS_EXPLICIT_FREE_MAXHEIGHT_HOST_BASIS',sources={p:sha(p)for p in paths},
        basis=B,gram=G,generators_in_basis_integer_pairs=Rint,basis_in_generators_integer_pairs=Tint,
        decomposition='E8 orthogonal sum binary Gram[[-2phi,-1],[-1,2]]',determinant=['-3','-2'],determinant_field_norm=-11,
        enumeration_coordinates='r=x*e0+sum(y_i*e_i)+t*w, x,t,y_i in Z[phi]; w=((phi−1)/2)*e0−e9. Norm=−2phi*x²−2xt+2t²+E8[y].',
        scope='Both module inclusions are checked by literal integral coefficient matrices; this is a free basis of the15-generator reflection module, not merely a rational basis or a lattice saturation. It prepares arithmetic enumeration but does not enumerate additional roots or establish reflectivity.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['determinant_field_norm'])
