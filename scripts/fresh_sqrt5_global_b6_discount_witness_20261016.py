#!/usr/bin/env python3
"""Exact compact local B6 quadrilateral missed by the global40 proposal."""
import argparse,json
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_adjacent_discount_witness_20261016 as A
K,F,Z,O=A.K,A.F,A.Z,A.O
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_global_b6_discount_witness_20261016'
SOURCE='results/fresh_sqrt5_global_b1_finite_first_20261016.json'
RULE='results/fresh_sqrt5_global_discount_20261016.json'

def run():
    prior=json.loads((ROOT/SOURCE).read_text());rec=prior['first_unpaid_singular'];row=rec['row'];n=row['nodes']
    assert n==8 and rec['index']==39 and row['white']==list(range(4,8))
    phi=K(F(1,2),F(1,2));G=[[K(2*(i==j))for j in range(n)]for i in range(n)];G[4][4]=K(4)
    bonds=[]
    for(i,j),m in zip(combinations(range(n),2),row['labels']):
        if m!=2:G[i][j]=G[j][i]=-{3:O,4:K(2),5:phi}[m];bonds.append((i,j,m))
    assert A.enc(G)==rec['scaled_gram']
    assert A.inertia(G)==[6,1,1]and A.inertia([[x.conj()for x in r]for r in G])==[7,0,1]
    ids=list(range(1,n));Q=A.sub(G,ids);Qi=A.inverse(Q)
    N=[A.mv(Qi,[G[j][i]for j in ids])for i in range(n)]
    assert[[A.pair(Q,a,b)for b in N]for a in N]==G
    H=A.mv(Qi,[-O]*7);hp=[A.pair(Q,H,v)for v in N];hq=A.pair(Q,H,H)
    assert hq.sign()<0 and max(x.sign()for x in hp)<0
    white=row['white'];assert A.inertia(A.sub(G,white))==[4,0,0]
    vertices=[]
    for a,b in product((0,2),(1,3)):
        S=white+[a,b];SG=A.sub(G,S);assert A.inertia(SG)==[6,0,0]
        c=A.mv(A.inverse(SG),[hp[i]for i in S])
        v=[H[j]-sum((d*N[i][j]for d,i in zip(c,S)),Z)for j in range(7)]
        ps=[A.pair(Q,v,r)for r in N]
        assert all(not ps[i]for i in S)and all(ps[i].sign()<0 for i in range(n)if i not in S)
        assert A.pair(Q,v,v).sign()<0 and A.pair(Q,v,H).sign()<0
        vertices.append(dict(black=[a,b],active=sorted(S),coordinates=[x.enc()for x in v],norm=A.pair(Q,v,v).enc(),products=[x.enc()for x in ps]))
    tests=[];accepted=[]
    for black in combinations(range(4),2):
        S=white+list(black);ker=A.kernel([A.mv(Q,N[i])for i in S]);assert len(ker)==1
        signs=[]
        for sign in(-1,1):
            v=[sign*x for x in ker[0]];ps=[A.pair(Q,v,r)for r in N]
            if max(x.sign()for x in ps)<=0:
                assert A.pair(Q,v,v).sign()<0 and A.pair(Q,v,H).sign()<0
                signs.append(sign);accepted.append(list(black))
        tests.append(dict(black=list(black),feasible_orientations=signs))
    assert sorted(accepted)==sorted([sorted(r['black'])for r in vertices])and len(accepted)==4
    # Explicit common orthogonal spectators: 8 A1 + 2 I2(5) + A3.
    S=[[K(2*(i==j))for j in range(15)]for i in range(15)]
    sbonds=[(8,9,5),(10,11,5),(12,13,3),(13,14,3)]
    for i,j,m in sbonds:S[i][j]=S[j][i]=-(phi if m==5 else O)
    assert A.inertia(S)==A.inertia([[x.conj()for x in r]for r in S])==[15,0,0]
    hs=[-x/K(1024)for x in A.mv(A.inverse(S),[O]*15)];fullnorm=hq+A.pair(S,hs,hs)
    assert fullnorm.sign()<0
    big=[[Z]*23 for _ in range(23)]
    for i in range(n):
        for j in range(n):big[i][j]=G[i][j]
    for i in range(15):
        for j in range(15):big[n+i][n+j]=S[i][j]
    assert A.inertia(big)==[21,1,1]and A.inertia([[x.conj()for x in r]for r in big])==[22,0,1]
    assert all((2*big[i][j]/big[j][j]).integral()for i in range(23)for j in range(23))
    # Quotient the one rational relation; the Z[phi] span of these roots is a
    # reflection-stable full lattice in the 22-dimensional admissible K space.
    import fresh_sqrt5_global_discount_20261016 as B
    B.install();M=B.M;X=B.X
    edges=bonds+[(n+i,n+j,m)for i,j,m in sbonds];fullwhite=list(range(4,23))
    compiled=[X.compile_corner(fullwhite+list(pair),edges,pair,[None])for pair in product((0,2),(1,3))]
    profile,tables=X.add_bonus(compiled,(0,)*8)
    for table in tables:profile=M.plus(profile,table[None])
    expected={-1:F(1),1:F(-1),14894:F(1,10),21728:F(1,10)}
    assert profile==expected
    assert profile=={int(j):F(z)for j,z in rec['payment']['worst_profile']}
    doc=json.loads((ROOT/RULE).read_text());weights=list(map(F,doc['certificate']['primal']))
    value=M.value(profile,weights);assert weights[14894]==weights[21728]==0 and value==1-weights[1]<1
    metadata=doc.get('variable_order')
    paths=[f'scripts/{TAG}.py',f'scripts/{A.TAG}.py',f'results/{A.TAG}.json',SOURCE,RULE,f'scripts/{B.TAG}.py']
    return dict(status='PASS_EXACT_GLOBAL_B6_LOCAL_UNDERPAYMENT',
      sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},field='a+b*sqrt(5)',
      original_minimal_class=3249,finite_first_class=39,row=row,scaled_gram=A.enc(G),basis_indices=ids,
      basis_gram=A.enc(Q),root_coordinates=A.enc(N),kernel=[[x.enc()for x in v]for v in A.kernel(G)],
      core_inertias=dict(physical=[6,1,1],conjugate=[7,0,1]),strict_center=[x.enc()for x in H],
      center_norm=hq.enc(),center_products=[x.enc()for x in hp],vertices=vertices,complete_face_cone=tests,
      spectator_forest=['A1']*8+['I2(5)']*2+['A3'],spectator_gram=A.enc(S),full_center_norm=fullnorm.enc(),
      full_inertias=dict(physical=[21,1,1],conjugate=[22,0,1]),integral_cartan_entries=529,
      literal_rank21_corner_profile=[[j,str(z)]for j,z in sorted(profile.items())],
      coefficient_values={str(j):str(weights[j])for j in(1,14894,21728)},payment=str(value),
      necessary_payment_inequality='(w14894+w21728)/10 >= delta_adjacent3',
      scope='Actual compact quadrilateral in a local six-dimensional arithmetic-admissible integral-reflection chamber, with a common rank15 spherical spectator forest giving rank21 corner data. All six candidate ray pairs and both orientations are checked, with exactly four finite feasible future rays. This disproves full coverage of the frozen global40 proposal, but does not invalidate its explicitly restricted certificate or the separate compact arithmetic dimension21 proof candidate. No compact ambient21 realization and no new polytope or dimension bound is claimed.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['payment'],flush=True)
