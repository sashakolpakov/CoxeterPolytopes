#!/usr/bin/env python3
"""Exact H4-squared terminal domain, duplicated9/10 contact, deletion debit.

Uses frozen octic physical arithmetic and the complete199 finite incidence;
this is not a new chamber or an all-cap obstruction.
"""
import argparse,json
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_escape_target_20261016 as E
import fresh_nonarithmetic_release30_66_hinge_20261017 as H
ROOT=E.ROOT;TAG='fresh_old9_terminal_identity_20261017'

def run():
    R,C=E.R,E.C
    C.K,C.Z,C.O=R.K,R.Z,R.O
    base=json.loads((ROOT/E.BASE).read_text());Q,N=R.m4(base['form']),R.m4(base['normals'])
    pair=lambda x,y:C.pair(x,y,Q)
    support=[0,1,2,3,5,6,7,8];phi=(R.O+R.BASIS[1])/2
    gram=[[pair(N[i],N[j])for j in support]for i in support]
    expected=[[2*R.O if i==j else R.Z for j in range(8)]for i in range(8)]
    for i,j,g in [(0,1,-phi),(1,2,-R.O),(2,3,-R.O),(4,5,-R.O),(5,6,-R.O),(6,7,-phi)]:expected[i][j]=expected[j][i]=g
    assert gram==expected and R.pd(gram)[0]
    # The second mandatory endpoint support has the same finite type.
    support2=[0,1,2,3,5,6,7,30]
    gram2=[[pair(N[i],N[j])for j in support2]for i in support2]
    assert gram2==expected
    difference=[x-y for x,y in zip(N[9],N[10])]
    coeff=C.mv(C.inverse(gram),[pair(difference,N[i])for i in support])
    reconstructed=[sum((a*N[i][k]for a,i in zip(coeff,support)),R.Z)for k in range(10)]
    assert reconstructed==difference
    assert pair(N[9],N[9])==pair(N[10],N[10])==2*R.O
    current,labels,ns,h,profiles,primary=H.setup()
    debits=[]
    for deleted in[[9],[9,10],[53],[12,47,52,64,65]]:
        d=H.deletion(primary,current,deleted)
        hs=[r for r in d['compact_faces_touched']if r['dimension']==8]
        debits.append(dict(deleted_walls=deleted,lost_finite_count=d['lost_finite_count'],lost_finite_supports=d['lost_finite_supports'],
          compact_faces_by_dimension=d['compact_faces_by_dimension'],compact_H8_certificates=[dict(support=r['support'],finite_vertices=len(r['finite_supports']))for r in hs],
          compact_S6_sides_touched=d['compact_S6_sides_touched']))
    assert [d['lost_finite_count']for d in debits]==[47,114,0,11]
    assert debits[0]['compact_H8_certificates']==[dict(support=[11],finite_vertices=41)]
    assert debits[1]['compact_H8_certificates']==debits[0]['compact_H8_certificates']
    assert debits[-1]['compact_faces_by_dimension']=={2:6}
    paths=[Path(__file__),Path(E.__file__),Path(R.__file__),Path(C.__file__),Path(H.__file__),ROOT/E.BASE,H.PARENT,H.BOUNDARY]
    return dict(status='PASS_EXACT_H4_SQUARED_TERMINAL_CONTACT_IDENTITY',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      support=support,second_endpoint_support=support2,finite_type='H4+H4',literal_bonds=[[0,1,5],[1,2,3],[2,3,3],[5,6,3],[6,7,3],[7,8,5]],
      field='Frozen physical octic field; coefficients in descending powers of its exact primitive element',relation='n9-n10=sum(coefficients[i]*n_support[i])',coefficients=R.enc(coeff),
      exact_coordinate_identities=10,native_squares=[2,2],deletion_debits=debits,
      theorem='Every spherical extension of the displayed H4+H4 support has the added root orthogonal to all eight roots, by the finite irreducible Coxeter classification. Consequently its contacts with native norm2 walls9 and10 are equal. Deleting9 alone cannot repair a failed10/9 Coxeter contact for this terminal duty.',
      scope='Finite terminal at S6+7+8 only; no restriction on walls missing that corner. The second H4+H4 support identifies the analogous rightness at S6+7+30. Deletion counts are losses of existing certificates, not a claim that a replacement cannot rebuild them. No new geometry, volume claim, or general obstruction after releasing support walls.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['exact_coordinate_identities'],[d['lost_finite_count']for d in out['deletion_debits']],flush=True)
