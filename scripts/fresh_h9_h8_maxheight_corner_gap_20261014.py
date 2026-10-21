#!/usr/bin/env python3
"""All-parameter obstruction for caps right to a source vertex's eight walls."""
import argparse
import json
from hashlib import sha256
from pathlib import Path
import sympy as s
import fresh_h9_mixed_d4_coupling_20261014_review as A
from fresh_h9_h8_all_high_20261014_review import K,decode,pair as independent_pair
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_h8_maxheight_corner_gap_20261014'
BASE='results/fresh_h9_h8_cap_height_20261014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    base=json.loads((ROOT/BASE).read_text())
    for p,h in base['sources'].items():assert sha(p)==h,p
    Q,N=A.dm(base['form']),A.dm(base['normals'])
    E=A.dv(base['strict_escape']['point']);d=Q[-1][-1]
    assert E[-1]==-5 and A.pair(E,E,Q).sign()>0
    Qk=[[decode(x)for x in row]for row in base['form']]
    Nk=[[decode(x)for x in row]for row in base['normals']]
    Ek=[decode(x)for x in base['strict_escape']['point']];dk=Qk[-1][-1]
    records=[]
    for i,row in enumerate(base['complete_finite_vertices'][:41]):
        v=A.dv(row['point']);assert not v[-1]
        contacts=[-A.pair(v,n,Q)for n in N[:11]]
        nonzero=[x for x in contacts if x];assert len(nonzero)==3
        minimum=nonzero[0]
        for x in nonzero:
            assert x.sign()>0
            if(x-minimum).sign()<0:minimum=x
        v=[x/minimum for x in v]
        contacts=[x/minimum for x in contacts]
        assert any(x==1 for x in contacts)
        assert all(not x or(x-1).sign()>=0 for x in contacts)
        b=-A.pair(v,v,Q);g=A.pair(v,E,Q)
        margin=g*g-25*d*(2+b)
        assert b.sign()>0 and g.sign()<0 and margin.sign()>0
        # Reconstruct the normalized literal point independently and compare
        # the exact polynomial margin in a separate algebraic-field engine.
        vk=[decode(x)for x in row['point']]
        h=[-independent_pair(vk,n,Qk)for n in Nk[:11]]
        hk=decode(A.enc(minimum))
        assert any(x==hk for x in h)
        for x in h:
            if x:assert s.simplify(K.to_sympy(x-hk)).is_nonnegative is True
        vk=[x/hk for x in vk]
        bk=-independent_pair(vk,vk,Qk);gk=independent_pair(vk,Ek,Qk)
        mk=gk*gk-25*dk*(2*K.one+bk)
        assert mk==decode(A.enc(margin))
        assert s.simplify(K.to_sympy(mk)).is_positive is True
        assert s.simplify(K.to_sympy(gk)).is_negative is True
        records.append(dict(source_vertex=i,source_support=[j for j,x in enumerate(contacts)if not x],
                            normalized_point=v,contacts=contacts,norm_magnitude=b,
                            escape_source_product=g,squared_strict_margin=margin))
    assert len(records)==41
    paths=[BASE,f'scripts/{TAG}.py','scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',
           'scripts/fresh_h9_h8_all_high_20261014_review.py',
           'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return A.enc(dict(status='PASS_ALL_MAXHEIGHT_SOURCE_CORNER_RIGHT_ENDPOINT_CAPS_MISS_ESCAPE',
        sources={p:sha(p)for p in paths},source_vertex_count=41,records=records,
        proof='For any source vertex normalized to minimum nonzero contact1, let b=-v² and g=v.E_source<0. Every norm2forward cap right to its eight source walls has r=(s*v,-sqrt((2+b*s²)/d)). Coxeter contacts force s>=1. All41 exact margins g²-25D(2+b)>0, so r.E=s*g+5sqrt(d(2+b*s²))<=s*(g+5sqrt(d(2+b)))<0. No angle cutoff or parameter sampling.',
        scope='All real scales and all source-corner directions for this right-support family; the complete41-source-vertex census is inherited. Every scalar margin is checked in separate Fraction and SymPy algebraic-field arithmetic. Does not treat nonright support profiles, arbitrary caps or general later edge supports. A failure to cut this one escape does not exclude adding these caps as intermediate steps.'))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['source_vertex_count'])
