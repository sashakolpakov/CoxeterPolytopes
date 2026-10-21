#!/usr/bin/env python3
"""Two quartic caps meeting at pi3, from rank7/two-anchor proposals."""
import argparse,json
from fractions import Fraction as F
import fresh_h9_quartet_polar_second_20261014 as B
A,D,R,C,ROOT=B.A,B.D,B.P,B.C,B.ROOT
TAG='fresh_h9_joint_pair_20261015_twocaps'
BASE='results/fresh_h9_nine_contact_extension_20261015.json'
SCOUT='results/fresh_h9_joint_pair_20261015_scout_escape.json'
def run():
    d=B.read(BASE);A.RAD=A.A.Q(*d['field']['new_positive_radical_square']);C.Golden=A.L;C.ZERO=A.Z
    Q,N,H,E=D.dm(d['form']),D.dm(d['normals']),D.dv(d['center']),D.dv(d['strict_escape']['point'])
    ids=[0,1,2,3,5,6,7,10,29];cov=[A.mv(Q,n)for n in N]
    assert all(A.pair(N[i],N[i],Q)==2 for i in ids)
    v=C._full_rank_kernel([cov[i]for i in ids]);v=[x/v[9]for x in v]
    phi=A.L(A.A.Q(F(1,2),F(1,2)));rows=[]
    recipes=[({7:5,29:3},(F(30471,18764),F(7413,18764))),({3:5,5:3,7:5},(F(7236,4691),F(4059,4691)))]
    for contacts,beta in recipes:
        rhs=[-phi if contacts.get(i)==5 else A.L(-1)if contacts.get(i)==3 else A.Z for i in ids]
        p=A.mv(A.inv([cov[i]for i in ids]+[A.mv(Q,v)]),rhs+[A.Z]);assert not A.pair(p,v,Q)
        delta=(A.L(2)-A.pair(p,p,Q))/A.pair(v,v,Q);root=A.L(0,A.A.Q(*beta));assert root.sign()>0 and root*root==delta
        options=[[a+sgn*root*b for a,b in zip(p,v)]for sgn in(-1,1)]
        accepted=[]
        for n in options:
            assert A.pair(n,n,Q)==2
            vals=[A.pair(n,m,Q)for m in N]
            if A.pair(n,H,Q).sign()>=0 or A.pair(n,E,Q).sign()<=0:continue
            relations=[]
            for i,g in enumerate(vals):
                kind,m,c=R.relation(g,A.L(2),A.pair(N[i],N[i],Q))
                if kind=='forbidden':break
                relations.append(dict(wall=i,kind=kind,label=m,invariant=c))
            else:accepted.append(dict(normal=n,norm=A.L(2),pairings=vals,relations=relations,center_pairing=A.pair(n,H,Q),old_escape_pairing=A.pair(n,E,Q),contacts=[[i,contacts.get(i,2)]for i in ids],reconstruction=dict(kernel=v,particular=p,delta=delta,positive_square_root=root)))
        assert len(accepted)==1;rows+=accepted
    mutual=A.pair(rows[0]['normal'],rows[1]['normal'],Q);assert mutual==-1
    allN=N+[r['normal']for r in rows]
    for r in rows:
        n=r['normal'];t=A.pair(H,n,Q)/2;V=[a-t*b for a,b in zip(H,n)];vp=[A.pair(V,m,Q)for m in allN]
        assert A.pair(V,V,Q).sign()<0 and sum(x==0 for x in vp)==1 and all(x.sign()<=0 for x in vp)
        r['facet_witness']=dict(point=V,norm=A.pair(V,V,Q),pairings=vp)
    paths=[BASE,SCOUT,f'scripts/{TAG}.py','scripts/fresh_h9_quartet_polar_second_20261014.py',
           'scripts/fresh_h9_quartet_polar_20261014.py','scripts/fresh_h9_h8_disjoint_walk_20261014_exact.py',
           'scripts/fresh_h9_h8_disjoint_walk_20261014_next_frontier.py']
    return A.enc(dict(status='PASS_EXACT_TWO_QUARTIC_ESCAPE_CAPS_WITH_MUTUAL_PI3',sources={p:R.sha(p)for p in paths},field=d['field'],caps=rows,mutual_pairing=mutual,mutual_label=3,
       scope='Two exact norm2 caps compatible with all30old walls and each other, both cutting the currentescape. Actualfacet witnesses checked. Alternative to the polarcap; no fullgeometry or finitevolume claim in this localcertificate.'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
