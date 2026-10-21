#!/usr/bin/env python3
"""Independent reconstruction of the nine-contact quadratic cap in the existing quartic field."""
import argparse,json
from pathlib import Path
from collections import Counter
from hashlib import sha256
import fresh_h9_quartet_polar_extension_20261014_review as R

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_nine_contact_cap_20261015_review'
BASE='results/fresh_h9_fourth_escape_extension_20261015.json'
PRIMARY='results/fresh_h9_nine_contact_cap_20261015.json'
FACES='results/fresh_h9_vertex_free_edges_20261015.json'
K,Z,O,C=R.K,R.Z,R.O,R.C
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):return list(map(str,x.to_list()))

def run():
    paths=[BASE,PRIMARY]
    b,d=[json.loads((ROOT/p).read_text())for p in paths]
    for doc in(b,d):
        for p,h in doc.get('sources',{}).items():assert sha(p)==h,p
    assert len(d['caps'])==1
    Q,N,H,E=R.dm(b['form']),R.dm(b['normals']),R.dv(b['center']),R.dv(b['strict_escape']['point'])
    assert len(N)==29
    pair=lambda x,y:C.pair(x,y,Q)
    cov=[C.mv(Q,n)for n in N];norms=[pair(n,n)for n in N]
    row=d['caps'][0];contacts=dict(row['finite_contacts']);phi=(O+R.BASIS[1])/2
    targets={2:Z,3:-O,5:-phi};basis=list(reversed(list(contacts)))
    assert len(basis)==9 and len(C.eliminate([cov[j]for j in basis])[1])==9
    v=C.null_one([cov[j]for j in basis]);assert v is not None and v[9]
    v=[a/v[9]for a in v]
    p=C.mv(C.inverse([cov[i]for i in basis]+[C.mv(Q,v)]),[targets[contacts[i]]for i in basis]+[Z])
    assert not pair(p,v)
    delta=(2*O-pair(p,p))/pair(v,v)
    root=(5*O+7*R.BASIS[1])*R.BASIS[2]/11
    assert root*root==delta and R.sg(root)>0
    n=[a-root*b for a,b in zip(p,v)]
    assert n==R.dv(row['normal'])and pair(n,n)==2*O
    vals=[C.dot(n,cc)for cc in cov];assert vals==R.dv(row['pairings'])
    assert all(vals[i]==targets[m]for i,m in contacts.items())
    counts=Counter()
    for i,g in enumerate(vals):
        kind,label,u=R.relation(g,2*O,norms[i]);assert kind!='forbidden'
        saved=row['relations'][i];assert(kind,label)==(saved['kind'],saved['label'])
        assert u==(None if saved['invariant']is None else R.dec(saved['invariant']))
        counts[kind+(str(label)if label else'')]+=1
    assert pair(n,H)==R.dec(row['center_pairing'])and R.sg(pair(n,H))<0
    assert pair(n,E)==R.dec(row['old_escape_pairing'])and R.sg(pair(n,E))>0
    S=[0,1,2,5,6,8,10,28];saved=row['finite_endpoint']
    assert S==saved['support']and len(S)==8
    V=C.null_one([C.mv(Q,n)]+[cov[i]for i in reversed(S)])
    if R.sg(pair(V,H))>0:V=[-v for v in V]
    sv=R.dv(saved['point']);assert R.canon(V)==R.canon(sv)
    assert pair(sv,sv)==R.dec(saved['norm'])and R.sg(pair(sv,sv))<0 and R.sg(pair(sv,H))<0
    support=[N[i]for i in S]+[n];G=[[pair(a,b)for b in support]for a in support]
    assert R.pd(G)and all(R.sg(x)>=0 for rr in C.inverse(G)for x in rr)
    vp=[C.dot(sv,cc)for cc in cov]
    assert vp==R.dv(saved['pairings'])and all(R.sg(g)==(0 if i in S else-1)for i,g in enumerate(vp))
    assert not pair(sv,n)
    assert all(R.canon(sv)!=R.canon(R.dv(v['point']))for v in b['finite_vertices'])
    octic=json.loads((ROOT/'results/fresh_h9_octic_escape_caps_20261015.json').read_text())
    common=R.dv(octic['common_strict_escape']['point'])
    assert pair(n,common)==R.dec(row['octic_common_escape_pairing'])and R.sg(pair(n,common))>0
    allN=N+[n];allcov=cov+[C.mv(Q,n)]
    assert R.sg(pair(H,H))<0 and all(R.sg(C.dot(H,cc))<0 for cc in allcov)
    for i,m in enumerate(allN):
        p=[h-C.dot(H,allcov[i])/pair(m,m)*x for h,x in zip(H,m)]
        assert R.sg(pair(p,p))<0 and R.sg(pair(p,H))<0
        assert all(R.sg(C.dot(p,cc))==(0 if i==j else-1)for j,cc in enumerate(allcov))
    paths += ['results/fresh_h9_octic_escape_caps_20261015.json',f'scripts/{TAG}.py','scripts/fresh_h9_nine_contact_cap_20261015.py',
              'scripts/fresh_h9_quartet_polar_extension_20261014_review.py',
              'scripts/fresh_h9_h8_all_high_20261014_review.py',
              'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return dict(status='PASS_INDEPENDENT_NINE_CONTACT_QUARTIC_CAP',sources={p:sha(p)for p in paths},
        arithmetic='Independent SymPy ANP/192bit outward signs, reversed contact solve and physical nullspace; no producer arithmetic imports.',
        independent_reconstruction_basis=basis,normal=[enc(x)for x in n],norm=enc(pair(n,n)),
        new_pairs=29,pair_counts=dict(counts),actual_facets=30,facet_sign_checks=900,
        quadratic_radicand=enc(delta),positive_square_root=enc(root),endpoint_support=S+[29],endpoint=[enc(x)for x in sv],endpoint_norm=enc(pair(sv,sv)),
        endpoint_sign_checks=30,endpoint_distinct_from_old88=True,old_escape_pairing=enc(pair(n,E)),
        scope='Exact reversed nine-contact solve and quadratic root identity reconstruct the norm2 wall in K. All29 pair relations and30 actualfacet witnesses pass; the wall cuts both recorded escapes and creates a finite endpoint. Global30-wall cone and newescape require the separate fullreview.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved independent replacement-cap replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'])
