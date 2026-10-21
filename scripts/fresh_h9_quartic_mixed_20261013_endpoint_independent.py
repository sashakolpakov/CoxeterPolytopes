#!/usr/bin/env python3
"""Fraction-only replay of Pell identity, roots, angles, endpoint, and escapes."""
from fresh_h9_quartic_mixed_20261013_target_independent import *

SOURCE=ROOT/'results/fresh_h9_quartic_mixed_20261013_endpoint.json'
OUT=ROOT/'results/fresh_h9_quartic_mixed_20261013_endpoint_independent.json'


def pd(matrix):
    rows=[list(row)for row in matrix]
    while rows:
        pivot=rows[0][0]
        assert sign(pivot)>0
        rows=[[sub(rows[i][j],divide(mul(rows[i][0],rows[0][j]),pivot))
               for j in range(1,len(rows))]for i in range(1,len(rows))]


def replay_endpoint():
    started=time.monotonic();data=json.loads(SOURCE.read_text());host=Host(2)
    pell=data['relative_pell'];y,s,x,t=(dec(pell[k])for k in('y','s','time','inner_covector'))
    assert sub(sq(y),scale(mul(ALPHA,sq(s)),6))==O
    assert x==scale(s,3)and t==scale(sub(y,O),Q(1,2))
    assert all(map(integral,(x,t)))and all(sign(z)>0 for z in(y,s,x,t))
    r=(x,*(add(mul(row[2],t),row[3])for row in FI),*(Z,)*5)
    new=tuple(vec(row)for row in data['new_roots'])
    assert set(new)==orbit(r)and len(new)==4
    source=json.loads((ROOT/data['source']).read_text());oldwalls=tuple(vec(r)for r in source['all_walls'])
    walls=tuple(vec(r)for r in data['all_walls']);assert walls==oldwalls+new
    assert divide(mul(ALPHA,sq(x)),scale(O,2))==dec(pell['height'])
    geometry=geometry_check(data,2);face=geometry_check(data['forced_H5_face'],1)
    assert tuple(vec(r)for r in data['forced_H5_face']['all_walls'])==tuple((*r[:5],r[9])for r in walls if all(z==Z for z in r[5:9]))
    endpoint=data['finite_endpoint'];v=vec(endpoint['vector'])
    c=divide(mul(ALPHA,x),scale(y,3));expected=(O,*(mul(row[2],c)for row in FI),*(Z,)*5)
    assert v==expected
    norm=host.pair(v,v);assert norm==dec(endpoint['norm'])==neg(divide(ALPHA,sq(y)))and sign(norm)<0
    center=vec(data['strict_interior']);assert host.pair(v,center)==dec(endpoint['center_pairing'])and sign(host.pair(v,center))<0
    values=tuple(host.pair(w,v)for w in walls);assert values==vec(endpoint['wall_pairings'])and all(sign(z)<=0 for z in values)
    active=[i for i,z in enumerate(values)if z==Z];assert active==endpoint['active_facets']and len(active)==9
    gram=tuple(tuple(host.pair(walls[i],walls[j])for j in active)for i in active);pd(gram)
    oldv=vec(source['forced_H5_face']['geometry']['vector'])
    cut=Host(1).pair((*r[:5],r[9]),oldv)
    assert cut==dec(data['previous_H5_escape_cut'])and sign(cut)>0
    profiles=set();alphabets={dec(row['norm']):row for row in data['finite_endpoint_alphabets']}
    for q in POOL:
        ws=sorted({r for label in ALLOWED if (r:=field_sqrt(scale(mul(q,label),Q(1,2))))is not None and integral(r)})
        zs=sorted({r for label in ALLOWED if (r:=field_sqrt(scale(mul(q,label),Q(1,4))))is not None and integral(r)})
        assert set(ws)==set(vec(alphabets[q]['F4_covectors']))and set(zs)==set(vec(alphabets[q]['scalar_coordinates']))
        for a,b,c,z in product(ws,ws,ws,zs):
            qp=sub(sub(sub(q,scale(add(add(sq(a),mul(a,b)),sq(b)),Q(2,3))),scale(sq(c),Q(1,2))),sq(z))
            if any(sign(qp,e)<=0 for e in range(4)):continue
            shift=add(add(scale(mul(SQRT2,a),2),scale(mul(SQRT2,b),4)),scale(c,3))
            profiles.add((q,(a,b,c,z),qp,shift))
    saved={(dec(row['norm']),vec(row['fixed_covectors']),dec(row['orthogonal_norm']),dec(row['completed_square_shift']))for row in data['finite_endpoint_profiles']}
    assert saved==profiles and len(profiles)==30
    forged=list(new[0]);forged[0]=add(forged[0],O);assert host.pair(tuple(forged),tuple(forged))!=scale(O,2)
    return {'status':'PASS_INDEPENDENT_PELL_ENDPOINT_AND_BOTH_CONES','source':str(SOURCE.relative_to(ROOT)),
            'pell_identity':True,'new_orbit_roots':4,'H9':geometry,'H5':face,
            'finite_endpoint_active_rank':9,'finite_endpoint_all_inequalities':len(walls),
            'finite_endpoint_norm':endpoint['norm'],'finite_label_profiles':len(profiles),
            'forged_root_rejected':True,'seconds':round(time.monotonic()-started,3),
            'scope':'No producer or PARI import. Exact Fraction replay of the four explicit roots and their actual finite endpoint, all cone root/angle gates, strict interior, and both surviving spacelike escapes. The finite-profile census is necessary endpoint data, not a full root census.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    result=replay_endpoint()
    if args.verify:
        saved=json.loads(OUT.read_text());saved.pop('seconds');result.pop('seconds');assert saved==result
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
