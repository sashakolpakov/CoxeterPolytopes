#!/usr/bin/env python3
"""Exact endpoint test for every unresolved numerical B2 closure recipe.

This does not certify the floating-point census or its rejected candidates.
It reconstructs each saved unresolved recipe in its exact quadratic tower.
"""
from fractions import Fraction as F
from collections import Counter
import argparse, hashlib, json, math
import fresh_h20_b2_edge_sections_20261014_exact as X
P=X.P; ROOT=X.ROOT; TAG=X.TAG
SUMMARY=f'results/{TAG}_scout.json'

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()

def run():
    scout=json.loads((ROOT/SUMMARY).read_text())
    for p,h in scout['sources'].items():assert sha(p)==h,p
    assert scout['all5048_families_complete']
    assert sha(f'results/{TAG}_scout.jsonl.gz')==scout['archive_sha256']
    raw=json.loads((ROOT/P.RAW).read_text())['coxeter_input']
    Q=[[-int(x)for x in r]for r in raw['basis_gram']]
    W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']]
    C=[P.mv(Q,w)for w in W];G=[[P.dot(w,c)for c in C]for w in W]
    h=[sum(w[j]for w in W)for j in range(22)];g=math.gcd(*h);h=[F(x,g)for x in h];hc=P.mv(Q,h)
    walk=json.loads((ROOT/P.WALK).read_text())['corner_orbits'];rows=[];counts=Counter()
    for r in scout['survivors']:
        fam=scout['families'][r['family']];corner=walk[fam['corner']];edge=corner['edge_orbits'][fam['edge_index']]
        E=fam['edge_walls'];a=fam['anchor'];b=next(i for i in fam['component']if i!=a);l=r['closure_wall']
        assert E==edge['edge_walls'] and l not in E and G[a][a]==G[b][b]==2 and G[a][b]==-1
        assert all(G[i][j]==0 for i in[a,b]for j in E if j not in[a,b])
        ga,gb,d=G[l][a],G[l][b],G[l][l]
        matches=[v for v in X.contacts(ga,gb,d)if v[0]==r['contact_type']];assert len(matches)==1
        kind,base,x=matches[0]
        left=list(map(F,corner['primitive_ray']));right=list(map(F,edge['endpoint_primitive_ray']))
        left=P.scale(-1/P.dot(left,hc),left);right=P.scale(-1/P.dot(right,hc),right)
        B=[[P.pair(v,w,Q)for w in[left,right]]for v in[left,right]]
        p=[P.dot(v,C[l])for v in[left,right]];pp=P.dot(p,p);assert pp>0
        target=X.add(x,(-F(2*ga+gb,3),F(0)));v0=[v/pp for v in p];y=[p[1],-p[0]]
        aa=P.dot(y,P.mv(B,y));bb=X.scale(target,2*P.dot(v0,P.mv(B,y)))
        dd=X.add(X.scale(X.mul(target,target,base),P.dot(v0,P.mv(B,v0))),(-F(1,3),F(0)))
        D=X.add(X.mul(bb,bb,base),X.scale(dd,-4*aa));branch=r['quadratic_branch']
        if aa:
            assert X.ksign(D,base)>=0 and branch in[-1,1]
            field=X.Field(base,D);root=(X.scale(bb,-1/(2*aa)),(F(branch,2)/aa,F(0)))
        else:
            assert X.ksign(bb,base) and branch==0
            field=X.Field(base,X.ZERO);root=field.k(X.scale(X.divide(dd,bb,base),-1))
        zc=[field.plus(field.k(X.scale(target,v0[i])),field.scaled(root,y[i]))for i in range(2)]
        vals=[field.rational_dot(B[i],zc)for i in range(2)]
        assert field.equal(field.dot(zc,vals),field.c(F(1,3)))
        assert field.equal(field.rational_dot(p,zc),field.k(target))
        signs=[field.sign(v)for v in vals]
        assert 0 in signs and signs!=[0,0],(r,signs)
        counts['BOUNDARY_ROOT']+=1
        if X.ksign(D,base)==0:counts['DOUBLE_ROOT']+=1
        rows.append({'scout_candidate':r['index'],'family':r['family'],'closure_wall':l,
            'contact_type':kind,'base_radicand':base,'outer_radicand':D,
            'quadratic_coefficients':[aa,bb,dd],'quadratic_branch':branch,
            'edge_plane_coefficients':zc,'endpoint_pairings':vals,'endpoint_signs':signs,
            'zero_endpoint':signs.index(0),'status':'EXACT_OLD_ENDPOINT_NOT_TRANSVERSE_INTERIOR'})
    sources=[SUMMARY,f'results/{TAG}_scout.jsonl.gz',P.RAW,P.WALK,
        f'scripts/{TAG}_remaining.py',f'scripts/{TAG}_exact.py','scripts/fresh_h20_cusp_sections_20261013.py']
    return X.enc({'status':'PASS_ALL_SAVED_B2_NUMERICAL_REMAINDERS_ARE_EXACT_BOUNDARY_ROOTS',
        'sources':{p:sha(p)for p in sources},'counts':dict(counts),'cases':rows,
        'scope':'Every saved unresolved numerical recipe is exactly at an old endpoint, not a transverse interior section. This does not exactify the numerical census, its deduplication, or its other ridge rejections.'})

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();z=run()
    out=ROOT/f'results/{TAG}_remaining.json'
    if args.verify:assert z==json.loads(out.read_text())
    else:out.write_text(json.dumps(z,indent=2)+'\n')
    print(z['status']);print(z['counts'])
