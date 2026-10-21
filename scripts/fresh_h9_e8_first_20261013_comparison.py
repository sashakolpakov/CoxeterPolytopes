#!/usr/bin/env python3
"""Exact comparison with D9 and its standard maximal overclosure."""
import argparse,hashlib,json
from fresh_h9_e8_first_20261013 import ROOT,PREFIX,BETA,G,K,Z,O,T,ALPHA,pair,inverse,mv,enc,menc

OUT=ROOT/'results'/f'{PREFIX}_comparison.json'
def run(verify=False):
    src=ROOT/'results'/f'{PREFIX}.json';d=json.loads(src.read_text())
    tri=ROOT/'results/fresh_h9_tripled_d9_20261013_e8_triage.json';t=json.loads(tri.read_text())
    oldp=ROOT/'results/fresh_h9_tripled_d9_20261013_edge1.json';old=json.loads(oldp.read_text())
    r=tuple(map(K.from_json,d['chosen_root']));q=pair(r,G,r)
    U=[tuple(O if i==j else Z for i in range(10))for j in range(10)]
    DB=[U[0]]+[tuple(a-b for a,b in zip(U[i],U[i+1]))for i in range(1,9)]+[tuple(a+b for a,b in zip(U[8],U[9]))]
    MB=tuple(tuple(map(K.from_json,s))for s in d['lattice_basis']);MI=inverse(tuple(zip(*MB)))
    def Rd(s):return tuple(a-2*pair(r,G,s)/q*b for a,b in zip(s,r))
    def d9(s):return all(x.is_integral for x in s)and(sum(s[1:],Z)/2).is_integral
    assert all(all(x.is_integral for x in Rd(s))for s in U)
    assert all(all(x.is_integral for x in mv(MI,Rd(s)))for s in MB)
    failures=[]
    for i,s in enumerate(DB):
        if not d9(Rd(s)):failures.append({'basis_index':i,'vector':enc(s),'image':enc(Rd(s)),'spatial_half_sum':(sum(Rd(s)[1:],Z)/2).to_json()})
    assert failures and failures[0]['basis_index']==0
    # Coordinate permutation preserves D9 and O^9. Move the new scalar
    # mirror to the old first axis to compare the certified edge intervals.
    rr=list(r);rr[1],rr[9]=rr[9],rr[1];rr=tuple(rr)
    p=tuple(map(K.from_json,old['endpoint']));cut=pair(rr,G,p)
    assert cut>0
    assert K.from_json(d['upper_parameter'])<K.from_json(old['upper_parameter'])
    # Primitive endpoint A1 norms lie in different square classes, so a
    # K-isometry cannot identify these two complete pointed sectors.
    assert q==3*(2+T)and q.norm()==162
    assert 12**2<162<13**2
    out={'status':'PASS_EXPLICIT_OVERGROUP_MIRROR_COMPARISON',
        'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(src,tri,oldp)},
        'shared_root':enc(r),'norm':q.to_json(),
        'preserves_e8_lattice':True,'preserves_standard_maximal_lattice':True,
        'd9_failures':failures,'axis_permuted_root':enc(rr),
        'old_d9_first_endpoint':enc(p),'strict_cut_pairing':cut.to_json(),
        'new_first_parameter':d['upper_parameter'],'old_d9_first_parameter':old['upper_parameter'],
        'initial_endpoint_A1_norm_squareclasses':['1','3*(2+T)'],
        'scope':'The standard O^9 maximal overclosure has additional mirrors that cut the old full D9 chamber. The D9 nonreflectivity theorem therefore does not transfer merely by saturation. The displayed E8 and O^9 lattices share this first binary reflection, but this neither identifies their full groups nor excludes another indefinite conjugacy. The two E8+A1 corners of the new first edge are not K-isometric as pointed full root sectors because the distinguished A1 norm ratio has nonsquare rational norm162.'}
    if verify:assert out==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'D9 basis failures',len(failures),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
