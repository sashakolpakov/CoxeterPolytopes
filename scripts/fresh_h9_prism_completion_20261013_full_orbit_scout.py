#!/usr/bin/env python3
"""Full16-seed diagram action, rather than the coordinate stabilizer."""
import json,sys,itertools
from fractions import Fraction
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scripts import fresh_h9_prism_completion_20261013_scout as S
V=S.V

def run():
    d=S.setup();G=np.array(d['gram']);old=np.array(d['normals']);U=[0,2,1,3,4,6,5,7,8,9]
    R=np.vstack([old,old[90:,U]]);names=d['names']+['mirror'+n[0].upper()+n[1:]for n in d['names'][90:]];cov=R@G
    assert names[:16]==['seed'+str(i)for i in range(16)]
    basis=list(range(8))+[8,12];assert np.max(abs(R[basis]-np.eye(10)))<1e-10
    source=S.PREFIX+'_caps.json';recipe=V.read(source);y=np.array([float(Fraction(x['a']))for x in recipe['escape']]);h=np.array(d['h'])
    target2=np.array(json.load(open('/tmp/fresh_h9_prism_all_orbit_geometry.json'))['point'])
    oldkeys={tuple(np.round(r,6)):i for i,r in enumerate(R)};seen={};rows=[];families=[];stats={'group_elements':0,'new_individual_images':0,'all102_compatible':0,'at_most_two_optional_conflicts':0}
    oldgram=R[:16]@G@R[:16].T
    for p,q,bs,uv in itertools.product(itertools.permutations(range(4)),itertools.permutations(range(4)),range(2),range(2)):
        perm=[]
        for i in range(16):
            layer=i//8;block=(i%8)//4;pos=i%4;j=(p if block==0 else q)[pos]
            perm.append((layer^uv)*8+(block^bs)*4+j)
        M=R[[perm[i]for i in basis]];assert np.max(abs(M@G@M.T-G))<1e-8
        assert np.max(abs(oldgram-oldgram[np.ix_(perm,perm)]))<1e-8
        images=old[90:]@M;family=[];allok=True
        for si,r in enumerate(images):
            key=tuple(np.round(r,6))
            if key in oldkeys:family.append({'existing':names[oldkeys[key]]});continue
            if key in seen:
                ri=seen[key];family.append({'candidate':ri});allok&=not rows[ri]['conflicts'];continue
            pp=cov@r
            valid=(pp<-2-1e-7)|np.logical_or.reduce([abs(pp+2*np.cos(np.pi/m))<1e-6 for m in(2,3,4,5,6,8,10,12)])
            bad=np.flatnonzero(~valid);assert not any(j<16 for j in bad),(perm,si,bad)
            stats['new_individual_images']+=1
            if len(bad)>2:allok=False;family.append({'omitted_conflicts':len(bad)});continue
            ri=len(rows);seen[key]=ri;stats['at_most_two_optional_conflicts']+=1;stats['all102_compatible']+=int(len(bad)==0)
            rows.append({'source':d['names'][90+si],'seed_permutation':perm,'normal':r.tolist(),'conflicts':[names[j]for j in bad],
                         'current_target_cut':float(-r@y),'coordinate_stabilizer_union_target_cut':float(-r@target2),'center':float(r@h),
                         'finite_relations':{names[j]:min((2,3,4,5,6,8,10,12),key=lambda m:abs(pp[j]+2*np.cos(np.pi/m)))for j in range(len(R))if valid[j]and pp[j]>-2+1e-6}})
            family.append({'candidate':ri});allok&=not len(bad)
        stats['group_elements']+=1
        if allok and any('candidate'in z for z in family):families.append({'seed_permutation':perm,'members':family})
    out={'status':'NUMERICAL_FULL_SEED_SYMMETRY_IMAGE_SCOUT','source':source,'source_sha256':V.sha(source),'stats':stats,'coherent_six_image_families':families,'candidates':rows,
         'second_target':target2.tolist(),'scope':'Full2304 seed-diagram action on six certified source normals. Exact action and accepted literal candidates still require replay. Numerical deduplication and pair gates; at most two optional conflicts are retained. No completeness beyond this finite image family.'}
    (V.ROOT/(S.PREFIX+'_full_orbit_scout.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(stats,'COHERENT_FAMILIES',len(families),flush=True)
    good=[(r['coordinate_stabilizer_union_target_cut'],i,r['source'])for i,r in enumerate(rows)if not r['conflicts']]
    print('BEST_REMOTE_CUTTERS',sorted(good,reverse=True)[:20],flush=True)

if __name__=='__main__':run()
