#!/usr/bin/env python3
"""Complete finite profiles for the outer edge leaving the central H8 facet."""
import hashlib,json
from fresh_h9_new_host_next_20261013_a9 import *
from fresh_h9_quartic_mixed_20261013_endpoint import square_root
SOURCE=ROOT/'results/fresh_h9_a9_next_20261013_geometry.json'
FIRST=ROOT/'results/fresh_h9_a9_next_20261013_independent.json'
PREVIOUS=ROOT/'results/fresh_h9_a9_next_20261013_profiles.json'
OUT=ROOT/'results/fresh_h9_a9_outer_20261013_axis_profiles.json'

def run():
    sd=json.loads(SOURCE.read_text());first=json.loads(FIRST.read_text());assert first['status']=='PASS_COMPLETE_FOUR_A9_FIRST_ENDPOINTS'
    corner_record=sd['actual_finite_vertices'][0];assert corner_record['case']==0 and not corner_record['reversed']
    walls=[tuple(map(K.from_json,r))for r in sd['all_walls']];corner=corner_record['certificate']['indices']
    start=tuple(map(K.from_json,corner_record['point']));assert first['cases'][0]['root']==sd['all_walls'][corner_record['wall_index']]
    scale=square_root(-20*ALPHA/pair(start,G,start));assert scale is not None and scale>0
    V=tuple(scale*t for t in start);assert pair(V,G,V)==-20*ALPHA and all(t.is_integral for t in V)
    alphabets=[]
    for q in POOL:
        values={Z}
        for ratio in ALLOWED:
            s=square_root(q*ratio/2)
            if s is not None and s.is_integral and(2*s/q).is_integral:values.add(s if s>=0 else-s)
        alphabets.append(tuple(sorted(values)))
    cases=[]
    for released in corner:
        if released!=10:continue
        ids=[i for i in corner if i!=released];edge=[walls[i]for i in ids];w=walls[released]
        eg=tuple(tuple(pair(r,G,s)for s in edge)for r in edge);ei=inverse(eg)
        cf=mv(ei,tuple(pair(w,G,r)for r in edge));b=tuple(w[j]-sum((c*r[j]for c,r in zip(cf,edge)),Z)for j in range(10));B=pair(b,G,b)
        assert B==2
        kind='axis';D=10*ALPHA;multiplier=50;den=10
        assert pair(V,G,b)==0 and all(v.is_integral for v in b)
        zcov=tuple(-10*t/B for t in mv(G,b));xcov=tuple(-den*t/(20*ALPHA)for t in mv(G,V))
        assert all(t.is_integral for t in zcov+xcov)and D==(20*ALPHA*100)/(den*den*B)
        rows=[];nodes=0
        for q,alphabet in zip(POOL,alphabets,strict=True):
            def walk(pre):
                nonlocal nodes
                nodes+=1;k=len(pre)
                if k and pair(pre,tuple(row[:k]for row in ei[:k]),pre)>=q:return
                if k==8:
                    co=mv(ei,tuple(-m for m in pre));r0=tuple(sum((c*r[j]for c,r in zip(co,edge)),Z)for j in range(10));R=q-pair(r0,G,r0)
                    if any(R.sign(e)<=0 for e in range(4)):return
                    N=multiplier*R;assert N.is_integral and N.norm().denominator==1
                    assert all(pair(r0,G,v)==0 for v in(b,V))
                    rows.append({'norm':q.to_json(),'magnitudes':enc(pre),'fixed_component':enc(r0),'residual':R.to_json(),'relative_target':N.to_json(),'absolute_norm':int(N.norm())});return
                for value in alphabet:walk((*pre,value))
            walk(())
        case={'released_wall':released,'kind':kind,'edge_indices':ids,'edge_roots':menc(edge),'edge_gram':menc(eg),
          'start':enc(start),'timelike_vector':enc(V),'released_residual':enc(b),'residual_square':B.to_json(),
          'relative_discriminant':D.to_json(),'target_multiplier':multiplier,'timelike_denominator':den,
          'integral_Z_covector':enc(zcov),'integral_X_covector':enc(xcov),'profiles':rows,'nodes':nodes,
          'absolute_norms':sorted({r['absolute_norm']for r in rows})}
        cases.append(case);print(kind,released,'profiles',len(rows),'targets',len({tuple(r['relative_target'])for r in rows}),flush=True)
    assert len(cases)==1 and len(cases[0]['profiles'])==50
    out={'status':'COMPLETE_OUTER_AXIS_EDGE_FINITE_PROFILES','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,FIRST,PREVIOUS)},
      'source_corner':corner,'source_corner_record':0,'source_corner_complete_mirrors':first['cases'][0]['full_endpoint_mirrors'],
      'norm_pool':enc(POOL),'alphabets':[enc(a)for a in alphabets],'cases':cases,
      'scope':'All finite Coxeter profiles in the complete safe norm alphabet on the edge releasing facet10 at the first outer vertex. Integral covectors Z=-5(r,b) and X=-(r,V)/(2alpha) cover every lattice root. The equation is Z²-10alpha X²=50R and the normal is r0-Zb/10+XV/10.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print(out['status'],flush=True)

if __name__=='__main__':run()
