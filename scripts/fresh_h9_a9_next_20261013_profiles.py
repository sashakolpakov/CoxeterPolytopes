#!/usr/bin/env python3
"""Complete finite profiles on four next A9 edges, reduced to two norm fields."""
import hashlib,json
from fresh_h9_new_host_next_20261013_a9 import *
from fresh_h9_quartic_mixed_20261013_endpoint import square_root
SOURCE=ROOT/'results/fresh_h9_new_host_next_20261013_a9_endpoint.json'
FIRST=ROOT/'results/fresh_h9_new_host_next_20261013_a9_first_independent.json'
OUT=ROOT/'results/fresh_h9_a9_next_20261013_profiles.json'

def run():
    sd=json.loads(SOURCE.read_text());first=json.loads(FIRST.read_text());assert first['status']=='PASS_COMPLETE_A9_FIRST_ENDPOINT_WITH_CERTIFIED_PARI_INTERFACE'
    assert first['sources'][str(SOURCE.relative_to(ROOT))]==hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    walls=[tuple(map(K.from_json,r))for r in sd['all_walls']];corner=sd['actual_finite_vertex']['indices'];z=K.from_json(sd['z']);x=K.from_json(sd['time'])
    start=tuple(map(K.from_json,sd['endpoint']));V=tuple(z*t for t in start);assert pair(V,G,V)==-20*ALPHA
    assert all(t.is_integral for t in V)and(z/10).is_integral
    alphabets=[]
    for q in POOL:
        values={Z}
        for ratio in ALLOWED:
            s=square_root(q*ratio/2)
            if s is not None and s.is_integral and(2*s/q).is_integral:values.add(s if s>=0 else-s)
        alphabets.append(tuple(sorted(values)))
    cases=[]
    for released in range(4):
        ids=[i for i in corner if i!=released];edge=[walls[i]for i in ids];w=walls[released]
        eg=tuple(tuple(pair(r,G,s)for s in edge)for r in edge);ei=inverse(eg)
        cf=mv(ei,tuple(pair(w,G,r)for r in edge));b=tuple(w[j]-sum((c*r[j]for c,r in zip(cf,edge)),Z)for j in range(10));B=pair(b,G,b)
        kind='end'if released in(0,3)else'inner';D=ALPHA if kind=='end'else 6*ALPHA;multiplier=20 if kind=='end'else 30;den=20 if kind=='end'else 10
        assert B==(K(F(5,4))if kind=='end'else K(F(5,6)))and b[0]==0 and pair(V,G,b)==0
        assert all(((4 if kind=='end'else 6)*v).is_integral for v in b)
        # For any integral root r, Z=−(5/B)*(r,b) and
        # X=−den*(r,V)/(20alpha) are integral, checked covector-wise.
        zcov=tuple(-5*t/B for t in mv(G,b));xcov=tuple(-den*t/(20*ALPHA)for t in mv(G,V))
        assert all(t.is_integral for t in zcov+xcov)
        assert D==(20*ALPHA*25)/(den*den*B)
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
        cases.append(case);print(kind,released,'profiles',len(rows),'targets',len({tuple(r['relative_target'])for r in rows}),'absolute norms',len(case['absolute_norms']),flush=True)
    out={'status':'COMPLETE_FOUR_NEXT_EDGE_FINITE_PROFILES','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(SOURCE,FIRST)},
      'norm_pool':enc(POOL),'alphabets':[enc(a)for a in alphabets],'cases':cases,
      'scope':'All finite Coxeter coupling profiles in the complete safe norm alphabet on the four edge directions modulo initial reversal. Every integral lattice normal induces integral norm-equation coordinates Z,X by the displayed covectors. No orbit search or first-endpoint claim is contained in this profile artifact.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print(out['status'],flush=True)

if __name__=='__main__':run()
