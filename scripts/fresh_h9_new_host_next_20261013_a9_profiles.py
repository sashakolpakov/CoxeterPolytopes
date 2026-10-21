#!/usr/bin/env python3
"""Complete finite Coxeter coupling profiles on the middle A4+A4 edge."""
from fresh_h9_new_host_next_20261013_a9 import *
from fresh_h9_quartic_mixed_20261013_endpoint import square_root
OUT=ROOT/'results/fresh_h9_new_host_next_20261013_a9_profiles.json'

def run():
    ids=[i for i in range(9)if i!=4];edge=[tuple(-O if j==i+1 else Z for j in range(10))for i in ids]
    eg=tuple(tuple(pair(r,G,s)for s in edge)for r in edge);ei=inverse(eg);rows=[];alphabets=[];nodes=0
    for q in POOL:
        alphabet={Z}
        for ratio in ALLOWED:
            s=square_root(q*ratio/2)
            if s is not None and s.is_integral and(2*s/q).is_integral:alphabet.add(abs(s)if hasattr(s,'__abs__')else s if s>=0 else-s)
        alphabet=sorted(alphabet);alphabets.append({'norm':q.to_json(),'magnitudes':enc(alphabet)})
        def walk(pre):
            nonlocal nodes
            nodes+=1;k=len(pre)
            if k and pair(pre,tuple(row[:k]for row in ei[:k]),pre)>=q:return
            if k==8:
                cf=mv(ei,tuple(-m for m in pre));r0=tuple(sum((c*r[j]for c,r in zip(cf,edge)),Z)for j in range(10));R=q-pair(r0,G,r0)
                if any(R.sign(e)<=0 for e in range(4)):return
                N=10*R;assert N.is_integral and N.norm().denominator==1
                rows.append({'norm':q.to_json(),'magnitudes':enc(pre),'fixed_component':enc(r0),'residual':R.to_json(),'relative_target':N.to_json(),'absolute_norm':int(N.norm())});return
            for v in alphabet:walk((*pre,v))
        walk(())
    out={'status':'COMPLETE_FINITE_COUPLING_PROFILES','edge_indices':ids,'edge_roots':menc(edge),'edge_gram':menc(eg),'norm_pool':enc(POOL),
      'alphabets':alphabets,'profiles':rows,'profile_count':len(rows),'nodes':nodes,
      'absolute_norms':sorted({r['absolute_norm']for r in rows}),
      'relative_targets':enc(sorted({K.from_json(r['relative_target'])for r in rows})),
      'scope':'All integral reflective finite-angle pairings from a primitive root in the safe norm pool to the middle A4+A4 edge. Positive-definite extension requires residual>0 at all embeddings. This is a necessary profile census for a finite first endpoint; roots and complete relative norm orbits are separate.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print(out['status'],len(rows),len(out['relative_targets']),out['absolute_norms'],flush=True)

if __name__=='__main__':run()
