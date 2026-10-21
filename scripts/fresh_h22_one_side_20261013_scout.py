#!/usr/bin/env python3
"""Numerical sign-unification gate for the 25 cap types of the 839 branch."""
import json
from collections import Counter
from pathlib import Path
import numpy as np
from fresh_h22_pi7_edgewalk_20261013_scout import load
import fresh_h22_incidence_prism_20261013_verify as V

P = Path(__file__).resolve().parents[1]
TAG = 'fresh_h22_one_side_20261013'


def run():
    M, W, h, G, norms, _, cert = load()
    c = V.read('results/fresh_h22_pi7_coupled_caps_20261013.json')
    extra = np.array([[(int(a)+int(b))/(2*V.S) for a,b in r['normal_intervals']]
                      for r in c['new_caps']])
    caps = np.vstack([W[812:], extra])
    q = np.r_[norms[812:], [2.,2.]]
    signs = np.sign(caps[:,-1]).astype(int)
    caps[:,-1] = abs(caps[:,-1])
    H = caps @ M @ caps.T
    bad, counts = [], Counter()
    for i in range(len(caps)):
        for j in range(i):
            z = -H[i,j]/np.sqrt(q[i]*q[j])
            if z >= 1-1e-8:
                label = 'parallel' if abs(z-1)<1e-8 else 'dotted'
            else:
                label = next((m for m in range(2,101) if abs(z-np.cos(np.pi/m))<1e-8), None)
            if label is None:
                bad.append({'i':i,'j':j,'negative_normalized_pair':float(z),
                            'original_signs':[int(signs[i]),int(signs[j])]})
            else:
                counts[str(label)] += 1
    out = {'status':'NUMERICAL_ONE_SIDE_COMPATIBILITY_SCOUT',
           'sources':{s:V.sha(s) for s in (V.OUT,V.REC,
               'results/fresh_h22_pi7_coupled_caps_20261013.json')},
           'cap_types':len(caps),'original_signs':signs.tolist(),
           'accepted_pair_counts':dict(counts),'bad_pairs':bad,
           'scope':'Floating-point proposal gate only; no exact angle or compatibility claim.'}
    (P/'results'/f'{TAG}_scout.json').write_text(json.dumps(out,indent=2)+'\n')
    print('types',len(caps),'bad',len(bad),'counts',dict(counts))
    print(bad)


if __name__ == '__main__':
    run()
