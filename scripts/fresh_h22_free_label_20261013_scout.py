#!/usr/bin/env python3
"""Numerical scan replacing T=2cos(pi/7) in the 25 coherent cap formulas."""
import json
from collections import Counter
import numpy as np
import fresh_h22_incidence_prism_20261013_verify as V

TAG='fresh_h22_free_label_20261013'

def run():
    M,R,h=V.setup();M=np.array(M,float);R=np.array(R,float);h=np.array(h,float)
    rn=np.diag(R@M@R.T);rec=V.read(V.REC)['recipes']
    c=V.read('results/fresh_h22_pi7_coupled_caps_20261013.json')
    h0=np.array(list(map(V.F,c['rational_center_projection'])),float)
    du,dv=np.array([[float(V.F(x))for x in r]for r in c['rational_duals']])
    records=[]
    for m in (3,4,5,6,7,8,10,12,16,24):
        T=2*np.cos(np.pi/m);caps=[];norms=[]
        for r in rec:
            u=np.array(list(map(V.F,r['u'])),float);v=np.array(list(map(V.F,r['v'])),float)
            q=float(V.F(r['norm']));qu=float(V.F(r['u_square']));qv=float(V.F(r['v_square']))
            d=(q-q*q*T*T*qu/4)/qv
            caps.append(-q*T*u/2+r['radical_branch']*np.sqrt(d)*v);norms.append(q)
        cc=caps[17];hp=h@M@cc;p=h0-hp*cc/2
        qp=float(V.F(c['rational_projection_norm']))-hp*hp/2
        Z=p*np.sqrt(2/-qp)
        caps.extend([(209*Z-241*du)/120,(241*Z-418*dv)/np.sqrt(29281)])
        norms.extend([2.,2.]);caps=np.array(caps);norms=np.array(norms)
        caps[:,-1]=abs(caps[:,-1])
        pp=caps@M@R.T;ss=caps@M@caps.T
        bad=[];cnt=Counter()
        def classify(x):
            if x>=1-1e-7:return 'parallel'if abs(x-1)<1e-7 else 'dotted'
            return next((str(k)for k in range(2,101)if abs(x-np.cos(np.pi/k))<1e-7),None)
        for i in range(25):
            for j in range(812):
                z=-pp[i,j]/np.sqrt(norms[i]*rn[j]);kind=classify(z)
                if kind is None:bad.append(['parent',i,j,float(z)])
                else:cnt[kind]+=1
            for j in range(i):
                z=-ss[i,j]/np.sqrt(norms[i]*norms[j]);kind=classify(z)
                if kind is None:bad.append(['cap',i,j,float(z)])
                else:cnt[kind]+=1
        records.append({'label':m,'bad_pair_count':len(bad),'first_bad_pairs':bad[:30],'accepted_counts':dict(cnt)})
        print(m,'bad',len(bad),'first',bad[:3],flush=True)
    out={'status':'NUMERICAL_FREE_LABEL_FORMULA_SCOUT','records':records,
         'sources':{p:V.sha(p)for p in(V.REC,'results/fresh_h22_pi7_coupled_caps_20261013.json')},
         'scope':'Numerical proposal scan only. All25 cap types retain the same formula parameters; no exact finite-angle, volume, or field claim is made.'}
    (V.P/'results'/f'{TAG}_scout.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':run()
