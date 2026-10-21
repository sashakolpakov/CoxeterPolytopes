#!/usr/bin/env python3
"""Revisit rejected algebraic branches with explicit optional-wall conflicts."""
import json,sys,math
from pathlib import Path
from collections import Counter
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import fresh_h9_global_h5_20261013_scout as S
def run():
 R,G,c,names,d,p,C,NR,refs,sources,counts=S.load();data=json.loads((S.P/'data/enumerations'/f'{S.TAG}_recipes.json').read_text());lookup=dict(zip(names,R));original=json.loads((S.P/'data/enumerations/fresh_h9_h5_closure_next_20261013_scout.json').read_text())['candidates']
 for spec in data['extra_anchors']:lookup[spec['name']]=np.array(original[spec['source_row']]['normal'])
 x=np.array(data['target_point']);wanted={(tag,i,sg)for tag,i,sg,cut in data['bad_branches']};out=[];stats=Counter();cache={};cosines=-2*np.cos(np.pi/np.arange(2,257));hard={n for n in names if n.startswith('seed')}|{'compact1','compact2','compact3'}
 for row in data['recipes']:
  for branch in(-1,1):
   if(row['source'],row['row'],branch)not in wanted:continue
   key=tuple(row['incident'])
   if key not in cache:
    E=np.array([lookup[n]for n in key]);MI=np.linalg.inv(E@G@E.T);z=c-(E@G@c)@MI@E;cache[key]=(E,MI,z,z@G@z)
   E,MI,z,qz=cache[key];b=-2*np.cos(np.pi/np.array(row['profile']));part=b@MI@E;rad=(2-part@G@part)/qz
   if rad<=0:continue
   r=part+branch*np.sqrt(rad)*z;ps=R@G@r;bad=[i for i,v in enumerate(ps)if v>-2+1e-6 and min(abs(v-cosines))>1e-6];stats['branches']+=1
   if any(names[i]in hard for i in bad):stats['hard_conflict']+=1;continue
   stats['hard_compatible']+=1
   out.append({'source':row['source'],'row':row['row'],'incident':row['incident'],'profile':row['profile'],'branch':branch,'normal':r.tolist(),'target_cut':float(r@G@x),'conflicts':[{'wall':names[i],'pairing':float(ps[i])}for i in bad],'removed_count':len(bad)})
 out.sort(key=lambda r:(r['removed_count'],-r['target_cut']));doc={'status':'NUMERICAL_OPPOSITE_BRANCH_REPLACEMENT_SCREEN','counts':dict(stats),'hard_retained':sorted(hard),'candidates':out,'scope':'Bounded rejected branches from the saved formula bank. Labels2..256 are a numerical convenience; no exact all-label rejection follows.'}
 (S.P/'data/enumerations'/f'{S.TAG}_opposite_scout.json').write_text(json.dumps(doc,indent=2)+'\n');print(dict(stats));print([(r['removed_count'],r['target_cut'],[x['wall']for x in r['conflicts']])for r in out[:15]],flush=True)
if __name__=='__main__':run()
