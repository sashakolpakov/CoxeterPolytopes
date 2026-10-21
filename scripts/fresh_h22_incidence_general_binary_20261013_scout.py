#!/usr/bin/env python3
"""Search binary point profiles with independently varying norm-two height.

For A integer and u a 21-bit vector, x=(A*vprime-sum(u_i P_i))/2
has x^2=(|u|-A^2)/2.  The lift (x,sqrt((A^2-|u|+4)/2))
has norm two.  The old210 gate below is exact integer arithmetic.  The
escape cut used to filter proposals is numerical; this is not a certificate.
"""
import json
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parents[1]
PRE='fresh_h22_incidence_general_binary_20261013'
def run():
 raw=json.loads((P/'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json').read_text())['coxeter_input']
 Q=-np.array(raw['basis_gram'],dtype=np.int64);W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64)
 core=json.loads((P/'results/fresh_h22_incidence_core_20261013.json').read_text());points,lines=core['incidence_parts'];vp=np.array(core['old_projections'][1],dtype=np.int64)
 inc=-(W[points]@Q@W[lines].T);wn=np.einsum('ij,ij->i',W,W@Q);long=np.flatnonzero(wn==4)
 escape=json.loads((P/'results/fresh_h22_transported_cusp_completion_20261013_escape.json').read_text());e=np.array(escape['numerator'],dtype=float)/escape['denominator'];ep=W[points]@Q@e[:22];ev=float(vp@Q@e[:22])
 records=[];counts={str(a):{'line_gate':0,'numerical_cut':0,'old210_gate':0}for a in range(2,13)}
 for start in range(0,1<<21,8192):
  masks=np.arange(start,min(start+8192,1<<21),dtype=np.int64);u=((masks[:,None]>>np.arange(21))&1).astype(np.int64);weights=u.sum(axis=1);c=u@inc;cuts=u@ep
  for a in range(2,13):
   pair=c-a;ok=np.all((pair==0)|(pair==-2)|(pair<=-4),axis=1);z2=(a*a-weights+4)/2;ok&=z2>0;counts[str(a)]['line_gate']+=int(ok.sum())
   val=a*ev/2-cuts/2+np.sqrt(np.maximum(z2,0))*abs(e[22]);ok&=val>1e-8;ix=np.flatnonzero(ok);counts[str(a)]['numerical_cut']+=len(ix)
   if not len(ix):continue
   rows=a*vp-u[ix]@W[points];pairs=rows@Q@W[long].T;good=np.all((pairs==0)|(pairs==-4)|(pairs<=-6),axis=1);counts[str(a)]['old210_gate']+=int(good.sum())
   for j,row in zip(ix[good],rows[good]):records.append({'A':a,'mask':int(masks[j]),'weight':int(weights[j]),'twice_old_projection':row.tolist(),'twice_height_square':int(a*a-weights[j]+4),'numerical_escape_pair':float(val[j])})
  if start%(1<<18)==0:print('profiles',start,'accepted',len(records),flush=True)
 out={'status':'NUMERICAL_ESCAPE_FILTER_EXACT_OLD210_BINARY_PROFILE_SCOUT','A_range':[2,12],'all_binary_profiles':1<<21,'counts':counts,'candidates':records,'scope':'Old210 Coxeter pairs exact; escape filtering numerical. Other retained banks and mutual pairs are not checked.'}
 (P/('data/enumerations/'+PRE+'_scout.json')).write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(records),counts,flush=True)
if __name__=='__main__':run()
