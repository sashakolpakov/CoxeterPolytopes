#!/usr/bin/env python3
"""Exact polar-cap census for the three unused released-wall H9 seeds."""
import json
from scripts.fresh_h8_alternate_caps_20261013 import *

def run():
 for label,deleted in ((3,5),(3,6),(4,4)):
  path=ROOT/f'data/enumerations/fresh_h8_alternate_{label}_{deleted}_20261013.json';data=json.loads(path.read_text());I=matrix(data['seed_gram_inverse']);G=matrix(data['seed_gram']);center=tuple(map(parse,data['strict_center_y']));old=[tuple(map(parse,r['seed_coefficients']))for r in data['caps']];polar=[]
  for i in range(10):
   if sign(I[i][i])<=0:continue
   angle=div(scale(O,4),mul(I[i][i],G[i][i]))
   if sign(sub(angle,scale(O,4)))<0 and angle not in (Z,O,scale(O,2),scale(O,3),add(O,P),add(scale(O,2),P)):continue
   c=tuple(neg(x)for x in I[i]);polar.append({'seed_index':i,'original_wall':data['retained_original_walls'][i],'norm':enc(I[i][i]),'four_cosine_squared':enc(angle),'seed_coefficients':vec(c)});old.append(c)
  out={'base_source':str(path.relative_to(ROOT)),'polar_caps':polar,'volume_gate':strict_escape(I,old,center)}
  target=ROOT/f'data/enumerations/fresh_h8_alternate_{label}_{deleted}_polar_20261013.json';target.write_text(json.dumps(out,indent=2)+'\n');print({'label':label,'deleted':deleted,'polar_caps':len(polar),'gate':out['volume_gate']['status']},flush=True)
if __name__=='__main__':run()
