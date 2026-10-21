#!/usr/bin/env python3
"""Adaptive coupled caps on the alternate seed with its exact pi/10 polar."""
import argparse
import json
from pathlib import Path
import numpy as np
import scripts.fresh_h8_alternate_coupled_20261013 as C
from scripts.fresh_h8_alternate_caps_20261013 import *
BASE=ROOT/'data/enumerations/fresh_h8_alternate_3_6_20261013.json'
POLAR=ROOT/'data/enumerations/fresh_h8_alternate_3_6_polar_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h8_alternate_polar_coupled_20261013.json'

def run(rounds=3,resume=False):
 data=json.loads(BASE.read_text());polar=json.loads(POLAR.read_text());I=matrix(data['seed_gram_inverse']);C.I_G=matrix(data['seed_gram']);old=[tuple(map(parse,r['seed_coefficients']))for r in data['caps']]+[tuple(map(parse,r['seed_coefficients']))for r in polar['polar_caps']];center=tuple(map(parse,data['strict_center_y']));y=tuple(map(parse,polar['volume_gate']['y']));caps=[];history=[]
 if resume:
  previous=json.loads(OUTPUT.read_text());caps=previous['caps'];history=previous['rounds'];y=tuple(map(parse,history[-1]['volume_gate']['y']))
 for iteration in range(len(history),rounds):
  # Select ten old walls closest to the current escape plus the pi/10 polar.
  margins=[real(dot(a,y))/np.linalg.norm([real(x)for x in a])for a in old[:-1]]
  C.ANCHORS=tuple(sorted(set(np.argsort(margins)[:10].tolist()+[len(old)-1])))
  new,count=C.proposals(I,old,center,y,max_support=4,per_anchor=6);caps.extend(new);gate=C.escape(I,old,caps,center)
  row={'round':iteration,'anchor_indices':list(C.ANCHORS),'numeric_admissible_solutions':count,'exact_new_caps':len(new),'total_new_caps':len(caps),'target_escape_y':vec(y),'volume_gate':gate};history.append(row)
  out={'schema':'alternate-H8-pi10-polar-coupled-caps-v1','base_source':str(BASE.relative_to(ROOT)),'polar_source':str(POLAR.relative_to(ROOT)),'fixed_pairing_maximum_support':4,'kept_proposals_per_anchor':6,'caps':caps,'rounds':history,'scope':'Finite adaptive coupled-cap dictionary with exact norm/angle and strict escape checks. Numeric proposal stage is not exhaustive; no global completion obstruction.'};OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print({'round':iteration,'proposals':count,'new_caps':len(new),'total_caps':len(caps),'gate':gate['status']},flush=True)
  if gate['status']!='EXACT_STRICT_SPACELIKE_ESCAPE':break
  y=tuple(map(parse,gate['y']))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--rounds',type=int,default=3);p.add_argument('--resume',action='store_true');a=p.parse_args();run(a.rounds,a.resume)
