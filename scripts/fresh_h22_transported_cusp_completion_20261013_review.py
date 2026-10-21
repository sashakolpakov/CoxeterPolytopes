#!/usr/bin/env python3
"""Replay the separate local H22 review and structural facet certificates.

The full native orbit/pair and escape maxima replays remain separately named
primary commands. This wrapper uses no numerical arithmetic.
"""
import argparse,json
from fresh_h22_transported_cusp_completion_20261013_local_independent import (
 run as local_run,P,F,read,sha,dot,mv,RAW,SOURCE)
from fresh_h22_transported_cusp_completion_20261013_facet_obstruction import run as facet_run
PRE='fresh_h22_transported_cusp_completion_20261013';OUT='results/'+PRE+'_review.json'
def run():
 local=local_run();lp='results/'+PRE+'_local_independent.json';assert local==read(lp)
 facets=facet_run();fp='results/'+PRE+'_facet_obstruction.json';assert facets==read(fp)
 ep='results/'+PRE+'_escape.json';escape=read(ep);assert all(sha(p)==h for p,h in escape['sources'].items())
 raw=read(RAW)['coxeter_input'];Q=[[-F(x)for x in row]for row in raw['basis_gram']];v=list(map(F,escape['numerator']));den=escape['denominator'];cv=mv(Q,v[:22]);q=(dot(v[:22],cv)+v[22]*v[22])/(den*den);assert q>0 and q==F(escape['norm'])
 count=0
 for row in escape['new_family_integer_maxima']:
  a=int(row['maximum_old_integer_pairing']);b=int(row['absolute_sqrt11_numerator']);assert a<0 and b>=0 and a*a-11*b*b==int(row['positive_square_margin'])>0;count+=row['covered_signed_walls']
 assert count==6289920 and escape['facets']==11688030
 sources=(SOURCE,lp,fp,ep,RAW)
 return {'status':'PASS_INDEPENDENT_LOCAL_AND_STRUCTURAL_H22_REVIEW','sources':{p:sha(p)for p in sources},'local_independent_status':local['status'],'structural_facet_status':facets['status'],'ambient_signature':[22,1,0],'local_normals':31,'positive_quotient_rank':21,'positive_fullsupport_delta_multiple':66,'repaired_specified_signed_cusps':483840,'facet_count':11688030,'exact_strict_escape_square':str(q),'scope':'Separate stdlibFraction reconstruction replays representative13norms,2730oldpairs,465localpairs,ranks andfullsupportrelation; separatefacetdirectionsandtheiropennessarereplayed. Theescapesquareandsavedinteger-rootmarginsarerecomputedhere. The factthesavedmaximaareglobal,andfullorbitcompleteness/distantpairgates,dependonthe separatelyreplayed primaryverifiers;theyarenotindependentlyenumerated bythiswrapper.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(P/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
