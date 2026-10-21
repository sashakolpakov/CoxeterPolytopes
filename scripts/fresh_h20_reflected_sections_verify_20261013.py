#!/usr/bin/env python3
"""Independent Fraction/Boolean replay of finite spherical-star obstructions."""
from copy import deepcopy
from fractions import Fraction as F
import json
import fresh_h20_joint_cusp_irrational_verify_20261013 as prior
K=prior.K;A=prior.A;B=prior.B;ROOT=prior.ROOT
SOURCE=ROOT/'data/enumerations/fresh_h20_reflected_sections_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h20_reflected_sections_verify_20261013.json'

def boolean_tree(clauses,tree,assignment=None):
 assignment={} if assignment is None else dict(assignment)
 def residual(index):
  assert 0<=index<len(clauses);c=clauses[index]
  assert not any(abs(l) in assignment and assignment[abs(l)]==(l>0) for l in c)
  return [l for l in c if abs(l) not in assignment]
 for literal,reason in tree['units']:
  assert residual(reason)==[literal];assignment[abs(literal)]=literal>0
 if 'conflict_clause' in tree:
  assert residual(tree['conflict_clause'])==[];return 1
 v=tree['branch_variable'];assert 1<=v<=210 and v not in assignment
 return 1+boolean_tree(clauses,tree['false'],{**assignment,v:False})+boolean_tree(clauses,tree['true'],{**assignment,v:True})

def check_case(row,source,Q,W,G,h):
 ci,ri=row['source_case'],row['source_candidate'];original=source['cases'][ci]['candidates'][ri]
 assert row['normal']==original['normal'];field=K.Field(A.trim((-2,0,1)),*row['field_interval']);assert field.sign(A.T)>0
 e=[A.readp(p) for p in row['normal']];sd=prior.products(e,Q,W,G,field);assert field.sign(sd[1])>0
 clauses=row['clauses'];proofs=row['clause_certificates'];assert len(clauses)==len(proofs)
 counts={'ridges':0,'nonelliptic_pairs':0,'nonspherical_subsets':0}
 for clause,proof in zip(clauses,proofs):
  ids=[x-1 for x in proof['walls']];assert len(set(ids))==len(ids) and all(0<=i<210 for i in ids)
  if proof['kind']=='actual_bad_ridge':
   assert clause==proof['walls'] and len(ids)==2
   cert=proof['certificate'];assert cert['ridge_walls_one_based']==proof['walls'] and cert['zero_walls_one_based']==sorted(proof['walls'])
   aa,ab,bb=B.actual(cert,Q,W,G,h,sd,field);N=field.rem(A.mul(ab,ab));den=field.rem(A.mul(aa,bb))
   if cert['angle_rejection']=='obtuse':assert field.sign(ab)>0
   else:
    assert field.sign(ab)<0 and field.sign(N)>0 and field.sign(A.sub(A.scale(den,3),A.scale(N,4)))>0
    for p in (A.sub(A.scale(N,4),den),A.sub(A.scale(N,2),den),A.add(A.sub(A.scale(A.mul(N,N),16),A.scale(A.mul(N,den),12)),A.mul(den,den))):assert field.sign(p)!=0
   counts['ridges']+=1
  elif proof['kind']=='nonelliptic_pair':
   assert clause==[-i-1 for i in ids] and len(ids)==2;a,b=ids
   assert G[a][a]*G[b][b]-G[a][b]**2<=0;counts['nonelliptic_pairs']+=1
  else:
   assert proof['kind']=='nonspherical_subset' and clause==[-i-1 for i in ids]
   v=list(map(F,proof['nonpositive_vector']));assert len(v)==len(ids) and any(v)
   square=sum(v[a]*v[b]*G[i][j] for a,i in enumerate(ids) for b,j in enumerate(ids))
   assert square==F(proof['square']) and square<=0;counts['nonspherical_subsets']+=1
 nodes=boolean_tree(clauses,row['boolean_unsat_tree']);assert row['status']=='NO_FINITE_SPHERICAL_REFLECTION_STAR'
 return {**counts,'boolean_tree_nodes':nodes}

def main():
 data=json.loads(SOURCE.read_text());assert data['complete'];source=json.loads((ROOT/data['source']).read_text());Q,W,G,h=A.source()
 expected={(i,j) for i,c in enumerate(source['cases']) for j in range(len(c['candidates']))}
 assert {(c['source_case'],c['source_candidate']) for c in data['cases']}==expected and len(data['cases'])==len(expected)==12
 results=[check_case(c,source,Q,W,G,h) for c in data['cases']]
 forged=deepcopy(data['cases'][0]);cert=next(p['certificate'] for p in forged['clause_certificates'] if p['kind']=='actual_bad_ridge');cert['point_coordinate_polynomials'][0][-1]=str(F(cert['point_coordinate_polynomials'][0][-1])+1)
 try:check_case(forged,source,Q,W,G,h)
 except AssertionError:pass
 else:raise AssertionError('forged actual ridge accepted')
 forged=deepcopy(data['cases'][0]);node=forged['boolean_unsat_tree']
 while 'conflict_clause' not in node:node=node['false']
 assert node['units'];node['units'][0][0]*=-1
 try:check_case(forged,source,Q,W,G,h)
 except AssertionError:pass
 else:raise AssertionError('forged Boolean contradiction accepted')
 result={'status':'PASS','standard_library_only':True,'fixed_section_planes':12,'case_counts':results,'forged_ridge_rejected':True,'forged_boolean_proof_rejected':True,
  'scope':data['scope'],'imports_no_sat_solver':True}
 OUTPUT.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
