#!/usr/bin/env python3
"""Source-bound coverage index for the cusp-section proof bridge.

This checks literal factor/deletion/target coverage and dependencies.  The
facet-count/fiber theorem is proved in the accompanying report, not inferred
from sampled slices.  The completed rank-three global proof is source-bound.
"""
import argparse,hashlib,json
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_cusp_coverage_20261014'
WALK='results/fresh_h21_source_corner_walk_20261013.json'
FILES={
 'mixed':'results/fresh_h20_retained_link_20261013.json',
 'absorption':'results/fresh_h20_cusp_sections_20261013_absorption.json',
 'all_retained':'results/fresh_h20_single_factor_classification_20261014.json',
 'rank_two':'results/fresh_h20_single_factor_links_20261013.json',
 'rank_four_local':'results/fresh_h20_omitted_factor_classification_20261014.json',
 'rank_four_global':'results/fresh_h20_omitted_factor_classification_20261014_global.json',
 'rank_three':'results/fresh_h20_rank3_omitted_sections_20261014.json',
 'rank_three_review':'results/fresh_h20_rank3_omitted_sections_20261014_review.json',
 'rank_three_global':'results/fresh_h20_rank3_global_sections_20261014.json',
 'rank_three_global_review':'results/fresh_h20_rank3_global_sections_20261014_review.json',
}
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def key(c,w):return int(c),tuple(sorted(w))
def run():
 walk=json.loads((ROOT/WALK).read_text())
 assert walk['status']=='PASS_COMPLETE_CORNER_ORBIT_CENSUS' and not walk['pending_orbit_ids']
 assert walk['processed_orbits']==len(walk['corner_orbits'])==395
 cusps=[c for c in walk['corner_orbits']if c['type']=='ideal'];assert len(cusps)==18
 factors={key(c['id'],f['walls']):f for c in cusps for f in c['affine_components']}
 assert len(factors)==53 and all(sum(f['rank']for f in c['affine_components'])==20 for c in cusps)
 for f in factors.values():assert len(f['walls'])==f['rank']+1
 docs={name:json.loads((ROOT/p).read_text())for name,p in FILES.items()}
 for doc in docs.values():
  assert doc['status'].startswith('PASS')
  for p,h in doc['sources'].items():assert sha(p)==h
 mixed={key(c['id'],f['walls'])for c in docs['mixed']['cusp_factors']for f in c['factors']}
 assert mixed==set(factors)
 retained={key(c['cusp'],c['walls'])for c in docs['all_retained']['factor_census']}
 assert retained==set(factors)
 rank1={k for k,f in factors.items()if f['rank']==1}
 assert rank1=={(0,(162,164))}
 absorb=docs['absorption'];assert absorb['cusp_id']==0 and absorb['only_rank_one_component_in_18_representatives']
 rank2={k for k,f in factors.items()if f['rank']==2}
 tw=docs['rank_two'];twkeys={key(r['cusp'],r['walls'])for r in tw['full_C2_spans']}
 twkeys|={key(r['cusp'],r['factor'])for r in tw['full_A2_spans']}
 assert twkeys==rank2 and len(rank2)==4
 # Every affine-A2 factor has all three choices of negative coefficient.
 for k in rank2:
  aa=[r for r in tw['full_A2_spans']if key(r['cusp'],r['factor'])==k]
  if aa:assert {r['negative_wall']for r in aa}==set(k[1])
 rank3={k for k,f in factors.items()if f['rank']==3}
 three=docs['rank_three'];assert {key(r['cusp'],r['walls'])for r in three['factors']}==rank3
 assert len(rank3)==3 and three['system_count']==120 and three['direction_count']==91
 three_global=docs['rank_three_global'];assert three_global['pending_endpoint_rows']==[]
 assert sorted(r['row']for r in three_global['families'])==list(range(91))
 endpoint_rows={r['row']for r in three_global['separate_endpoint_certificates']}
 assert len(endpoint_rows)==22
 assert endpoint_rows=={r['row']for r in three_global['families']if not r['u1_is_old_facet_hyperplane']and not r['u1_bad_angle_verified']}
 assert docs['rank_three_global_review']['counts']['families']==91
 for k in rank3:
  for omit in k[1]:
   rows=[r for r in three['systems']if key(r['cusp'],r['factor'])==k and r['omitted_wall']==omit]
   assert len(rows)==10 and Counter(r['target']for r in rows)=={'A2':1,'C2':3,'G2':6}
 rank4={k for k,f in factors.items()if f['rank']>=4}
 lo=docs['rank_four_local'];deletions={(*key(r['cusp'],r['factor']),r['omitted_wall'])for r in lo['deletions']}
 expected={(cid,ws,i)for cid,ws in rank4 for i in ws}
 assert deletions==expected and len(deletions)==387 and len(rank4)==45
 glob=docs['rank_four_global'];covered={r['row']for r in glob['old_difference_plane_matches']}|{r['row']for r in glob['new_plane_global_certificates']}
 assert covered==set(range(len(lo['local_directions']))) and len(covered)==175
 assert rank1|rank2|rank3|rank4==set(factors)
 # No mixed active set here can give the exceptional affine triangle:
 # there is only one rank1 factor, and its partners have ranks5,7,7.
 for c in cusps:
  rr=[f['rank']for f in c['affine_components']]
  if 1 in rr:assert c['id']==0 and sorted(rr)==[1,5,7,7]
 paths=[WALK,*FILES.values(),str(Path(__file__).resolve().relative_to(ROOT))]
 return {'status':'PASS_COMPLETE_CUSP_CONTAINING_INTERIOR_SECTION_COVERAGE_INDEX',
  'sources':{p:sha(p)for p in paths},'cusp_orbits':18,'affine_factors':53,
  'rank_partition':{'rank1':len(rank1),'rank2':len(rank2),'rank3':len(rank3),'rank_at_least4':len(rank4)},
  'one_omitted_rank_at_least4_deletions':len(deletions),'rank3_target_systems':120,'rank3_local_records':91,
  'rank3_factors':[{'cusp':c,'walls':list(w)}for c,w in sorted(rank3)],
  'conclusion':'No Coxeter H20 hyperplane section meeting the source strict interior can contain a source cusp. All91 rank3 records and all interior endpoints are now covered by the source-bound global proof and independent review.',
  'remaining_section_scope':'Every remaining interior Coxeter H20 section would have to be compact and avoid all source cusps. Boundary sections in source-facet supporting hyperplanes are outside this particular theorem.',
  'scope':'Coverage index and source/proof bridge only. Underlying arithmetic exclusions and complete cusp-orbit walk are separately archived dependencies; no duplicate replay of those full computations.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run()
 p=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print(out['rank_partition'])
