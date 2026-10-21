#!/usr/bin/env python3
"""Logical field coverage binder; the compact dimension21 case stays open."""
import argparse
import json
from functools import lru_cache
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_compact21_field_reduction_20261016'
TAGS=['fresh_arithmetic_high_label_20261016','fresh_arithmetic_high_label_20261016_review',
      'fresh_sqrt2_quad_exhaustion_20261015_review','fresh_sqrt2_weight_bound_20261015_review',
      'fresh_compact_global_bound_20261015_review',
      'fresh_arithmetic_dimension21_composition_20261016_bound21']


@lru_cache(None)
def sha(path): return sha256((ROOT/path).read_bytes()).hexdigest()


def run():
    documents={tag:json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in TAGS}
    assert all(d['status'].startswith('PASS') for d in documents.values())
    high,review,q2coverage,q2weights,k7,previous=[documents[t] for t in TAGS]
    assert high['status']=='PASS_EXACT_DISTINGUISHED_BOND_REDUCTION_AND_EIGHT_FIELD_DIMENSION20_CANDIDATE'
    assert review['status']=='PASS_INDEPENDENT_EIGHT_FIELD_HIGH_BOND_REVIEW'
    assert q2coverage['status']=='PASS_INDEPENDENT_UNION_FIELD_PROFILE_REPLAY_AND_ALL_SIZE_CLOSURE_PROOF_REVIEW'
    assert k7['status']=='PASS_INDEPENDENT_FRACTION_AND_EXHAUSTIVE_QUAD_SCOPE_REVIEW_CUBIC_SEVEN_BOUND20'
    q2at21=next(row for row in q2weights['hvector_bounds'] if row['dimension']==21)
    assert q2at21['strict_weight_lower']=='40' and q2at21['upper']==40 and q2at21['gap']=='0'
    assert previous['compact_arithmetic_dimension_upper_bound']==21
    assert previous['finite_volume_arithmetic_reflection_dimension_upper_bound']==21
    fields=['Q(sqrt2)','Q(sqrt3)','Q(sqrt5)','Q(sqrt6)','Q(sqrt2,sqrt3)',
            'Q(sqrt2,sqrt5)','K7','K9','K11','K15','K16','K20']
    assert fields==high['fields']['vinberg_compact_dimension_at_least14_fields']
    new=high['fields']['new_fields']
    canonical={'Qsqrt3':'Q(sqrt3)','Qsqrt6':'Q(sqrt6)','Qsqrt2sqrt3':'Q(sqrt2,sqrt3)'}
    assert new==[canonical.get(row['field'],row['field']) for row in review['fields']]
    assert len(new)==8 and review['counting']['strict_threshold']==review['counting']['vertex_maximum']==40
    branches=[dict(field='Q',dimension=21,status='EXCLUDED_PUBLISHED',
                   reason='Meyer isotropy in rank>=5 rules out compact rational simplest-type arithmetic lattices in dimension>=4.')]
    for field in fields:
        if field in new:
            row=dict(status='EXCLUDED_RESEARCH_CANDIDATE',upper_bound=20,
                     evidence=TAGS[:2])
        elif field=='Q(sqrt2)':
            row=dict(status='EXCLUDED_RESEARCH_CANDIDATE',upper_bound=20,evidence=TAGS[2:4])
        elif field=='K7':
            row=dict(status='EXCLUDED_RESEARCH_CANDIDATE',upper_bound=20,evidence=[TAGS[4]])
        else:
            row=dict(status='OPEN',reason='Not covered by the new field hypothesis or existing bound20 certificates.')
        branches.append(dict(field=field,dimension=21,**row))
    remaining=[row['field'] for row in branches if row['status']=='OPEN']
    assert remaining==['Q(sqrt5)','Q(sqrt2,sqrt5)']
    checked=set();references=0
    def check(path):
        nonlocal references
        if path in checked or not path.endswith('.json'): return
        checked.add(path)
        data=json.loads((ROOT/path).read_text())
        for field in ('sources','source_hashes'):
            for source,expected in data.get(field,{}).items():
                assert isinstance(expected,str) and sha(source)==expected,(path,source)
                references+=1;check(source)
    paths=[f'results/{tag}.json' for tag in TAGS]+[f'scripts/{TAG}.py']
    for path in paths: check(path)
    return dict(status='PASS_COMPACT_DIMENSION21_REDUCTION_TO_TWO_OPEN_FIELDS',
                sources={path:sha(path) for path in paths},branches=branches,
                remaining_fields=remaining,dimension21_excluded=False,
                compact_arithmetic_upper_bound_unchanged=21,
                new_field_bound20_candidates=8,recursive_manifest_count=len(checked),source_reference_count=references,
                published_inputs=[dict(result='Dimension>=14 arithmetic field list, including Q; compact rational case removed separately.',
                                       source='Vinberg1984; Nikulin0708.3991 introduction and Section3',url='https://arxiv.org/pdf/0708.3991'),
                                  dict(result='Arithmetic reflection lattices are of simplest type; rational isotropy prevents cocompactness in dimension>=4.',
                                       source='Vinberg simplest-type criterion; Meyer theorem; arithmetic compactness criterion',url='https://arxiv.org/pdf/1506.03111')],
                scope='Dependency and logical-coverage composition of internally reviewed research proof candidates with published inputs. Not a new compact upper bound20, not independent reruns of inherited large enumerations, not external validation or priority.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();result=run();path=ROOT/'results'/f'{TAG}.json'
    if args.verify: assert result==json.loads(path.read_text())
    else: path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['remaining_fields'],result['recursive_manifest_count'])
