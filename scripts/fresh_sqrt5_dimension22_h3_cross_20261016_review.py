#!/usr/bin/env python3
"""Independent graph/Fraction review of the three-hub rank22 weight rule.

No optimizer or producer imports. Graph types and all-root canonical keys
come from the frozen separate review engine; the three-count payments and
reverse-type forest dynamic program are implemented here.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import fresh_sqrt5_enlarged_a1_weights_20261015_review as R
import fresh_sqrt5_parameter_weight_stress_20261015_cross_review as X

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = 'fresh_sqrt5_dimension22_h3_cross_20261016'
TAG = PRIMARY + '_review'
BOUND, CAP, H3CAP = 22, 2, 1
TARGET = F(44)
KINDS = ('connected', 'a1cross', 'hubcross', 'h3cross')
HUBS = ('A1', 'I2(5)', 'H3')


def read(tag):
    return json.loads((ROOT / f'results/{tag}.json').read_text())


def contexts(rank):
    """All exact hub counts, with residual rank zero or at least two."""
    return [(a,b,c) for a in range(rank+1) for b in range((rank-a)//2+1)
            for c in range((rank-a-2*b)//3+1) if rank-a-2*b-3*c != 1]


def corners(row):
    matrix = R.label_matrix(row)
    edges = [(i,j,matrix[i][j]) for i,j in combinations(range(row['nodes']),2)]
    out = []
    for left,right in product((0,2),(1,3)):
        parts = R.components(R.graph(row['white']+[left,right],edges))
        lc = next(g for g in parts if left in g)
        rc = next(g for g in parts if right in g)
        keys = {kind:None for kind in KINDS}
        base = F(0)
        if right in lc:
            keys['connected'] = R.tree_key(lc,[left,right])
            base = R.baseline(lc,left,right)
        else:
            lt,rt = R.type_of(lc),R.type_of(rc)
            for kind,hub,excluded in [('a1cross','A1',('A1',)),
                                       ('hubcross','I2(5)',('A1','I2(5)')),
                                       ('h3cross','H3',HUBS)]:
                if lt == hub and rt not in excluded:
                    keys[kind] = R.tree_key(rc,[right])
                elif rt == hub and lt not in excluded:
                    keys[kind] = R.tree_key(lc,[left])
        assert sum(keys[k] is not None for k in KINDS) <= 1
        forest = [R.type_of(g) for g in parts]
        out.append(dict(baseline=base,**keys,a=forest.count('A1'),
                        b=forest.count('I2(5)'),c=forest.count('H3')))
    return out


class Rule:
    def __init__(self,primary):
        assert primary['status'] == 'PASS_EXACT_DIMENSION22_H3_CROSS_MODEL'
        assert primary['dimension'] == BOUND and F(primary['target']) == TARGET
        assert primary['I25_count_context_cap'] == CAP and primary['I25_context_mode'] == 'capped'
        assert primary['H3_count_context_cap'] == H3CAP and primary['H3_context_mode'] == 'presence'
        assert primary['certificate']['kind'] == 'exact_feasible_rule'
        self.x = list(map(F,primary['certificate']['exact_primal']))
        assert len(self.x) == primary['variables']+1 and min(self.x)>=0
        assert all((v*10000).denominator == 1 for v in self.x[:-1])
        labels = {}
        self.orbits = {kind:{} for kind in KINDS}
        for kind,records in [('connected',primary['connected_orbits']),
                             ('a1cross',primary['cross_orbits']),
                             ('hubcross',primary['new_cross_orbits']),
                             ('h3cross',primary['h3_cross_orbits'])]:
            for rec in records:
                assert len(rec['labels']) == len(rec['nodes'])*(len(rec['nodes'])-1)//2
                graph = R.graph(rec['nodes'],[(i,j,m) for (i,j),m in
                                zip(combinations(rec['nodes'],2),rec['labels'])])
                marks = rec['marked'] if kind == 'connected' else [rec['marked_root']]
                typ = R.type_of(graph)
                assert typ == R.alias(rec['type']) and len(graph) == rec['rank'] <= BOUND
                assert set(marks)<=set(graph) and len(set(marks))==len(marks)
                if kind == 'connected':
                    assert len(marks)==2 and R.baseline(graph,*marks)==0
                else:
                    excluded = {'a1cross':('A1',),'hubcross':('A1','I2(5)'),
                                'h3cross':HUBS}[kind]
                    assert len(marks)==1 and typ not in excluded
                key = R.tree_key(graph,marks)
                assert rec['key'] not in labels and key not in self.orbits[kind]
                labels[rec['key']] = kind,key,len(graph)
                self.orbits[kind][key] = dict(rank=len(graph),label=rec['key'])
        self.index = {}
        assert len(primary['variable_order']) == primary['variables']
        for i,var in enumerate(primary['variable_order']):
            kind,key,rank = labels[var['orbit']]
            a,b,c = var['A1_count'],var['I25_context'],var['H3_context']
            assert isinstance(a,int) and 0<=a<=BOUND-rank
            assert b in range(CAP+1) and c in range(H3CAP+1)
            entry = (kind,key,a,b,c)
            assert entry not in self.index
            self.index[entry] = i
        wanted = set()
        for kind,orbits in self.orbits.items():
            hub_rank = {'connected':0,'a1cross':0,'hubcross':2,'h3cross':3}[kind]
            for key,rec in orbits.items():
                for a in range(int(kind=='a1cross'),BOUND-rec['rank']+1-hub_rank):
                    for b in range(int(kind=='hubcross'),CAP+1):
                        for c in range(int(kind=='h3cross'),H3CAP+1):
                            wanted.add((kind,key,a,b,c))
        assert wanted == set(self.index) and len(wanted)==primary['variables']

    def bonus(self,kind,key,a,b,c):
        idx = self.index.get((kind,key,a,min(b,CAP),min(c,H3CAP)))
        return F(0) if idx is None else self.x[idx]

    def payment(self,row):
        remaining = BOUND+2-row['nodes']
        assert remaining>=0
        cs = corners(row)
        result = []
        for ea,eb,ec in contexts(remaining):
            values = []
            for corner in cs:
                a,b,c = corner['a']+ea,corner['b']+eb,corner['c']+ec
                value = corner['baseline'] + self.bonus('connected',corner['connected'],a,b,c)
                for kind,multiplicity in [('a1cross',a),('hubcross',2*b),('h3cross',3*c)]:
                    key = corner[kind]
                    if key in self.orbits[kind]:
                        assert multiplicity>0
                        value += self.bonus(kind,key,a,b,c)/multiplicity
                values.append(value)
            total = sum(values,F(0))
            assert min(values)>=0 and total>=1,(row['nodes'],ea,eb,ec,values)
            result.append(dict(external_A1=ea,external_I25=eb,external_H3=ec,
                               corner_weights=list(map(str,values)),payment=str(total)))
        return result


def component_bounds(rule):
    models = {name:g for name,g in R.MODELS.items() if len(g)<=BOUND}
    assert len(models)==3*BOUND+5==71
    profiles = {}
    observed = {kind:Counter() for kind in KINDS}
    for name,g in models.items():
        assert R.type_of(g)==name
        base = sum((R.baseline(g,u,v) for u,v in combinations(g,2)),F(0))
        connected = Counter(R.tree_key(g,[u,v]) for u,v in combinations(g,2))
        rooted = Counter(R.tree_key(g,[u]) for u in g)
        profiles[name] = base,connected,rooted
        for key,count in connected.items():
            if key in rule.orbits['connected']:observed['connected'][key]+=count
        for kind in KINDS[1:]:
            for key,count in rooted.items():
                if key in rule.orbits[kind]:observed[kind][key]+=count
    assert all(set(observed[kind])==set(rule.orbits[kind]) for kind in KINDS)
    fronts,totals = [],{}
    for a,b,c in contexts(BOUND):
        remaining = BOUND-a-2*b-3*c
        values = {}
        for name,(base,connected,rooted) in profiles.items():
            value = base + sum((count*rule.bonus('connected',key,a,b,c)
                               for key,count in connected.items()),F(0))
            for kind,multiplicity,excluded in [('a1cross',a,('A1',)),
                                               ('hubcross',2*b,('A1','I2(5)')),
                                               ('h3cross',3*c,HUBS)]:
                if multiplicity and name not in excluded:
                    for key,count in rooted.items():
                        bonus = rule.bonus(kind,key,a,b,c)
                        assert sum((bonus/multiplicity for _ in range(multiplicity)),F(0))==bonus
                        value += count*bonus
            values[name]=value
        assert values['A1']==0
        dp={0:F(0)}
        # Types, not total ranks, form the outer loop. This is an independent
        # unbounded-coin recurrence, with all three exact hub types removed.
        for name in sorted(models,reverse=True):
            if name in HUBS:continue
            r=len(models[name])
            for k in range(r,remaining+1):
                if k-r in dp:
                    dp[k]=max(dp.get(k,F(-1)),dp[k-r]+values[name])
        assert remaining in dp
        maximum = dp[remaining]+a*values['A1']+b*values['I2(5)']+c*values['H3']
        totals[a,b,c]=values
        fronts.append(dict(context=[a,b,c],maximum=str(maximum),
                           residual_rank_maxima=[str(dp[k]) if k in dp else None for k in range(remaining+1)]))
    assert len(fronts)==410
    return models,totals,fronts,observed


def run():
    primary=read(PRIMARY+'_cap2');rule=Rule(primary)
    old=read('fresh_sqrt5_parameter_weight_stress_20261015')
    queue=read('fresh_sqrt5_nonsingular_rank29_fast_20261015')
    assert queue['complete_queue'] and queue['state']['next_core']==queue['marked_states']==len(queue['records'])==5172
    assert not queue['counts'].get('PENDING',0)
    singular=[row for row in queue['records'] if row['result']['singular']]
    assert len(singular)==2070 and [r['index'] for r in singular]==queue['weight_obligations']
    faces=[('prior_catalogue',i,r['row']) for i,r in enumerate(old['faces']) if r['row']['nodes']<=BOUND+2]
    old_count=len(faces)
    faces += [('finite_queue',r['index'],r['row']) for r in singular if r['row']['nodes']<=BOUND+2]
    assert len(faces)==primary['faces']==2151
    expected={kind:set() for kind in KINDS}
    for rec in old['orbits']:
        if rec['rank']<=BOUND:
            g=R.graph(rec['nodes'],[(i,j,m) for (i,j),m in zip(combinations(rec['nodes'],2),rec['labels'])])
            expected['connected'].add(R.tree_key(g,rec['marked']))
    payments=[];case_count=0;minimum=None
    for family,index,row in faces:
        for corner in corners(row):
            for kind in KINDS:
                if corner[kind] is not None and (kind!='connected' or corner['baseline']==0):
                    expected[kind].add(corner[kind])
        cases=rule.payment(row)
        assert len(cases)==len(contexts(BOUND+2-row['nodes']))
        for case in cases:
            slack=F(case['payment'])-1
            minimum=slack if minimum is None else min(minimum,slack)
        case_count+=len(cases)
        payments.append(dict(family=family,index=index,nodes=row['nodes'],cases=cases))
    assert all(expected[kind]==set(rule.orbits[kind]) for kind in KINDS)
    assert case_count==primary['face_cases']==89255
    raw=X.inputs();assert len(raw)==len(old['literal_maps'])==2397
    literals=0
    for (family,index,row,ratio),mp in zip(raw,old['literal_maps']):
        assert (family,index)==(mp['family'],mp['index'])
        target=old['faces'][mp['representative']]
        assert ratio==target['dotted_ratio']
        R.check_bijection(R.label_matrix(row),R.label_matrix(target['row']),mp['to_representative'])
        assert (row['nodes']<=BOUND+2)==(target['row']['nodes']<=BOUND+2)
        literals+=row['nodes']<=BOUND+2
    models,totals,fronts,multiplicities=component_bounds(rule)
    saved=primary['verification']['fronts']
    assert len(saved)==len(fronts)==primary['verification']['contexts']==410
    saved_by_context={tuple(row['context']):row for row in saved}
    assert len(saved_by_context)==len(saved)
    assert set(saved_by_context)==set(contexts(BOUND))
    for ours in fronts:
        context=tuple(ours['context']);theirs=saved_by_context[context]
        assert F(ours['maximum'])==F(theirs['maximum'])
        forest=[R.alias(name) for name in theirs['types']]
        assert sum(len(models[name]) for name in forest)==BOUND
        assert tuple(forest.count(name) for name in HUBS)==context
        assert sum((totals[context][name] for name in forest),F(0))==F(ours['maximum'])
    maximum=max(F(row['maximum']) for row in fronts)
    assert maximum==rule.x[-1]==F(primary['certificate']['maximum'])==F(primary['verification']['maximum'])==F(21977,500)<TARGET
    triangle=read('fresh_sqrt5_triangle_weights_20261015_review')
    assert len(triangle['markings'])==76
    assert all(sum(c['h4_only_weight'] for c in row['corners'])>=2 for row in triangle['markings'])
    paths=[f'scripts/{TAG}.py',f'scripts/{R.TAG}.py',f'scripts/{X.TAG}.py',
           f'results/{R.TAG}.json',f'results/{X.TAG}.json',
           f'scripts/{PRIMARY}.py',f'results/{PRIMARY}_cap2.json',
           'results/fresh_sqrt5_nonsingular_rank29_fast_20261015.json',
           'results/fresh_sqrt5_nonsingular_rank29_fast_20261015_review.json',
           'results/fresh_sqrt5_parameter_weight_stress_20261015.json',
           'results/fresh_sqrt5_triangle_weights_20261015_review.json']
    for path in paths:
        if path.endswith('.json'):
            for dep,h in json.loads((ROOT/path).read_text()).get('sources',{}).items():
                assert sha256((ROOT/dep).read_bytes()).hexdigest()==h,(path,dep)
    return dict(status='PASS_INDEPENDENT_EXACT_THREE_HUB_FACE_AND_FOREST_REVIEW',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        dimension=BOUND,I25_count_context_cap=CAP,H3_count_context_cap=H3CAP,
        original_literal_inputs=len(raw),dimension_feasible_literals=literals,
        old_face_rows=old_count,all_finite_singular_stops=len(singular),
        dimension_feasible_finite_stops=len(faces)-old_count,face_rows=len(faces),
        face_cases=case_count,minimum_face_slack=str(minimum),variables=primary['variables'],
        orbit_counts={kind:len(orbits) for kind,orbits in rule.orbits.items()},
        finite_component_types=len(models),exact_A1_I25_H3_contexts=len(fronts),
        vertex_maximum=str(maximum),strict_face_threshold=str(TARGET),threshold_margin=str(TARGET-maximum),
        forest_maxima=fronts,face_payments=payments,
        multiplicities={kind:[dict(label=rule.orbits[kind][key]['label'],multiplicity=count)
                              for key,count in sorted(counts.items())]
                        for kind,counts in multiplicities.items()},
        proof='Exactly a,2b,3c marked source roots cancel the actual A1,I2(5),H3 denominators. Each cross family excludes all earlier hub types and itself; I2(5) may target H3 but H3 may not target I2(5), so no pair is double counted. Coefficients use exact a,capped b,and H3 presence; denominators always use actual counts. Fixed exact(a,b,c) yields additive components and a reverse-type unbounded Fraction DP excluding exactly the three hub types. Common external triples are enumerated at all four corners. Nonnegative additions preserve the independently reviewed triangular and short-path/H4 baseline.',
        scope='Independent literal graph and exact rational weight verification, with no optimizer or producer imports. Complete geometric branch coverage and field/dimension composition remain separate obligations. No optimality, global bound, external validation or priority claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['vertex_maximum'],result['face_cases'],result['exact_A1_I25_H3_contexts'])
