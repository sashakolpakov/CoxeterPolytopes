#!/usr/bin/env python3
"""SOURCE ONLY: exact selected D13 negative-guard completion bridge.

Only ambient58 and this guarded selected domain are checked. Numerical
projection/radical certificates and their source-proved geometric implications
remain explicit inputs; this is not a full-terminal or global finite57 PASS.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import types

ROOT=Path(__file__).resolve().parents[1]
ROUTES='scripts/index_finite57_current_terminal_routes_coverage_safe_20261025.py'
ROUTES_SHA='65c7a7e018800fa7bd259a984ec24249656410569aee578f668a6e04f5fde5d1'
SHARDS=8
PAIRS={(0,11),(0,12)}
MAX_OUTPUT=500*1024


def need(ok,message):
    if not ok:raise ValueError(message)


def raw(relative,cap=2*1024**2):
    path=(ROOT/relative).resolve();need(path.is_relative_to(ROOT)and path.is_file(),'bounded source path')
    with path.open('rb')as stream:data=stream.read(cap+1)
    need(0<len(data)<=cap,'bounded source read')
    return data


def load():
    data=raw(ROUTES);need(sha256(data).hexdigest()==ROUTES_SHA,'route-index source pin')
    obj=types.ModuleType('finite57_selected_guard_routes');obj.__file__=str(ROOT/ROUTES)
    exec(compile(data,str(ROOT/ROUTES),'exec'),obj.__dict__)
    return obj,obj.prepare()


def canon(c,values):return sorted(values,key=c.canonical)


def domain(c,r):
    cases,cores,partition=r.domain(c)
    eligible=[p for p in cases if(p.get('component_order_bound')or 0)>=13 and len(p['marked_component'])<=13]
    doc=c.document(c.FINAL,('selected_support',),('selected_proofs',),('ambient_dimension',))
    need(doc['ambient_dimension']==58,'ambient58 final guard premise')
    support=doc['selected_support'];saved=support['rows']
    need(len(eligible)==len(saved)==801,'exact eligible801 guard domain')
    need([p['index']for p in eligible]==[p['index']for p in saved],'full guard801 high-index order')
    expected={p['index']:p for p in saved}
    proofs={p['index']:p for p in doc['selected_proofs']}
    need(len(proofs)==len(doc['selected_proofs'])==33,'unique33 selected guarded proofs')
    affected={p['index']for p in saved if p['affected_maps']}
    need(set(proofs)==affected and sum(len(p['affected_maps'])for p in saved)==76,'exact33/76 guarded provider domain')
    need(support['patterns']==801 and support['maps']==862 and support['affected_contexts']==33 and support['affected_maps']==76,
         'saved guard aggregate metadata agrees with exact key domain')
    return eligible,cores,expected,proofs,partition


def proof(c,p,high,core,maps):
    need(p['core_key']==high['core_key']and p['gram']==core['gram'],'selected guarded proof original physical core')
    need(p['selected_nodes']==high['marked_component']and p['selected_marks']==high['preferred_corner'],
         'selected guarded proof base/mark association')
    g=c.E.matrix(core['gram']);base=high['marked_component'];marks=high['preferred_corner']
    left=high['left'][1];right=[high['right_copy_map'][i]for i in high['right'][1]]
    white=sorted(set(range(len(g)))-set(left+right))
    need(p['white']==white,'selected guard original white set')
    method=p['method'];result=dict(high_index=high['index'],core_key=high['core_key'],method=method,
                                  full_map_keys=maps,guard='three distinct whole A15/D15 spectators at this actual finite corner')
    if method=='genuine_new_fork_leaf_projection_exclusion':
        need(len(base)==12,'proper selected D13 proof base has order12')
        saved=p['maps']
        need(canon(c,[q['images']for q in saved])==canon(c,[q['images']for q in maps]),
             'every actual proper induced map has an exact exclusion provider')
        for q in saved:
            mp=q['images'];missing=sorted(set(range(13))-set(mp))
            need(q['missing']==missing and len(missing)==1,'exact required missing target node')
            target=c.W.classical_graph('D',13)
            contacts=[base[i]for i,t in enumerate(mp)if target[missing[0]].get(t,0)]
            need(contacts==[q['forced_contact']],'forced old contact is the actual full-map target neighbor')
            contact=q['forced_contact'];need(contact in high['possible_contacts'],'forced contact belongs to actual contact domain')
            need(core['inverse']is not None and F(q['projection_lower'])==F(core['inverse'][contact][contact])>2,
                 'strict source projection bound joined to actual core inverse entry')
        result['conclusion']='Every listed proper completion is impossible without an independently paid short8 escape.'
        result['computed_premise']='Saved inverse/contact/projection entries are correct; their accepted geometric proof applies to every native completion.'
        return result
    need(len(base)==13,'whole guarded selected factor has order13')
    need(method in('radical_full_core_seal','untouched_small_factor','one_contact_inverse_rank_bound'),
         'recognized exact whole guarded source proof')
    # The inverse-rank branch allows the displayed component to grow. Its
    # source calls this old_small_factor; it is not assumed already whole.
    small=p['old_small_factor']if method=='one_contact_inverse_rank_bound'else p['whole_small_factor']
    nodes=small['nodes']
    parts=c.E.components(g,sorted(set(white)|set(marks)),c.FUEL)
    need(nodes in parts and not set(nodes)&set(marks)and 1<=len(nodes)<=4,'actual displayed disjoint nonempty small factor')
    maximum=p['maximum_small_rank'];need(type(maximum)is int and len(nodes)<=maximum<=4,'source-proved whole-small-factor rank bound')
    if method=='radical_full_core_seal':
        need(high['method']=='radical','exact radical branch source domain')
    else:
        need(core['inverse']is not None,'actual inverse data required for nonradical guard branches')
        free=sorted(set(high['possible_contacts'])-set(base))
        expected=[[i,str(F(core['inverse'][i][i]))]for i in free]
        actual=[[i,str(F(x))]for i,x in p['remaining_inverse_diagonals']]
        need(actual==expected,'every remaining inverse/contact entry joined to actual high domain')
        possible={i for i in free if F(core['inverse'][i][i])<=2}
        need(p['possible_contacts']==sorted(possible),'exact remaining feasible-contact subset, not an arbitrary saved subset')
    if method=='radical_full_core_seal':
        need(p['remaining_positions']==sorted(set(high['possible_contacts'])-set(base)),
             'radical free positions joined to actual high contact domain')
        need(p['radical']==core['radical'][0],'radical vector joined to actual physical core')
    elif method=='untouched_small_factor':
        need(not set(nodes)&possible,'whole-small-factor seal avoids every remaining feasible contact')
    else:
        rho=F(p['inverse_diagonal']);extra=p['maximum_added_rank']
        contact=p['contact']
        need(contact in nodes and set(nodes)&possible=={contact},'unique actual small-factor contact')
        need(rho==F(core['inverse'][contact][contact]),'one-contact rho joined to actual core inverse')
        need(rho>1,'one-contact strict inverse diagonal')
        need(type(extra)is int and extra>=0 and len(nodes)+extra==maximum==4,'one-contact complete attached-rank cap')
        need(extra==int(1/(rho-1)),'exact complete source added-rank cap')
        need(rho>1 and rho*F(extra+1,extra+2)>1,'strict first-forbidden inverse-energy inequality')
        need(F(p['first_forbidden_rank_energy'])==rho*F(extra+1,extra+2),'exact saved forbidden-rank scalar')
    # No subtraction of unsealed displayed ranks. The accepted theorem gives
    # one distinct ACTUAL whole factor of positive rank, bounded above by4.
    # It is neither D13 nor a rank15 witness; it therefore exceeds the already
    # full58-rank budget under the guard, regardless of other spectators.
    need(13+3*15+1>58,'actual disjoint whole factors exceed ambient58')
    result.update(whole_small_nodes=nodes,maximum_whole_small_rank=maximum,
                  conclusion='Guard TRUE is impossible: marked D13 plus three rank15 whole factors plus one distinct positive-rank whole small factor exceeds rank58.',
                  computed_premise='Accepted source theorem proves the named physical small factor stays whole of rank at most the stated cap; saved radical/inverse/contact data are correct.')
    return result


def run(shard,c,r):
    eligible,cores,saved,proofs,partition=domain(c,r);target=c.W.classical_graph('D',13)
    need(all(c.W.finite_same_weight(('D',13),*p)==1 for p in PAIRS),'guard-FALSE current ordinary value is one')
    result=[];allmaps=0;affectedmaps=0
    for p in eligible:
        if p['index']%SHARDS!=shard:continue
        core=cores[p['core_key']];nodes=p['marked_component'];marks=p['preferred_corner'];g=c.E.matrix(core['gram'])
        lab=c.E.labels(g,nodes);positions=[nodes.index(m)for m in marks]
        maps=[dict(images=mp,marked_images=sorted(mp[i]for i in positions))for mp in r.embeddings(lab,target,c.FUEL,c.E)]
        need(canon(c,maps)==canon(c,saved[p['index']]['maps']),'exact complete native D13 map-list equality')
        affected=[q for q in maps if tuple(q['marked_images'])in PAIRS]
        need(canon(c,affected)==canon(c,saved[p['index']]['affected_maps']),'exact changed-role map-list equality')
        allmaps+=len(maps);affectedmaps+=len(affected)
        if affected:result.append(proof(c,proofs[p['index']],p,core,affected))
    return dict(shard=shard,shards=SHARDS,full_domain801=[p['index']for p in eligible],
                complete_shard_maps=allmaps,complete_shard_affected_maps=affectedmaps,
                guarded_providers=result,guard_false_payment='ordinary current D13 value1',
                scope='Only selected conditional negative-guard exclusions; the unselected universal quiet/one-sided receipts are separate.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',choices=('index',)+tuple(f'guard-{i}'for i in range(SHARDS)),required=True)
    parser.add_argument('--parent-lock',action='store_true',required=True)
    args=parser.parse_args();relative=str(Path(__file__).resolve().relative_to(ROOT));source=raw(relative)
    r,c=load()
    if args.batch=='index':
        eligible,cores,saved,proofs,partition=domain(c,r)
        out=dict(high_indices=[p['index']for p in eligible],guarded_case_keys=sorted(proofs),
                 full_induced_maps=862,changed_maps=76,shards=SHARDS,
                 scope='Pinned801/33/76 DOMAIN INDEX only; no new isomorphism replay in this index batch')
    else:out=run(int(args.batch[-1]),c,r)
    for path in c.READ:need(c.digest(path)==c.PINS[path],'consumed pinned source/input changed')
    need(raw(relative)==source and sha256(raw(ROUTES)).hexdigest()==ROUTES_SHA,'guard/index source changed')
    need(sha256(c.raw(c.MANIFEST)).hexdigest()==c.MANIFEST_SHA,'collector manifest changed')
    need(sha256(c.raw(r.COLLECTOR)).hexdigest()==r.COLLECTOR_SHA,'collector source changed')
    payload=json.dumps(dict(status='PASS_NAMED_SELECTED_D13_GUARD_BINDING',batch=args.batch,
        source_sha256=sha256(source).hexdigest(),route_source_sha256=ROUTES_SHA,
        consumed_inputs={p:c.PINS[p]for p in sorted(c.READ)},fuel_used=c.FUEL[0],result=out,
        scope='Source-geometric and numerical-premise conditional; no unselected/universal, all-terminal or global finite57 PASS'),sort_keys=True)
    need(len(payload.encode())<=MAX_OUTPUT,'guard receipt byte cap; refusal rather than truncation')
    print(payload)


if __name__=='__main__':main()
