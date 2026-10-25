#!/usr/bin/env python3
"""SOURCE ONLY: exact original/later/newer/universal current-terminal adapters.

Named batches have separate receipts. No combined full-terminal PASS is emitted.
Use the existing guarded bounded runner with its shared parent lock. Numerical
classification, complete geometric routing and whole-component seals are inputs.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
from collections import Counter
import fcntl
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from pathlib import Path
import types

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = 'results/finite57_terminal_adapter_inputs_small_affine_20261025.json'
MANIFEST_SHA = '8a0973d7b7ef08231e199bee9f04f45c6a8dc01ba66ca8de8acb2d8e513ee7b7'
ENGINE = 'scripts/finite57_terminal_literal_engine_small_affine_20261025.py'
LIB = 'scripts/finite57_literal_weights_audit_20261024.py'
PROJECTOR = 'scripts/check_d13_transfer_batches_20261020.py'
HIGH = 'fresh_distance8_high_overlap_20261017'
OLD16 = 'fresh_order16_finite_cross_payments_20261017'
OLD17 = 'fresh_ideal7_finite_cross_payments_20261017'
D13 = 'fresh_affine_d13_long_pair_transfer_20261017'
D14 = 'fresh_affine_d14_long_pair_transfer_20261017'
P17 = 'fresh_d17_partial_weight_transfer_20261017'
P16 = 'fresh_d16_partial_weight_transfer_20261017'
D11 = 'fresh_affine_d11_leaf_pair_transfer_20261017'
B13 = 'fresh_b13_endpoint_retention_20261018'
FINAL = 'fresh_d13_three_rank15_context_transfer_20261018'
COMPOSE = 'fresh_d12_half_face_composition_20261018'
VARIANTS = 'fresh_d14_long_leaf_a16_half_selected_transfer_20261018'
A13 = 'fresh_a13_endpoint_a1_context_retention_20261018'
LATER = {
 'fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018': ('radical_cases','proper_completion_cases'),
 'fresh_finite_d17_five_unit_reduction_20261018': ('radical_transfers',),
 'fresh_rank16_existing_payment_exchange_20261018': ('selected_proofs',),
 'fresh_finite_rank15_existing_half_transfer_20261018': ('selected_proofs',),
 'fresh_finite_b15_d15_existing_payment_transfer_20261018': ('selected_proofs',),
 'fresh_finite_d15_three_singleton_half_transfer_20261018': ('selected_proofs',),
 'fresh_finite_d16_fork_exchange_20261018': ('selected_proofs','proper_completion_proofs'),
 'fresh_finite_d16_second_pair_zero_transfer_20261018': ('radical_proofs','proper_completion_proofs'),
 'fresh_d16_d17_eighth_balance_20261018': ('reweighted_proofs',),
 'fresh_finite_d16_alternate_half_zero_transfer_20261018': ('selected_proofs',),
 'fresh_finite_d17_leaf_a16_half_transfer_20261018': ('selected_proofs',),
 'fresh_finite_d17_a4_d13_half_transfer_20261018': ('selected_proofs',),
 'fresh_finite_d17_d5_d12_unit_transfer_20261018': ('selected_proofs',),
 'fresh_a16_a17_full_payment_exchange_20261018': ('selected_proofs',),
 'fresh_a17_five_pair_alternate_transfer_20261018': ('selected_proofs',),
 'fresh_finite_a16_endpoint_sparse_half_transfer_20261018': ('selected_proofs',),
 'fresh_d17_existing_three_corner_transfer_20261018': ('selected_proofs',),
 'fresh_finite_rank14_existing_payment_transfer_20261018': ('selected_proofs',),
 'fresh_finite_a14_near_endpoint_zero_transfer_20261018': ('selected_proofs',),
 'fresh_finite_a15_three_coupon_zero_transfer_20261018': ('selected_proofs',),
 'fresh_finite_d15_a4_a11_half_zero_transfer_20261018': ('selected_proofs',),
 'fresh_finite_d16_d17_two_unit_exchange_20261018': ('radical_exchange_proofs',),
}
NEWER = (
 ('fresh_rank16_a1_trade_composed_ambient58_frontier_20261018',('local_proof','whole_proofs'),False),
 ('fresh_d16_a1_long_leaf_d17_zero_ambient58_frontier_20261018',('local_proof','whole_proofs'),False),
 ('fresh_a16_a1_seven_eighth_second_coupon_20261018',('local_proof','whole_proofs'),False),
 ('fresh_d17_a9_a2_context_zero_ambient58_frontier_20261018',('local_proof','proofs'),False),
 ('fresh_rank13_variable_context_retention_20261018',('local_seals',),False),
 ('fresh_a14_endpoint_a1_context_half_transfer_20261018',('local_proof',),True),
 ('fresh_d17_leaf_a16_context_coupled_transfer_20261018',('local_proof','old_ideal17_seals'),False),
 ('fresh_d17_leaf_a16_context_coupled_transfer_20261018',('local_proof','A16_seals'),False),
)
BATCHES = ('index',) + tuple(f'{family}-{i}' for family in ('original','later','newer') for i in range(4)) + (
 'a13-guards','d14-variants','selected-d15','selected-e6','selected-d17-contexts',
 'universal-bc13','universal-d13','universal-d14-near','universal-d15-basic',
 'universal-d11-repairs','universal-d13-repairs','final-b13-repairs',
 'final-one-sided','final-quiet')
PINS = {}; READ = set(); W = E = J = None; FUEL = [0]; WITNESS_POOL = None
MAX_INPUT_BYTES=40*1024**2
MAX_OUTPUT_BYTES=500*1024


def need(ok, message):
    if not ok:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'))


def raw(relative, cap=2*1024**2):
    path=(ROOT/relative).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(),'bounded source path')
    with path.open('rb') as stream: data=stream.read(cap+1)
    need(0<len(data)<=cap,'bounded source read')
    return data


def digest(relative):
    path=(ROOT/relative).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(),'pinned input path')
    need(0<path.stat().st_size<=MAX_INPUT_BYTES,'pinned input byte cap')
    h=sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(65536),b''):h.update(block)
    return h.hexdigest()


def unique(pairs):
    result={}
    for key,value in pairs:
        need(key not in result,'duplicate JSON key')
        result[key]=value
    return result


def module(relative,name):
    data=raw(relative)
    need(sha256(data).hexdigest()==PINS[relative],'module source pin')
    obj=types.ModuleType(name);obj.__file__=str(ROOT/relative)
    exec(compile(data,str(ROOT/relative),'exec'),obj.__dict__)
    READ.add(relative)
    return obj


def document(stem,*paths):
    relative=f'results/{stem}.json'
    need(relative in PINS and digest(relative)==PINS[relative],'scientific input pin: '+relative)
    READ.add(relative)
    return J.JsonProjection(relative).load(paths)


def at(doc,path):
    for key in path:doc=doc[key]
    return doc


def field(stem,path):return at(document(stem,path),path)


def key(stem,index):return stem,canonical(index)


def part(values,shard):return [value for i,value in enumerate(values) if i%4==shard]


def gram(n,edges):
    g=[['2' if i==j else '0' for j in range(n)]for i in range(n)]
    for edge in edges:
        i,j,*magnitude=edge;q=str(magnitude[0]) if magnitude else '1'
        need(i!=j and g[i][j]=='0','duplicate constructed Gram edge')
        g[i][j]=g[j][i]='-('+q+')'
    return g


def affine(kind,n):
    graph=W.classical_graph(kind,n,affine=True)
    return gram(n,[(i,j,'sqrt(2)' if label==4 else '1')
                   for i,row in enumerate(graph)for j,label in row.items()if i<j])


def extend(base,extra,edges):
    n=len(base);out=[row+['0']*extra for row in base]
    out.extend([['2' if i==j else '0' for j in range(n+extra)]for i in range(n,n+extra)])
    for i,j,*mag in edges:
        q=str(mag[0])if mag else'1';need(not E.number(out[i][j]).terms,'extension edge already present')
        out[i][j]=out[j][i]='-('+q+')'
    return out


def d13_gram(row):
    left,right=row['leaf_attachments']
    edges=[(i,i+1)for i in range(9)]+[(0,10),(0,11),(9,12),(9,13)]
    return gram(19,edges+[(14,left),(15,right),(16,row['A3_attached_to']),(16,17),(17,18)])


def original():
    core=field(HIGH,('isolation','cores','(25, 25, ((9, 9),))'))['gram']
    a=field(OLD16,('exceptional_D15',))
    b=field(OLD17,('order17',))
    for stem,g,rows in ((OLD16,a['gram'],a['cases']),(OLD17,core,b['replacement_cases'])):
        for row in rows:
            white=sorted(set(range(len(g)))-set(row['black']))
            yield stem,row['case'],g,white,row['corners'],'stored-original'
    for stem,selector in ((D13,('extreme_equality','rows')),(D14,('equality','rows'))):
        for i,row in enumerate(field(stem,selector)):
            g=d13_gram(row)if stem==D13 else row['gram']
            yield stem,i,g,row['white'],row['corners'],'d13_gram(recipe)'if stem==D13 else'stored-equality'
    for row in field(P17,('cases',)):
        yield P17,row['index'],row['gram'],row['white'],row['corners'],'stored-partial17'
    for row in field(P16,('cases',)):
        for i,proof in enumerate(row['completions']):
            if 'corners'in proof:
                yield P16,[row['index'],i],proof['gram'],proof['white'],proof['corners'],'stored-completed16'
    for i,row in enumerate(field(D11,('two_A1_A3_spectators','rows'))):
        yield D11,i,row['gram'],row['white'],row['corners'],'stored-universal-D11'


def value(g,white,marks,extra=(),complete=False,sealed_witnesses=()):
    if not complete:
        known=witness_pool().get(context_key(g,white,marks),[])
        sealed_witnesses=[list(p)for p in sorted({tuple(sorted(p))for p in list(sealed_witnesses)+known})]
    return E.corner(g,white,marks,W,FUEL,extra,complete,sealed_witnesses)


def context_key(g,white,marks):
    parsed=E.matrix(g)
    normalized=[[[[d,str(c)]for d,c in sorted(x.terms.items())]for x in row]for row in parsed]
    return canonical([normalized,sorted(white),sorted(marks)])


def same_images(encoded,white,marks,kind,order):
    """Exact complete same-factor images, independently of a saved role name."""
    g=E.matrix(encoded);ids=sorted(set(white)|set(marks))
    parts=E.components(g,ids,FUEL)
    shared=next((p for p in parts if set(marks)<=set(p)),None)
    if shared is None or len(shared)!=order:return []
    lab=E.labels(g,shared);target=W.classical_graph(kind,order)
    return [list(p)for p in E.embeddings(lab,target,[shared.index(m)for m in marks],FUEL)]


def current_collector_join(g,white,marks,before,after):
    """Bind both final role-consumers, not merely their occurrence counts."""
    need(before['marks']==after['marks']==marks,'final collector ordered mark association')
    b=same_images(g,white,marks,'B',13);d=same_images(g,white,marks,'D',13)
    need(before['images']==b and after['images']==d,'final B13/D13 canonical image equality')
    need([0,12]not in b and not any(p in([0,11],[0,12])for p in d),
         'collector domain unexpectedly contains a final changed role')
    return dict(B13_images=b,D13_images=d,
                scope='Exact final role-consumer join; not a substitute for current payment recomputation')


WITNESS_FIELDS={'whole_A1_witnesses','whole_A1','whole_A2','whole_A3','whole_A9','whole_A10','whole_A11',
                'whole_witnesses','whole_unmarked_witnesses'}


def declared_witnesses(row):
    """Only explicit whole-witness metadata; a generic components list is ignored."""
    out=[]
    def take(x):
        if isinstance(x,dict):
            if 'nodes'in x:
                need(type(x['nodes'])is list,'whole-witness node schema')
                out.append(x['nodes'])
            else:
                for child in x.values():take(child)
        elif isinstance(x,list):
            for child in x:take(child)
        else:raise ValueError('whole-witness metadata is not a component collection')
    def find(x):
        E.charge(FUEL)
        if isinstance(x,dict):
            for name,child in x.items():
                if name in WITNESS_FIELDS:take(child)
                elif name not in('old_seal','full_seal','old_complete_seal','full_old_seal'):
                    find(child)
        elif isinstance(x,list):
            for child in x:find(child)
    find(row)
    return [list(p)for p in sorted({tuple(sorted(p))for p in out})]


def all_witnesses(g,white,corners,geometric_premise):
    """Explicit FULL-core-seal premise, never inferred from a paying-factor seal."""
    need(type(geometric_premise)is str and geometric_premise,'named full-core-seal source premise required')
    physical=E.matrix(g);out={}
    for row in corners:
        if row is None:continue
        marks=row['marks'];parts=E.components(physical,sorted(set(white)|set(marks)),FUEL)
        out[tuple(sorted(marks))]=[p for p in parts if not set(p)&set(marks)]
    return out


def witness_pool():
    """Join explicit sealed-witness fields to exact physical corner contexts.

    No component is made whole merely because it is displayed or classified.
    The sole complete-forest insertion is the four separately accepted
    rank13 FULL-core seals. All other masks come from named whole-witness
    output fields in their original proof records.
    """
    global WITNESS_POOL
    if WITNESS_POOL is not None:return WITNESS_POOL
    WITNESS_POOL={}
    def add(g,w,marks,nodes):
        key=context_key(g,w,marks)
        physical=E.matrix(g);parts=E.components(physical,sorted(set(w)|set(marks)),FUEL)
        for p in nodes:need(sorted(p)in parts and not set(p)&set(marks),'whole witness metadata must be an unmarked full physical component')
        old=WITNESS_POOL.get(key,[])
        WITNESS_POOL[key]=[list(p)for p in sorted({tuple(sorted(p))for p in old+nodes})]
    for stem,path,single in NEWER:
        found=field(stem,path)
        for p in [found]if single else found:
            g,w,gp,wp=old_context(p)
            need(g is not None and w is not None,'witness source physical context')
            if stem=='fresh_rank13_variable_context_retention_20261018':
                need(p['name']in('order16','cycle14','D15_1147','D15_1269'),'exact four rank13 full-core seal premises')
                full=all_witnesses(g,w,p['corners'],'accepted rank13 local full-core seals: radical or no remaining contact below norm2')
            elif stem=='fresh_d17_a9_a2_context_zero_ambient58_frontier_20261018':
                need(p['index']in(1062,2525),'exact two A9/A2 full-core D17 cases')
                full=all_witnesses(g,w,p['corners'],'accepted two full radical D17 seals; whole A9/A2 witnesses at the A1/D5 corner')
            else:full={}
            for c in p['corners']:
                add(g,w,c['marks'],declared_witnesses(c)+full.get(tuple(sorted(c['marks'])),[]))
    a13=field(A13,('local_proof',))
    for p in a13['radical_proofs']:
        add(p['gram'],p['white'],p['selected_marks'],[b['nodes']for b in p['whole_A1_witnesses']])
    p=a13['saturated_two_witness_proof']
    for v in p['variants']:add(v['gram'],v['white'],p['selected_marks'],[b['nodes']for b in v['whole_A1_witnesses']])
    guards=field('fresh_d14_small_component_context_retention_20261018',('selected_whole_variants',))
    physical={(p['index'],canonical(v['extra_contact'])):v for p in field(VARIANTS,('selected_proofs',))for v in p.get('variants',[])}
    need(len(guards)==len(physical)==24,'complete24 D14 witness join')
    for p in guards:
        v=physical[p['index'],canonical(p['extra_contact'])]
        for c in p['guarded_outputs']:add(v['gram'],v['white'],c['marks'],[b['nodes']for b in c['whole_witnesses']])
    # The historical universal D13 equalities have a separately proved wholeA3.
    hist=field('fresh_d14_small_component_context_retention_20261018',('historical_interfaces','records'))
    originals={key(s,i):(g,w)for s,i,g,w,cs,p in original()}
    for p in hist:
        g,w=originals[key(p['source'],p['case'])]
        for c in p['guarded_outputs']:add(g,w,c['marks'],[b['nodes']for b in c['whole_witnesses']])
    # Only these four explicitly reviewed radical full-core models supply
    # the A11 witnesses for the old long/short D5/D12 exchange. This is not
    # a generic inference from displayed component size.
    stem='fresh_d5_d12_alternate_coupon_rebalance_20261018'
    for p in field(stem,('local_proofs','proofs')):
        need(p['index']in(1978,2371,2535,2539),'exact four long/short full-core seals')
        old=p['old_seal'];g,w=old['gram'],old['white']
        full=all_witnesses(g,w,p['corners'],'accepted whole19-root affine/radical seal for four D5/D12 long/short cases')
        for c in p['corners']:add(g,w,c['marks'],full[tuple(sorted(c['marks']))])
    stem='fresh_order17_center_contextual_zero_transfer_20261018'
    for p in field(stem,('local_proof','proofs')):
        old=p['inherited_full_seal']
        for c in p['corners']:add(old['gram'],old['white'],c['marks'],declared_witnesses(c))
    return WITNESS_POOL


def receipt(stem,index,g,white,corners,provenance,lower=None,extra=(),complete=False,
            sealed_witnesses=None):
    need(type(corners)is list and corners,'nonempty actual corner list')
    out=[];seen=set();pointwise_shortfalls=[]
    for i,row in enumerate(corners):
        if row is None:continue
        marks=row.get('marks',row.get('preferred_corner'))
        need(type(marks)is list and len(marks)==2,'actual corner marks')
        unordered=tuple(sorted(marks));need(unordered not in seen,'corner counted twice')
        witnesses=declared_witnesses(row)if sealed_witnesses is None else sealed_witnesses.get(unordered,[])
        seen.add(unordered);actual=value(g,white,marks,extra,complete,witnesses)
        if lower is not None:
            old=lower[i]
            need(old['marks']==marks,'literal ledger corner association')
            chosen=F(old['new_lower']);need(chosen>=0,'saved legacy lower is nonnegative')
            actual['saved_legacy_lower']=str(chosen)
            if chosen>F(actual['weight']):
                pointwise_shortfalls.append(dict(source=stem,index=index,marks=marks,
                    saved_legacy_lower=str(chosen),conservative_current_lower=actual['weight'],
                    literal_corner=actual,status='UNUSED_IN_DIRECT_RECEIVER_PROOF'))
        out.append(actual)
    total=sum(F(r['weight'])for r in out)
    need(total>=1,('current distinct-corner receipt below one',stem,index,str(total)))
    return dict(source=stem,index=index,provenance=provenance,
                context_sha256=sha256(canonical([g,white]).encode()).hexdigest(),
                corners=out,total=str(total),direct_current_receiver_condition=True,
                legacy_pointwise_bounds_used=False,legacy_pointwise_shortfalls=pointwise_shortfalls,
                scope='Direct CURRENT distinct-corner total >=1; saved legacy pointwise lowers are diagnostic only; geometric corner/whole-component premises retained')


def old_context(x,g=None,white=None,gpath=None,wpath=None,path=''):
    for name in ('old_seal','old_complete_seal','full_seal','full_old_seal','full_old_complete_seal','inherited_full_seal'):
        inherited=x.get(name)
        if isinstance(inherited,dict):
            g,white,gpath,wpath=old_context(inherited,g,white,gpath,wpath,path+'/'+name)
    for name in ('core_gram','gram'):
        if name in x:g=x[name];gpath=path+'/'+name
    if 'white'in x:white=x['white'];wpath=path+'/white'
    return g,white,gpath,wpath


def walk(x,path='',g=None,white=None,gpath=None,wpath=None,depth=0):
    E.charge(FUEL);need(depth<=64,'collector nesting cap')
    if isinstance(x,dict):
        g,white,gpath,wpath=old_context(x,g,white,gpath,wpath,path)
        yield path,x,g,white,gpath,wpath
        for name,child in x.items():
            gg,ww,gp,wp=g,white,gpath,wpath
            if name in('after','extended_corners'):
                if 'extended_gram'in x:gg=x['extended_gram'];gp=path+'/extended_gram'
                if 'extended_white'in x:ww=x['extended_white'];wp=path+'/extended_white'
            yield from walk(child,path+'/'+name,gg,ww,gp,wp,depth+1)
    elif isinstance(x,list):
        for i,child in enumerate(x):yield from walk(child,path+'/'+str(i),g,white,gpath,wpath,depth+1)


COMPLETE_CORNER_LISTS=('corners','new_corners','before','after','extended_corners')


def complete_models(x,path='',g=None,white=None,gpath=None,wpath=None,depth=0):
    """Every explicitly named complete receipt list, retaining its actual Gram.

    These names come from the22 pinned source packets. A list is not made a
    complete receipt merely because a dictionary somewhere contains marks.
    The caller must check every required interface path lies in a returned
    list, and every returned list has three/four distinct actual corner pairs.
    The source theorem supplies that this list funds its receiving face at
    threshold one. This function does not prove the separate triangle bound2.
    """
    E.charge(FUEL);need(depth<=64,'complete model nesting cap')
    if isinstance(x,dict):
        g,white,gpath,wpath=old_context(x,g,white,gpath,wpath,path)
        for name,child in x.items():
            gg,ww,gp,wp=g,white,gpath,wpath
            if name in('after','extended_corners'):
                if 'extended_gram'in x:gg=x['extended_gram'];gp=path+'/extended_gram'
                if 'extended_white'in x:ww=x['extended_white'];wp=path+'/extended_white'
            if name in COMPLETE_CORNER_LISTS and isinstance(child,list)and child and all(isinstance(c,dict)and'marks'in c for c in child):
                need(gg is not None and ww is not None,'complete receipt physical context')
                need(len(child)in(3,4),'only explicit three/four-corner receipt lists may fund an interface')
                need(len({tuple(sorted(c['marks']))for c in child})==len(child),'complete receipt has repeated corner')
                yield path+'/'+name,gg,ww,child,gp,wp
            yield from complete_models(child,path+'/'+name,gg,ww,gp,wp,depth+1)
    elif isinstance(x,list):
        for i,child in enumerate(x):yield from complete_models(child,path+'/'+str(i),g,white,gpath,wpath,depth+1)


def original_batch(shard):
    ledger=field(B13,('current_interfaces','literals'))
    saved={key(row['source'],row['case']):row for row in ledger}
    need(len(saved)==len(ledger)==216,'unique complete216 ledger keys')
    records=list(original());seen=[key(s,i)for s,i,*_ in records]
    need(len(seen)==len(set(seen))==216 and set(seen)==set(saved),'exact216 physical/ledger join')
    need(sum(len(r[4])for r in records)==832,'complete832 ordered corners')
    final=field(FINAL,('current_interfaces','literal_corner_rows'))
    final_by_key={(r['source'],canonical(r['case']),tuple(r['marks'])):r for r in final}
    wanted={(s,canonical(i),tuple(c['marks']))for s,i,g,w,cs,p in records for c in cs}
    need(len(final_by_key)==len(final)==832 and set(final_by_key)==wanted,'exact final832 source/case/mark join')
    result=[]
    for stem,index,g,white,corners,provenance in part(records,shard):
        before=saved[key(stem,index)]
        need(len(before['corners'])==len(corners),'original corner-domain identity')
        row=receipt(stem,index,g,white,corners,dict(origin=provenance,whole_witness_premise='exact separately indexed source-declared witness fields'),before['corners'])
        row['final_collector_joins']=[current_collector_join(g,white,c['marks'],b,
            final_by_key[stem,canonical(index),tuple(c['marks'])])for c,b in zip(corners,before['corners'],strict=True)]
        result.append(row)
    return dict(shard=shard,rows=result,complete_domain_keys=seen)


def later_batch(shard):
    expected=field(COMPOSE,('later_finite_output_interfaces',))
    reference={row['source']:row for row in expected}
    need(set(reference)==set(LATER)and len(expected)==22,'exact22 expected packet domain')
    b13=field(B13,('current_interfaces','later_corners'))
    d13=field(FINAL,('current_interfaces','later_corners'))
    bindex={(r['source'],r['path']):r for r in b13}
    dindex={(r['source'],r['path']):r for r in d13}
    wanted={(r['source'],c['path'])for r in expected for c in r['corners']}
    need(len(bindex)==len(b13)==len(dindex)==len(d13)==564 and set(bindex)==set(dindex)==wanted,
         'exact final564 source/path joins')
    result=[]
    for stem in part(list(LATER),shard):
        doc=document(stem,*[(section,)for section in LATER[stem]])
        keys={};models={};covered={}
        for section in LATER[stem]:
            for path,row,g,white,gp,wp in walk(doc[section],section):
                if 'marks'in row and ('components'in row or'whole_components'in row):
                    need(g is not None and white is not None,'later physical context missing')
                    need(path not in keys,'duplicate later path')
                    got=value(g,white,row['marks'])
                    keys[path]=dict(marks=row['marks'],literal=got,gram_path=gp,white_path=wp,
                        context_sha256=sha256(canonical([g,white]).encode()).hexdigest(),
                        final_collectors=current_collector_join(g,white,row['marks'],bindex[stem,path],dindex[stem,path]))
            for path,g,white,cs,gp,wp in complete_models(doc[section],section):
                name=canonical([g,white,[c['marks']for c in cs]])
                models.setdefault(name,(path,g,white,cs,gp,wp))
                for i,corner in enumerate(cs):
                    occurrence=path+'/'+str(i)
                    covered[occurrence]=dict(model_key=sha256(name.encode()).hexdigest(),marks=corner['marks'],
                        context_sha256=sha256(canonical([g,white]).encode()).hexdigest())
        wanted={r['path']:r for r in reference[stem]['corners']}
        need(len(wanted)==len(reference[stem]['corners']) and set(keys)==set(wanted),'exact later path domain')
        for path,row in keys.items():need(row['marks']==wanted[path]['marks'],'later key/mark association')
        need(set(keys)<=set(covered),('required later interface lacks a full current-funded polygon list',stem,sorted(set(keys)-set(covered))))
        for path,row in keys.items():
            cover=covered[path]
            need(cover['marks']==row['marks'],'funded-list/required-interface ordered mark association')
            need(cover['context_sha256']==row['context_sha256'],
                 'funded-list/required-interface exact physical context')
        funded=[]
        for name,(path,g,w,cs,gp,wp)in models.items():
            funded.append(dict(model_key=sha256(name.encode()).hexdigest(),
                receipt=receipt(stem,path,g,w,cs,dict(gram=gp,white=wp))))
        need(funded,'later packet has no complete physical receipts')
        result.append(dict(source=stem,interfaces=keys,required_path_to_funded_model={p:covered[p]for p in keys},receipts=funded))
    return dict(shard=shard,packets=result,complete_packet_keys=list(LATER),complete_expected_occurrences=564)


def newer_batch(shard):
    reference=field(B13,('current_interfaces','newer_corners'))
    rows=[];domain=[]
    for stem,path,single in NEWER:
        proofs=field(stem,path)
        for pindex,p in enumerate([proofs]if single else proofs):
            identifier=p.get('index',p.get('case',p.get('name','A14')))
            base='/'+ '/'.join(path)+(''if single else'/'+str(pindex))
            g,w,gp,wp=old_context(p,path=base)
            need(g is not None and w is not None,'newer physical seal missing')
            for corner in p['corners']:domain.append((stem,canonical(identifier),tuple(corner['marks'])))
            rows.append((stem,identifier,g,w,p['corners'],dict(gram=gp,white=wp)))
    expected=[(r['source'],canonical(r['index']),tuple(r['marks']))for r in reference]
    need(len(domain)==len(set(domain))and set(domain)==set(expected)and len(expected)==len(domain),'exact newer key domain')
    result=[]
    for row in part(rows,shard):
        stem,identifier,g,w,cs,provenance=row
        witnesses=None
        if stem=='fresh_rank13_variable_context_retention_20261018':
            need(identifier in('order16','cycle14','D15_1147','D15_1269'),'exact four full-core rank13 seals')
            premise='rank13 variable local four full-core seals, inherited from complete original radical/inverse exclusions'
            witnesses=all_witnesses(g,w,cs,premise);provenance=dict(provenance,whole_witness_premise=premise)
        result.append(receipt(stem,identifier,g,w,cs,provenance,sealed_witnesses=witnesses))
    return dict(shard=shard,receipts=result,complete_newer_keys=domain)


def variants():
    expected=field(B13,('current_interfaces','D14_variants'))
    wanted={(r['index'],canonical(r['extra_contact'])):r for r in expected}
    result=[]
    for proof in field(VARIANTS,('selected_proofs',)):
        for v in proof.get('variants',[]):
            k=proof['index'],canonical(v['extra_contact'])
            need(k in wanted,'unexpected D14 variant')
            result.append(receipt(VARIANTS,list(k),v['gram'],v['white'],v['corners'],'stored-D14-variant'))
    keys={(r['index'][0],r['index'][1])for r in result}
    need(len(result)==len(wanted)==24 and keys==set(wanted),'exact24 D14 variant domain')
    return dict(receipts=result)


def a13_guards():
    proof=field(A13,('local_proof',));result=[]
    for p in proof['radical_proofs']:
        g,w,gp,wp=old_context(p);need(g is not None and w is not None,'A13 radical physical context')
        row=dict(marks=p['selected_marks'],whole_A1_witnesses=p['whole_A1_witnesses'])
        result.append(receipt(A13,p['index'],g,w,[row],dict(gram=gp,white=wp)))
    p=proof['saturated_two_witness_proof']
    for v in p['variants']:
        result.append(receipt(A13,[p['index'],v['extra_contact']],v['gram'],v['white'],
                              [dict(marks=p['selected_marks'],whole_A1_witnesses=v['whole_A1_witnesses'])],'saturated-A13-context'))
    need({r['index']if isinstance(r['index'],int)else r['index'][0]for r in result}=={623,1349,1554,1627,1640},'exact five A13 source cases')
    return dict(receipts=result)


def selected_d15():
    stem='fresh_order15_affine_d_half_coupon_transfer_20261017'
    result=[];domain=set();whole=set()
    for p in field(stem,('cases',)):
        for v in p['proof'].get('variants',[]):
            k=p['index'],canonical(v['optional_core_contact'])
            need(k not in domain,'duplicate selectedD15 physical variant key')
            domain.add(k);whole.add(p['index'])
            result.append(receipt(stem,[p['index'],v['optional_core_contact']],v['gram'],v['white'],v['corners'],'stored-selected-D15-variant'))
    need(len(result)==20 and whole=={691,795,825,827,1874,1976,2267,2369,2481,2485,2506,2513},
         'complete20 selectedD15 whole variants over exact12 whole cases')
    return dict(receipts=result)


def selected_e6():
    stem='fresh_order17_affine_d_three_quarter_transfer_20261017';result=[]
    for p in field(stem,('cases',)):
        if p['index']not in(1981,2374,2549,2551):continue
        q=p['proof'];result.append(receipt(stem,p['index'],q['core_gram'],q['white'],q['corners'],'stored-whole-E6-D11'))
    need({r['index']for r in result}=={1981,2374,2549,2551},'exact four E6 whole cases')
    return dict(receipts=result)


def selected_d17():
    stem='fresh_d5_d12_alternate_coupon_rebalance_20261018';result=[]
    for p in field(stem,('local_proofs','proofs')):
        q=p['old_seal'];result.append(receipt(stem,p['index'],q['gram'],q['white'],p['corners'],'final-eighth-half-three-eighth'))
    center='fresh_order17_center_contextual_zero_transfer_20261018'
    for p in field(center,('local_proof','proofs')):
        q=p['inherited_full_seal'];result.append(receipt(center,p['index'],q['gram'],q['white'],p['corners'],'current-two-finite-halves'))
    need({r['index']for r in result}=={1978,2371,2535,2539,1314,2529},'exact six D17 contextual cases')
    return dict(receipts=result)


def single_corner(stem,index,g,ids,marks,required,provenance):
    white=sorted(set(ids)-set(marks));r=value(g,white,marks)
    need(F(r['weight'])>=F(required),'universal finite output below current required value')
    return dict(source=stem,index=index,corner=r,required=str(required),provenance=provenance)


def universal(name):
    rows=[]
    if name=='bc13':
        stem='fresh_order13_bc_inner_terminal_transfer_20261018'
        for p in field(stem,('whole_finite_corners',)):
            ids=sorted(set().union(*(set(x)for x in p['components'])))
            rows.append(single_corner(stem,p['profile'],p['gram'],ids,p['paying_corner'],p['payment'],'stored-full-BC13-graph'))
        need({r['index']for r in rows}=={'isolated','mate3'},'complete BC13 finite profiles')
    elif name=='d13':
        stem='fresh_order13_d_center_leaf_half_transfer_20261018'
        for p in field(stem,('finite_corner_representatives',)):
            leaf=p['attached_leaf'];extra=int(p.get('spectator_A1',False))
            edges=[(13,leaf)]+([(13,14)]if extra else[])
            g=extend(affine('D',13),1+extra,edges)
            omitted=0 if p['side']=='center'else 11
            ids=[i for i in range(len(g))if i!=omitted]
            rows.append(single_corner(stem,[p['side'],leaf,extra],g,ids,p['marks'],p['existing_payment'],'literal_coupons exact base/leaf recipe'))
        need(len(rows)==6,'complete six D13 center/leaf outputs')
    elif name=='d14-near':
        stem='fresh_order14_d_near_fork_half_transfer_20261018'
        recipes={'near_extreme_empty':(1,10),'far_isolated_empty':(12,12),'far_mate_empty':(12,13)}
        for p in field(stem,('finite_whole_component_corners',)):
            omit,attach=recipes[p['case']]
            need(p['omitted_mark']==omit,'D14 finite recipe association')
            g=extend(affine('D',14),1,[(14,attach)])
            rows.append(single_corner(stem,p['case'],g,[i for i in range(15)if i!=omit],p['marks'],p['alternative_payment'],'finite_literals exact base/attachment recipe'))
        need({r['index']for r in rows}==set(recipes),'complete D14 near finite profiles')
    elif name=='d15-basic':
        stem='fresh_order15_bd_basic_fractional_transfer_20261018'
        for p in field(stem,('D_whole_finite_corners',)):
            ids=sorted(set().union(*(set(b['nodes'])for b in p['components'])))
            rows.append(single_corner(stem,p['case'],p['gram'],ids,p['marks'],p['distinct_payment'],'stored-D15-basic-finite-output'))
        need({r['index']for r in rows}=={'near_empty','far_isolated_empty','far_mate_empty'},'complete D15 basic finite profiles')
    else:raise ValueError('unknown universal physical family')
    return dict(outputs=rows,scope='Finite output bindings only; arbitrary native-profile and spectator proofs retained as premises')


def universal_repairs(kind):
    stem='fresh_d12_universal_ideal_output_repairs_20261018';result=[]
    selector='affine_D11'if kind=='d11'else'affine_D13'
    for i,p in enumerate(field(stem,(selector,'rows'))):
        if kind=='d11':g,w=p['gram'],p['white'];identifier=[p['old_marks'],p['spectator_free_extreme_side'],p['branch']]
        else:
            a,b=p['old_marks'];leaf=p['attachment'];side=p['side']
            edges=[(14,1),(15,leaf)]if side=='left'else[(14,leaf),(15,8)]
            g=extend(affine('D',14),2,edges);w=sorted(set(range(14))-{a,b});identifier=[p['old_marks'],side,leaf]
        result.append(receipt(stem,identifier,g,w,p['corners'],'exact-universal-repair-construction'))
    need(len(result)==(32 if kind=='d11'else 12),'complete universal D12 repair domain')
    need(len({canonical(r['index'])for r in result})==len(result),'duplicate universal repair key')
    return dict(receipts=result)


def b13_repairs():
    result=[]
    for p in field(B13,('whole_component_repairs','seals')):
        q=p['full_old_seal'];result.append(receipt(B13,p['index'],q['gram'],q['white'],p['corners'],'stored-B13-existing-unit-repair'))
    need({r['index']for r in result}=={166,1744},'exact B13 endpoint consumers')
    return dict(receipts=result)


def final_one_sided():
    result=[]
    for p in field(FINAL,('one_sided_whole_models',)):
        need(p['active_attachment']in(9,10),'one-sided active attachment')
        for pair in combinations_with_replacement([('A',15),('D',15)],2):
            # The stored30-root model already contains the active D15. The two
            # inactive whole factors are actual orthogonal witnesses, not padding.
            use=[c for c in p['corners']if c['marks']in(p['quiet_corner'],p['new_paying_corner'])]
            need(len(use)==2,'two distinct one-sided paying corners')
            premise='final one-sided energy-one classification: stored active D15 and two inactive rank15 factors are whole at both finite corners'
            witnesses=all_witnesses(p['gram'],p['white'],use,premise)
            result.append(receipt(FINAL,[p['active_attachment'],pair],p['gram'],p['white'],use,
                                  premise,extra=pair,complete=True,sealed_witnesses=witnesses))
    need(len(result)==6,'two active models times three inactive pairs')
    return dict(receipts=result)


def final_quiet():
    stem='fresh_affine_d12_leaf_pair_transfer_20261017'
    originals=field(stem,('nonextreme_zero_cases',));repairs=field(FINAL,('quiet_repairs',));result=[]
    wanted={r['index']:r for r in repairs};need(len(wanted)==len(repairs)==24,'unique24 quiet repair keys')
    for index,repair in wanted.items():
        p=originals[index];need(p['extreme_profile']=='mate3','quiet repair extreme profile')
        need(p['nonextreme_profile']==repair['profile'],'quiet profile exact join')
        marks=[repair['quiet_corner'],repair['paying_corner']]
        need(marks[0]!=marks[1],'quiet payment counted twice')
        # All guard-active quiet cases have exactly three whole rank15 factors.
        # Eight ordered choices preserve different component copies.
        for c0 in ('A','D'):
            for c1 in ('A','D'):
                for c2 in ('A','D'):
                    extra=[(c0,15),(c1,15),(c2,15)]
                    if repair['method']=='independent_current_affine_B13_priority':
                        actual=value(p['gram'],p['white'],repair['paying_corner'])
                        need(actual['kind']=='ideal'and actual['forest']==[['B',13]],'quiet independent priority must be the actual whole affineB13 factor')
                        need(actual['marked_images']==repair['images'],'quiet priority complete affineB13 marked-image association')
                        result.append(dict(source=FINAL,index=[index,c0,c1,c2],independent_priority='whole-affine-B13 theorem',
                                           physical_corner=actual,current_intrinsic_weight=actual['weight'],
                                           required_provider=dict(theorem='current affineB13 independent whole-face replacement',
                                               kind='B',order=13,marked_images=actual['marked_images'],
                                               quiet_source_index=index,status='UNBOUND_G4_G6_PROVIDER'),
                                           scope='priority theorem is a separate universal receipt, not an intrinsic unit'))
                    else:
                        cs=[dict(marks=m)for m in marks]
                        premise='final quiet-profile zero-energy full-core model and three actual whole rank15 spectators'
                        witnesses=all_witnesses(p['gram'],p['white'],cs,premise)
                        result.append(receipt(FINAL,[index,c0,c1,c2],p['gram'],p['white'],cs,
                                              premise,extra=extra,complete=True,sealed_witnesses=witnesses))
    need(len(result)==192,'all24 quiet keys and eight ordered spectator choices')
    return dict(rows=result,external_current_priority='Independent affineB13 zero-pair replacement must be bound by universal BC13 route adapter; never credit it as an intrinsic one.')


def index():
    original_keys=[(s,i,[c['marks']for c in cs])for s,i,g,w,cs,p in original()]
    need(len(original_keys)==216 and sum(len(r[2])for r in original_keys)==832,'original index domain')
    later=field(COMPOSE,('later_finite_output_interfaces',))
    need({r['source']for r in later}==set(LATER)and sum(len(r['corners'])for r in later)==564,'later exact source/path index')
    return dict(original=original_keys,later=[dict(source=r['source'],paths=[c['path']for c in r['corners']])for r in later],
                newer_selectors=NEWER,batches=BATCHES,scope='Exact collector index; no numerical receipt or full-coverage PASS')


def parent_lock():
    path=ROOT/'formal/.proof-batches/exclusive.lock';need(path.is_file(),'shared lock missing')
    with path.open('r+')as stream:
        try:fcntl.flock(stream.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return
        fcntl.flock(stream.fileno(),fcntl.LOCK_UN)
    raise ValueError('shared parent lock must already be held')


def pending_provider_occurrences(value):
    if isinstance(value,dict):
        here=int(type(value.get('status'))is str and value['status'].startswith('UNBOUND_'))
        return here+sum(pending_provider_occurrences(x)for x in value.values())
    if isinstance(value,list):return sum(pending_provider_occurrences(x)for x in value)
    return 0


def main():
    global PINS,W,E,J
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',choices=BATCHES,required=True)
    parser.add_argument('--parent-lock',action='store_true',required=True)
    args=parser.parse_args();parent_lock()
    data=raw(MANIFEST);need(sha256(data).hexdigest()==MANIFEST_SHA,'adapter input manifest identity')
    meta=json.loads(data,object_pairs_hook=unique);PINS={r['path']:r['sha256']for r in meta['inputs']}
    need(len(PINS)==len(meta['inputs']),'duplicate manifest path')
    W=module(LIB,'finite57_collector_literal');E=module(ENGINE,'finite57_collector_engine');J=module(PROJECTOR,'finite57_collector_projection')
    # This module object is private to the adapter. Strengthen its selected
    # JSON decoder without changing the stdlib module or frozen source bytes.
    J.json=types.SimpleNamespace(loads=lambda encoded:json.loads(encoded,object_pairs_hook=unique,
        parse_constant=lambda x:(_ for _ in()).throw(ValueError('nonfinite JSON number'))))
    primary='manuscripts/finite_volume57.tex'
    need(digest(primary)==PINS[primary]and PINS[primary]==W.TEX_SHA256,'literal/primary source binding')
    READ.add(primary)
    source=raw(str(Path(__file__).resolve().relative_to(ROOT)))
    batch=args.batch
    if batch=='index':out=index()
    elif batch.startswith('original-'):out=original_batch(int(batch[-1]))
    elif batch.startswith('later-'):out=later_batch(int(batch[-1]))
    elif batch.startswith('newer-'):out=newer_batch(int(batch[-1]))
    elif batch=='a13-guards':out=a13_guards()
    elif batch=='d14-variants':out=variants()
    elif batch=='selected-d15':out=selected_d15()
    elif batch=='selected-e6':out=selected_e6()
    elif batch=='selected-d17-contexts':out=selected_d17()
    elif batch in('universal-d11-repairs','universal-d13-repairs'):out=universal_repairs(batch.split('-')[1])
    elif batch.startswith('universal-'):out=universal(batch[len('universal-'):])
    elif batch=='final-b13-repairs':out=b13_repairs()
    elif batch=='final-one-sided':out=final_one_sided()
    elif batch=='final-quiet':out=final_quiet()
    else:raise ValueError('unimplemented named adapter batch')
    for path in READ:need(digest(path)==PINS[path],'consumed source/input changed during batch')
    need(raw(str(Path(__file__).resolve().relative_to(ROOT)))==source,'adapter source changed during batch')
    need(sha256(raw(MANIFEST)).hexdigest()==MANIFEST_SHA,'input manifest changed during batch')
    pending=pending_provider_occurrences(out)
    payload=json.dumps(dict(status='PASS_NAMED_SCOPE_WITH_UNBOUND_PROVIDERS'if pending else'PASS_NAMED_CURRENT_TERMINAL_ADAPTER',
          pending_provider_occurrences=pending,current_terminal_bridge_complete=False,
          batch=batch,source_sha256=sha256(source).hexdigest(),
          fuel_used=FUEL[0],consumed_inputs={p:PINS[p]for p in sorted(READ)},result=out,
          scope='Only the named collector/model bridge. Complete route-exhaustion, remaining universal outputs, geometry and computed input validity are separate; no full finite57 PASS.'),sort_keys=True)
    need(len(payload.encode())<=MAX_OUTPUT_BYTES,'named receipt output cap; refuse rather than truncate')
    print(payload)


if __name__=='__main__':main()
