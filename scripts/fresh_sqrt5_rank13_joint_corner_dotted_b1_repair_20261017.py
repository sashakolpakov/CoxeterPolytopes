#!/usr/bin/env python3
"""Add every certified dotted-B1 spectator row to the selected-triangle model."""
import argparse,hashlib,json
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
from fractions import Fraction as F
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_joint_corner_dotted_b1_repair_20261017'
PREDECESSOR='fresh_sqrt5_rank13_joint_corner_all_nonade_repair_20261017'
FAMILY='fresh_sqrt5_rank13_dotted_b1_spectators_20261017'
FULL_DUAL_BOUND=B.dual_bound
RECOVER_PRIMAL=B.recover_primal
RECOVER_DUAL=B.recover_dual

def logged_primal(d,x):
    print('EXACT_PRIMAL_RECOVERY_START',flush=True)
    out=RECOVER_PRIMAL(d,x)
    print('EXACT_PRIMAL_RECOVERY_END',None if out is None else out['maximum'],flush=True)
    return out

def logged_dual(d,x):
    print('EXACT_SPARSE_DUAL_RECOVERY_START',flush=True)
    out=RECOVER_DUAL(d,x)
    print('EXACT_SPARSE_DUAL_RECOVERY_END',out['bound'],len(out['face_multipliers']),flush=True)
    return out

def sparse_dual_bound(d,alphas):
    loads=defaultdict(F)
    for i,a in alphas.items():
        assert a>=0
        if a:
            for j,c in d['faces'][i]['coefficients']:
                assert c>0;loads[j]+=a*c
    beta={};maximizers={}
    for j in sorted(loads):
        var=d['variables'][j];assert var['index']==j and var['multiplicity']>0
        vertex=var['vertex'];value=loads[j]/var['multiplicity'];assert value>0
        if value>beta.get(vertex,F(0)):
            beta[vertex]=value;maximizers[vertex]=[j]
        elif value==beta[vertex]:maximizers[vertex].append(j)
    beta=dict(sorted(beta.items()));maximizers={i:maximizers[i]for i in beta}
    numerator=sum((d['faces'][i]['required']*a for i,a in alphas.items()),F(0))
    denominator=sum(beta.values(),F(0));assert denominator>0
    for j,load in loads.items():
        var=d['variables'][j];assert load<=beta[var['vertex']]*var['multiplicity']
    # Every omitted column has zero load, and every beta and multiplicity is nonnegative.
    return numerator/denominator,beta,maximizers,numerator,denominator

def domain():
    prior=B.read(PREDECESSOR)
    for certificate in prior['certificates']:
        if certificate['kind']=='EXACT_RESTRICTED_MODEL_DUAL':
            alpha={i:F(a)for i,a in certificate['face_multipliers']}
            assert sparse_dual_bound(prior,alpha)==FULL_DUAL_BOUND(prior,alpha)
    B.dual_bound=sparse_dual_bound
    for certificate in prior['certificates']:B.verify_certificate(prior,certificate)
    keys=['sources','counts','variables','vertices','faces','spectator_context_rows',
          'double_h4_triangle_rows','all_selected_nonade_triangle_rows']
    d={k:prior[k]for k in keys};family=B.read(FAMILY);graphs=dict(B.C.models())
    lookup={(tuple(v['forest']),B.dump(v['marked_forest_code'])):v['index']for v in d['variables']}
    present={tuple(v['forest'])for v in d['vertices']}
    required={tuple(c['forest'])for c in family['corner_types']}
    entries=[(name,len(g))for name,g in graphs.items()if len(g)<=13]
    def forests(rank,start=0):
        if rank==0:yield ();return
        for j in range(start,len(entries)):
            name,size=entries[j]
            if size<=rank:
                for tail in forests(rank-size,j):yield tuple(sorted((name,)+tail))
    universe=set(forests(13));assert len(universe)==6789 and present|required<=universe
    new_forests=sorted(universe-present)
    for forest in new_forests:
        g={}
        for name in forest:
            h=graphs[name];off=len(g)
            for i,row in h.items():g[off+i]={off+j:m for j,m in row.items()}
        assert len(g)==13;pieces=B.C.pieces(g);orbits=defaultdict(list)
        for pair in combinations(range(13),2):
            code=B.dump(sorted(B.C.code(p,pair)for p in pieces));orbits[code].append(list(pair))
        vertex=len(d['vertices']);variables=[]
        for code,pairs in sorted(orbits.items()):
            index=len(d['variables']);lookup[forest,code]=index;variables.append([index,len(pairs)])
            d['variables'].append(dict(index=index,vertex=vertex,forest=list(forest),marked_forest_code=json.loads(code),
                multiplicity=len(pairs),literal_pairs=pairs))
        assert sum(m for _,m in variables)==78
        d['vertices'].append(dict(index=vertex,forest=list(forest),variables=variables))
    corner_lookup={}
    for i,c in enumerate(family['corner_types']):
        assert c['index']==i
        corner_lookup[i]=lookup[tuple(c['forest']),B.dump(c['marked_forest_code'])]
    face_lookup={(f['required'],tuple(map(tuple,f['coefficients']))):f['index']for f in d['faces']}
    mappings=[]
    for i,row in enumerate(family['face_rows']):
        assert row['index']==i and row['required']==1
        coeffs=Counter()
        for j,m in row['corner_coefficients']:
            assert m>0;coeffs[corner_lookup[j]]+=m
        assert sum(coeffs.values())==4
        coeff=tuple(sorted(coeffs.items()));key=(1,coeff)
        provenance=dict(kind='complete_dotted_b1_spectator_family',source=FAMILY,face=i)
        if key in face_lookup:d['faces'][face_lookup[key]]['provenance'].append(provenance)
        else:
            face_lookup[key]=len(d['faces'])
            d['faces'].append(dict(index=len(d['faces']),required=1,coefficients=[list(x)for x in coeff],provenance=[provenance]))
        mappings.append(dict(source_face=i,face=face_lookup[key]))
    # Check every literal common-spectator context against the supplied deduplication.
    for i,context in enumerate(family['contexts']):
        assert context['index']==i
        actual=Counter(context['corner_type_indices'])
        assert tuple(sorted(actual.items()))==tuple(map(tuple,family['face_rows'][context['face']]['corner_coefficients']))
        assert i in family['face_rows'][context['face']]['contexts']
    d['counts'].update(vertices=len(d['vertices']),variables=len(d['variables']),
        literal_face_rows=d['counts']['literal_face_rows']+len(family['contexts']),distinct_face_constraints=len(d['faces']),
        dotted_b1_spectator_contexts=len(family['contexts']),dotted_b1_distinct_input_rows=len(family['face_rows']),
        dotted_b1_added_whole_vertex_forests=len(new_forests),all_finite_rank13_vertex_forests=len(universe),
        added_forests_not_required_by_prior_or_dotted_faces=len(universe-(present|required)))
    d['dotted_corner_type_mapping']=[[i,j]for i,j in sorted(corner_lookup.items())]
    d['dotted_spectator_face_rows']=mappings
    d['inherited_domain']='The frozen predecessor literal domain is inherited verbatim; any saved exact primal and dual certificates are rechecked before adding a new row. Every new corner-type lookup, context deduplication and whole-forest78-pair budget is reconstructed.'
    for p in[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{B.C.TAG}.py',
             f'results/{PREDECESSOR}.json',f'results/{FAMILY}.json']:
        d['sources'][p]=hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
    return d

def run(verify):
    d=domain();path=ROOT/f'results/{TAG}.json'
    if verify:
        saved=json.loads(path.read_text())
        for k,v in d.items():assert saved[k]==v,k
        for c in saved['certificates']:B.verify_certificate(d,c)
        return saved
    B.recover_primal=logged_primal;B.recover_dual=logged_dual
    solver,certificates=B.optimize(d)
    feasible=any(c.get('reaches_target12',False)for c in certificates)
    excluded=any(c.get('excludes_target12',False)for c in certificates);assert not(feasible and excluded)
    status='PASS_EXACT_DOTTED_B1_REPAIRED_PARTIAL_RULE_AT12'if feasible else'PASS_EXACT_DOTTED_B1_REPAIRED_PARTIAL_MODEL_OBSTRUCTION'if excluded else'UNRESOLVED_DOTTED_B1_REPAIRED_PARTIAL_MODEL'
    d.update(status=status,dimension=13,target='12',solver=solver,certificates=certificates,
        change='All119698 certified dotted-B1 orthogonal spectator contexts, deduplicated into82263 input quad rows, are added together. Every one of all6789 finite rank13 forests receives its complete78-pair budget and free nonnegative corner-orbit variables, including forests occurring in no face row. All prior constraints and variable indices are retained.',
        scope='Exact model with a complete finite rank13 vertex universe and partial face coverage. The stated dotted-B1 spectator batch is complete, but general rank13 quadrilateral and triangle coverage is not. Any primal requires the remaining local domains before a dimension claim; any dual excludes only the stated corner-payment family.')
    path.write_text(json.dumps(d,indent=2)+'\n');return d

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run(a.verify)
    print(d['status'],d['solver'],[(c['kind'],c.get('maximum',c.get('bound')))for c in d['certificates']],flush=True)
