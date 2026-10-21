#!/usr/bin/env python3
"""All represented spectator contexts for the raw23/raw25 forced-quad rule.

Checks literal actual provenance maps, one selected cap per profile/coface,
all target forest labels, and a simultaneous monotone vertex-mass repair.
No optimizer and no repeated large local/joint scan.
"""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_sqrt5_multispectator_cap_capacity_20261017 as M

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_core23_25_spectator_incidence_20261017'
BASE='fresh_sqrt5_core23_25_incidence_repair_20261017'
GATE='fresh_sqrt5_core23_25_forced_quad_20261017'
FORESTS='fresh_sqrt5_rank13_combined_disjoint_ports_20261017'
FULL='fresh_sqrt5_cross_coface_cap_full_model_20261017'


def run():
    D={key:M.read(tag)for key,tag in M.TAGS.items()if key!='model'}
    S,I,A,Q,R=[D[k]for k in ('slices','incidence','actual','quad','raw')]
    base,gate,forests,full=map(M.read,[BASE,GATE,FORESTS,FULL])
    assert gate['status']=='PASS_EXACT_CORE23_25_FORCED_QUADRILATERAL_INCIDENCE'
    alpha={i:F(a)for i,a in base['face_multipliers']};mu={i:F(a)for i,a in base['vertex_multipliers']}
    assert all(a>=0 for a in alpha.values())and all(a>=0 for a in mu.values())
    pmap=dict(full['prepared']['provenance_to_face']);assert len(set(pmap.values()))==len(pmap)
    profiles={p['base_type']:p for p in S['coface_slice_demands']}
    context={row[0]:row for row in I['contexts']};assert len(context)==len(I['contexts'])
    templates={t['canonical_core']:t for t in I['prism_templates']}
    frozen_profiles={p['raw_source']:p for p in gate['selected_profiles']}
    vlookup={tuple(v['forest']):v['index']for v in forests['vertices']}
    groups=defaultdict(dict);selected=[]
    for B in S['actual_bases']:
        t=B['base_type'];typ=I['base_types'][t];Fst=typ['spectator_forest'];k=Fst.count('A1')
        if typ['core']!=1 or k<2:continue
        pr=B['provenance'];prov=A['provenance_masses'][pr];ci=prov['canonical_prism_core']
        C=Q['canonical_cores'][ci];raw_source=C['source_core']
        if raw_source not in (23,25):continue
        frozen=frozen_profiles[raw_source];T=next(z for z in templates[ci]['bases']if z['deleted_normal']==B['deleted'])
        assert B['coface_location']==['finite','A1',0] and T['finite_extra']==['A1']
        assert ci==frozen['canonical_prism_core'] and B['deleted']==frozen['deleted']and B['coface_root']==frozen['coface']
        assert T['core_isomorphism']==frozen['core_isomorphism']and B['canonical_vertices']==frozen['canonical_base_vertices']
        assert T['base_core']==1 and T['opposite_base_facet']==B['coface_root']
        assert B['face']==pmap[pr] and B['prism_dimension']==prov['prism_dimension']==len(B['canonical_vertices'])
        assert Counter(Fst)['A1']==k and sum(map(M.V.rank,Fst))==9
        assert profiles[t]['core']==1 and profiles[t]['simplex_dimension']==4
        # The old incidence catalogue predates the forced-base patch; its
        # represented_in_lp flag may be false for a now-present simplex row.
        assert typ['simplex_dimension']==4
        assert 0<=profiles[t]['simplex_face']<82806
        source=context[prov['source_context']]
        assert source[2:4]==[ci,prov['prism_dimension']]
        assert dict(source[4])[t]==1
        assert len([z for z in S['actual_bases']if z['provenance']==pr and z['base_type']==t])==1
        assert raw_source not in groups[t]
        rec=dict(base_index=B['index'],base_type=t,raw_source=raw_source,provenance=pr,face=B['face'],
            source_context=prov['source_context'],canonical_prism_core=ci,deleted=B['deleted'],
            coface=B['coface_root'],core_isomorphism=T['core_isomorphism'],canonical_base_vertices=B['canonical_vertices'])
        groups[t][raw_source]=rec;selected.append(rec)
    assert len(groups)==52 and len(selected)==65
    assert Counter(tuple(sorted(g))for g in groups.values())=={(23,):39,(23,25):13}
    constraints=[];violations=[];updated=dict(mu);targets=[]
    for t,g in sorted(groups.items()):
        typ=I['base_types'][t];Fst=typ['spectator_forest'];k=Fst.count('A1')
        si=profiles[t]['simplex_face'];remaining=list(Fst);remaining.remove('A1');remaining.remove('A1')
        Y=tuple(sorted(remaining+['H3','H3']));v=vlookup[Y];h=Y.count('H3')
        assert h>=2 and sum(map(M.V.rank,Y))==13
        assert v not in targets;targets.append(v)
        positive=[]
        if 23 in g:positive.append([g[23]['face'],k])
        if 25 in g:positive.append([g[25]['face'],1])
        assert len({f for f,c in positive})==len(positive)
        load=sum((c*alpha.get(f,F())for f,c in positive),F())-k*alpha.get(si,F())
        old=h*mu.get(v,F());gap=load-old
        rec=dict(index=len(constraints),base_type=t,simplex_face=si,simplex_dimension=4,
            spectator_forest=Fst,A1_count=k,selected_raw_sources=sorted(g),
            actual_profiles=[g[r]for r in sorted(g)],positive_faces=positive,negative_faces=[[si,k]],
            vertex_capacity=[[v,h]],target_forest=list(Y),forced_quad_spectators=sorted(remaining+['H3']),
            input_load=str(load),input_capacity=str(old),gap=str(gap))
        constraints.append(rec)
        if gap>0:
            violations.append(rec['index']);updated[v]=max(updated.get(v,F()),load/h)
    assert len(violations)==1
    witness=constraints[violations[0]]
    assert witness['base_type']==55293 and witness['simplex_face']==80185
    assert witness['positive_faces']==[[26406,2],[26876,1]] and witness['vertex_capacity']==[[5413,2]]
    assert F(witness['gap'])==F(50071,2500000000)
    assert witness['target_forest']==['A3','H3','H3','I2(5)','I2(5)']
    assert all(F(c['input_load'])<=sum((h*updated.get(v,F())for v,h in c['vertex_capacity']),F())for c in constraints)
    assert {v for v in set(mu)|set(updated)if mu.get(v,F())!=updated.get(v,F())}=={5413}
    assert updated[5413]-mu[5413]==F(50071,5000000000)
    # Recheck the earlier non-family row explicitly under the same alpha.
    assert 2*alpha[26410]+alpha[26880]-2*alpha[580]<=3*updated[5572]
    # Complete count-level proof for each k actually represented here.
    count_tables=[]
    for k in sorted({c['A1_count']for c in constraints}):
        for na in (0,1):
            for nb in range(k-na+1):
                assert na*nb>=k*na+nb-k
                count_tables.append([k,na,nb,na*nb,k*na+nb-k])
    normalizer=sum(updated.values(),F());objective=F(base['objective_numerator']);bound=objective/normalizer
    assert sum(mu.values(),F())==F(base['normalizer'])and objective/F(base['normalizer'])==F(base['bound'])
    assert 12<bound<F(base['bound'])
    paths=[Path(__file__),Path(M.__file__),Path(M.V.__file__)]
    paths += [ROOT/f'results/{tag}.json'for key,tag in M.TAGS.items()if key!='model']
    paths += [ROOT/f'results/{tag}.json'for tag in (BASE,GATE,FORESTS,FULL)]
    paths += [ROOT/f'results/{GATE}.md',ROOT/f'results/{GATE}_review.md']
    return dict(status='PASS_EXACT_RAW23_25_REPRESENTED_SPECTATOR_INCIDENCE_AND_REPAIR',dimension=13,
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        counts=dict(common_simplex_types=len(groups),actual_profile_occurrences=len(selected),paired_types=13,
            only_selected_A_types=39,positive_source_types=sum(any(alpha.get(f,F())for f,c in r['positive_faces'])for r in constraints),
            violations=len(violations),changed_vertex_masses=1),
        literal_profile_occurrences=selected,constraints=constraints,violation_indices=violations,integer_count_tables=count_tables,
        input_bound=base['bound'],bound=str(bound),objective_numerator=str(objective),normalizer=str(normalizer),
        face_multipliers=base['face_multipliers'],vertex_multipliers=[[v,str(a)]for v,a in sorted(updated.items())],
        vertex_updates=[[v,str(mu.get(v,F())),str(updated[v])]for v in sorted(updated)if updated[v]!=mu.get(v,F())],
        face_masses_unchanged=True,alpha_changed=False,optimizer_invocations=0,repeated_large_scans=0,
        inherited_preservation='Every face load is unchanged. The old complete local covers, signed/cap inequalities and17 joint families persist because only mu5413 increases; the old forced-quad row is explicitly rechecked.',
        scope='Complete matching raw23/25 family within the frozen actual-provenance domain,52 represented S types. Missing profile coefficients are zero, giving a weaker valid necessary inequality. The same local geometry applies with orthogonal finite spectators. Not every ambient profile or global incidence is covered; repaired partial-model obstruction only, no optimum or dimension theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    data=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(path.read_text())
    else:path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['counts'],float(F(data['bound'])),data['vertex_updates'],flush=True)
