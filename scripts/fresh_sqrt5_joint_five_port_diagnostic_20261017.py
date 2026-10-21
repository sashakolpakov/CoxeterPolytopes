#!/usr/bin/env python3
"""All fifteen exact five-port facets after the frozen four-port repair."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path
import json

import fresh_sqrt5_joint_triangle_square_ports_20261017 as J

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_joint_five_port_diagnostic_20261017'
THEOREM='fresh_sqrt5_joint_five_port_domain_20261017'
BASE='fresh_sqrt5_joint_four_port_diagnostic_20261017'
C=J.C


def dump(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def zeros():return [F() for _ in range(4)]


def run():
    base=J.read(BASE);theorem=J.read(THEOREM);local=J.read(J.LOCAL)
    assert theorem['status']=='PASS_EXACT_FIVE_PORT_SELECTED_PRISM_FACETS'
    facets=theorem['necessary_domain']['facets']
    assert len(facets)==15
    for i,f in enumerate(facets):
        assert f['index']==i and len(f['coefficients'])==8 and min(f['coefficients'])>=0 and f['capacity']>0
    t=next(x for x in local['tests'] if x['name']=='legacy_export_seed')
    assert t['complete_positive_simplex_forest_coverage'] and not t['open_forests']
    assert base['face_multipliers']==t['face_multipliers']
    alpha={i:F(a) for i,a in base['face_multipliers']};mu={v:F(a) for v,a in base['vertex_multipliers']}
    assert len(alpha)==2219 and len(mu)==593
    model=json.loads((ROOT/f'results/{J.MODEL}.json').read_text())
    for p,h in model['prepared']['sources'].items():assert J.B.P.sha(p)==h,p
    inverse={fi:j for j,fi in model['prepared']['provenance_to_face']}
    selected_q=sorted(set(alpha)&set(inverse));assert len(selected_q)==355
    actual=J.read(J.ACTUAL);quad=J.read(J.B.Q.TAG);raw=J.read(J.B.Q.CORES);spectators=J.read(J.B.Q.SPECTATORS)
    forests=J.read(J.FORESTS)['vertices'];lookup={tuple(v['forest']):v['index'] for v in forests}
    models=dict(C.models());finite={C.code(g):name for name,g in models.items()}
    graphs={};orbits={};code_to_oi={};representatives={}

    def ensure(v):
        if v in graphs:return
        g={}
        for name in forests[v]['forest']:
            offset=len(g)
            for i,row in models[name].items():g[offset+i]={offset+j:m for j,m in row.items()}
        assert set(g)==set(range(13));graphs[v]=g;pieces=C.pieces(g);groups=defaultdict(list)
        for U in combinations(range(13),5):groups[tuple(sorted(C.code(p,U) for p in pieces))].append(U)
        assert sum(map(len,groups.values()))==1287
        rows=[]
        for oi,(code,sets) in enumerate(sorted(groups.items())):
            rows.append(dict(orbit=oi,code=json.loads(dump(code)),representative=list(sets[0]),multiplicity=len(sets)))
            code_to_oi[v,code]=oi;representatives[v,oi]=sum(1<<p for p in sets[0])
        orbits[v]=rows

    blockcols={c['index']:c for c in t['literal_block_orbits']};weights=defaultdict(F);sfeatures=defaultdict(list)
    for f in t['positive_simplex_features']:
        touched=set()
        for j,m in f['block_orbits']:weights[j]+=m*alpha[f['face']];touched.add(blockcols[j]['vertex'])
        for v in touched:sfeatures[v].append(f)
    assert set(weights)==set(blockcols) and len(blockcols)==3485
    scoeff={};sload=defaultdict(zeros);sidentities=[]
    for serial,(j,c) in enumerate(sorted(blockcols.items())):
        v=c['vertex'];ensure(v);masks=c['literal_masks'];k=c['block_size']
        assert masks and len(masks)==len(set(masks)) and all(m.bit_count()==k for m in masks)
        totals=zeros()
        for row in orbits[v]:
            oi=row['orbit'];U=representatives[v,oi];hist=Counter((m&U).bit_count() for m in masks)
            coeff=tuple(F(hist[a],len(masks)) for a in range(2,6))
            if any(coeff):
                scoeff[j,oi]=coeff
                for h,a in enumerate(coeff):
                    sload[v,oi][h]+=weights[j]*a;totals[h]+=a*row['multiplicity']
        expected=[comb(k,a)*comb(13-k,5-a) if a<=k else 0 for a in range(2,6)]
        assert totals==expected;sidentities.append([j,k,list(map(str,totals))])
        if serial%400==0:print('FIVE_PORT_SIMPLEX',serial+1,'/',len(blockcols),'forests',len(graphs),flush=True)

    qload=defaultdict(zeros);qcoeff=defaultdict(lambda:defaultdict(zeros));qrecords=[];qidentities=[]
    oldq={r['face']:r for r in base['actual_quad_features']};assert set(oldq)==set(selected_q)
    for serial,fi in enumerate(selected_q):
        j=inverse[fi];rec=actual['provenance_masses'][j];ci=rec['canonical_prism_core'];core=quad['canonical_cores'][ci]
        g0,_=J.B.Q.graph(raw['records'][core['source_core']]);g={i:{k:m[0] for k,m in row.items()} for i,row in g0.items()}
        sp=spectators['contexts'][rec['source_context']]['spectator_forest']
        for name in sp:
            off=len(g)
            for i,row in models[name].items():g[off+i]={off+k:m for k,m in row.items()}
        assert len(g)==15
        corners=J.X.quad_features(dict(quadrilateral_black_partitions=core['quad_partitions']));assert len(corners)==2*rec['prism_dimension']
        previous=oldq[fi];assert (previous['provenance'],previous['canonical_core'],previous['spectator_forest'])==(j,ci,sp)
        cr=[]
        for corner,old in zip(corners,previous['corners']):
            omitted=set(corner['vertex_omitted']);h={i:{k:m for k,m in row.items() if k not in omitted} for i,row in g.items() if i not in omitted}
            assert len(h)==13 and all(m>0 for row in h.values() for m in row.values())
            pieces=C.pieces(h);forest=tuple(sorted(finite[C.code(p)] for p in pieces));v=lookup[forest];ensure(v)
            edges=set(map(tuple,corner['occupied_pairs']));e=len(edges)
            assert old['omitted']==sorted(omitted) and old['vertex']==v and old['occupied_pairs']==corner['occupied_pairs']
            assert e==rec['prism_dimension']-1 and 1<=e<=4
            hubs=set.intersection(*(set(edge) for edge in edges));assert hubs;hub=min(hubs)
            leaves=set().union(*(set(edge) for edge in edges))-{hub};assert len(leaves)==e
            hits=defaultdict(Counter);totalcounts=Counter()
            # Exactly zero if U omits the hub or contains no leaf. Every
            # positive case is among these C(12,4) subsets containing hub.
            for rest in combinations(sorted(set(h)-{hub}),4):
                count=len(leaves&set(rest))
                if not count:continue
                U=tuple(sorted((hub,)+rest));code=tuple(sorted(C.code(p,U) for p in pieces));oi=code_to_oi[v,code]
                hits[oi][count]+=1;totalcounts[count]+=1
            totals=zeros()
            for oi,counts in hits.items():
                multiplicity=orbits[v][oi]['multiplicity']
                assert sum(counts.values())<=multiplicity
                for z in range(1,5):
                    value=F(counts[z],multiplicity)
                    if value:
                        qload[v,oi][z-1]+=alpha[fi]*value;qcoeff[v,oi][fi][z-1]+=value
                        totals[z-1]+=value*multiplicity
            expected=[comb(e,z)*comb(12-e,4-z) if z<=e else 0 for z in range(1,5)]
            assert totals==expected and [totalcounts[z] for z in range(1,5)]==expected
            qidentities.append([fi,sorted(omitted),e,list(map(str,totals))])
            cr.append(dict(omitted=sorted(omitted),vertex=v,occupied_pairs=corner['occupied_pairs'],
                positive_five_port_subsets=sum(totalcounts.values()),histogram_sum=list(map(str,totals))))
        qrecords.append(dict(face=fi,provenance=j,canonical_core=ci,prism_dimension=rec['prism_dimension'],spectator_forest=sp,corners=cr))
        if serial%25==0:print('FIVE_PORT_ACTUAL_Q',serial+1,'/',len(selected_q),flush=True)

    updates=dict(mu);checks=[];critical=[];facetcounts=Counter();facetforests=defaultdict(set)
    for v,rows in sorted(orbits.items()):
        for row in rows:
            oi=row['orbit'];hist=sload[v,oi]+qload[v,oi]
            checks.append([v,oi,list(map(str,hist))])
            for f in facets:
                co=f['coefficients'];cap=f['capacity'];load=sum((a*b for a,b in zip(co,hist)),F())
                updates[v]=max(updates[v],load/cap)
                if load<=cap*mu[v]:continue
                facetcounts[f['index']]+=1;facetforests[f['index']].add(v)
                terms=Counter()
                for fi,qq in qcoeff[v,oi].items():
                    val=sum((a*b for a,b in zip(co[4:],qq)),F())
                    if val:terms[fi]+=val
                for sf in sfeatures[v]:
                    val=F()
                    for j,m in sf['block_orbits']:
                        if blockcols[j]['vertex']!=v:continue
                        cc=scoeff.get((j,oi))
                        if cc:val+=m*sum((a*b for a,b in zip(co[:4],cc)),F())
                    if val:terms[sf['face']]+=val
                assert sum((alpha[fi]*a for fi,a in terms.items()),F())==load
                critical.append(dict(index=len(critical),vertex=v,five_port_orbit=oi,facet=f['index'],
                    representative=row['representative'],orbit_multiplicity=row['multiplicity'],vertex_cost=cap,
                    positive_face_coefficients=[[fi,str(a)] for fi,a in sorted(terms.items())],histogram_loads=list(map(str,hist)),
                    load=str(load),old_vertex_mass=str(mu[v]),required_vertex_mass=str(load/cap)))
    assert set(graphs)==set(mu) and len(graphs)==593
    for v,oi,values in checks:
        hist=list(map(F,values))
        for f in facets:assert sum((a*b for a,b in zip(f['coefficients'],hist)),F())<=f['capacity']*updates[v]
    assert all(updates[v]>=a for v,a in mu.items())
    numerator=F(base['objective_numerator']);assert numerator==F(t['objective_numerator']);den=sum(updates.values(),F())
    paths=[f'scripts/{TAG}.py',f'scripts/{THEOREM}.py',f'results/{THEOREM}.json',f'results/{THEOREM}.md',
        f'scripts/{BASE}.py',f'results/{BASE}.json',f'scripts/{J.TAG}.py',f'results/{J.LOCAL}.json',
        f'results/{J.MODEL}.json',f'results/{J.ACTUAL}.json',f'results/{J.FORESTS}.json',
        f'scripts/{J.X.TAG}.py',f'results/{J.B.Q.TAG}.json',f'results/{J.B.Q.CORES}.json',f'results/{J.B.Q.SPECTATORS}.json']
    return dict(status='PASS_EXACT_FIFTEEN_FIVE_PORT_FACETS_AFTER_FOUR_PORT_REPAIR',sources={p:J.B.P.sha(p) for p in paths},
        dimension=13,target=12,optimizer_invocations=0,initial_bound=base['bound'],bound=str(numerator/den),
        objective_numerator=str(numerator),normalizer=str(den),facet_coordinate_order=theorem['necessary_domain']['coordinate_order'],
        facets=[dict(index=f['index'],coefficients=f['coefficients'],capacity=f['capacity']) for f in facets],
        counts=dict(simplex_faces=1864,simplex_block_orbits=3485,actual_quad_provenances=len(qrecords),
            quad_vertex_occurrences=sum(len(r['corners']) for r in qrecords),tested_forests=len(graphs),
            literal_five_port_subsets=1287*len(graphs),five_port_orbits=len(checks),facet_orbit_inequalities=15*len(checks),
            violated_columns=len(critical),increased_vertex_masses=sum(updates[v]>mu[v] for v in updates)),
        facet_violation_counts=[dict(facet=f['index'],violated_columns=facetcounts[f['index']],
            violating_forests=len(facetforests[f['index']])) for f in facets],
        finite_forest_five_port_orbits=[dict(vertex=v,forest=forests[v]['forest'],orbits=rows) for v,rows in sorted(orbits.items())],
        actual_quad_features=qrecords,simplex_histogram_sum_identities=sidentities,quad_histogram_sum_identities=qidentities,
        all_orbit_histogram_loads=checks,critical_columns=critical,face_multipliers=t['face_multipliers'],
        initial_vertex_multipliers=base['vertex_multipliers'],vertex_multipliers=[[v,str(a)] for v,a in sorted(updates.items())],
        scope='Same2219face masses and355actual selected dotted-Q provenances, after all frozen simplex-only, flag, strengthened triple and four-port repairs. Exact intersection-size histograms are used, including nonmonotone simplex facet coefficients. All1287five-subsets and all15facets are covered for every contributing forest. Per-vertex maxima raise only mu, preserving all inherited positive and signed constraints. This is a complete scan of these fifteen facets, not a complete joint local cone, global LP, realization or dimension theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=json.loads(dump(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(out.read_text())
    else:out.write_text(dump(result)+'\n')
    print(result['status'],result['counts'],result['bound'],float(F(result['bound'])),flush=True)
