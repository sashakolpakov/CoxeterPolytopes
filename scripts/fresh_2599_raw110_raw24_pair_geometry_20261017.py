#!/usr/bin/env python3
"""Eight genuinely unowned raw110/raw24 pairs after the actual118 successor.

The tested forest is unchanged by the preceding eight-forest patch.
All prior exclusions and unknowns are inherited. No global solve.
"""
import argparse,json,time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import fresh_290_class_shortlist_pair_geometry_20261017 as T

N=T.N;U=T.U;W=T.W;ROOT=T.ROOT
TAG='fresh_2599_raw110_raw24_pair_geometry_20261017'
TARGET=(2599,6,9,8)


def run(verify=False):
    start=time.monotonic();p=N.prepare();base=W.read(T.TAG);domain=W.read(N.TAG);priority=W.read(U.TAG)
    v,oi,oj,n=TARGET
    assert v not in {r['vertex']for r in base['forest_records']}
    assert base['baseline_source']==N.TAG
    for key in('face_multipliers','restricted_face_multipliers'):assert base[key]==domain[key]
    for old,new in zip(domain['repairs'],base['repairs']):
        om={i:F(x)for i,x in old['vertex_multipliers']};nm={i:F(x)for i,x in new['vertex_multipliers']}
        assert om.get(v,F())==nm.get(v,F())and all(nm.get(i,F())>=x for i,x in om.items())
    chosen=next(r for r in priority['ranked_candidate_domains']if r['vertex']==v and r['corner_orbits']==[oi,oj])
    assert chosen['pairs']==n and chosen['raw_families']==[[110],[24]]
    assert chosen['face_coefficients']==[[[48279,1]],[[26695,1]]]
    local_priority=dict(priority,shortlist=dict(chosen=[chosen]))
    sources={}
    for tag in(T.TAG,N.TAG,U.TAG):
        doc=W.read(tag)
        for path,h in doc['sources'].items():assert W.digest(path)==h,(tag,path)
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=W.digest(path)
    for mod in(T,N,U,T.B,T.R,T.H,W):
        path=str(Path(mod.__file__).relative_to(ROOT));sources[path]=W.digest(path)
    sources[f'scripts/{TAG}.py']=W.digest(f'scripts/{TAG}.py')
    saved=W.read(TAG)if verify else None
    if saved:assert saved['sources']==sources
    changed={r['vertex']:r for r in domain['forest_records']}
    row=T.analyze(v,oi,oj,n,p,changed,local_priority,saved['forest_records'][0]if saved else None)
    mus=[{i:F(x)for i,x in r['vertex_multipliers']}for r in base['repairs']];newmus=[dict(m)for m in mus]
    alphas=[p['alpha'],p['ralpha']];columns=[];seen=set()
    for k,c in enumerate(row['certificates']):
        newmus[k][v]=max(newmus[k].get(v,F()),F(c['optimal_mass']));co=Counter()
        for w,g in zip(map(F,c['weights']),row['corner_orbits']):
            for i,x in g['positive_face_coefficients']:co[i]+=w*F(x)
        terms=sorted((i,x)for i,x in co.items()if x);loads=[sum((a.get(i,F())*x for i,x in terms),F())for a in alphas]
        assert loads[k]==F(c['optimal_mass'])and all(loads[j]<=F(row['certificates'][j]['optimal_mass'])for j in range(2))
        key=tuple(terms)
        if key not in seen and any(loads[j]>mus[j].get(v,F())for j in range(2)):
            seen.add(key);columns.append(dict(index=len(columns),vertex=v,vertex_cost='1',orbit_weights=c['weights'],
                positive_face_coefficients=[[i,str(x)]for i,x in terms],loads=list(map(str,loads)),charge_domain=base['repairs'][k]['domain']))
    repairs=[]
    for k,old in enumerate(base['repairs']):
        assert all(newmus[k].get(i,F())>=x for i,x in mus[k].items())
        num=F(old['objective_numerator']);den=sum(newmus[k].values(),F())
        repairs.append(dict(domain=old['domain'],baseline_bound=old['bound'],objective_numerator=str(num),normalizer=str(den),bound=str(num/den),
            vertex_updates=[[i,str(mus[k].get(i,F())),str(x)]for i,x in sorted(newmus[k].items())if x>mus[k].get(i,F())],
            vertex_multipliers=[[i,str(x)]for i,x in sorted(newmus[k].items())if x]))
    out=T.normal(dict(status='PASS_2599_RAW110_RAW24_ACTUAL_PAIR_GEOMETRY',sources=sources,baseline_source=T.TAG,selection_source=U.TAG,
        selected_record=chosen,forest_records=[row],columns=columns,repairs=repairs,counts=row['counts'],
        face_multipliers=base['face_multipliers'],restricted_face_multipliers=base['restricted_face_multipliers'],
        changed_forests=[v]if row['counts']['new_exclusions']else[],inherited_forests=571-bool(row['counts']['new_exclusions']),
        preservation='The target forest is disjoint from every forest changed by the preceding118-pair successor. All earlier owners and explicit unknowns, including1763, are inherited unchanged. Fixed face masses and monotone mu preserve439410 model columns and prior capacities.',
        scope='One complete eight-pair unowned raw110/raw24 native domain, exact physical gates and actual local covering certificates. No hypothetical exclusions assumed, global optimizer, full physical realization, implicit coverage completion or dimension theorem. Fresh replay has no optimizer.'))
    path=ROOT/f'results/{TAG}.json'
    if saved:assert out==saved
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'columns',len(columns),'bounds',[float(F(x['bound']))for x in repairs],'seconds',round(time.monotonic()-start,2),flush=True)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');z=a.parse_args();run(z.verify)
