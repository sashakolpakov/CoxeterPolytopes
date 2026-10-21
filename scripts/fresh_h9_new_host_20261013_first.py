#!/usr/bin/env python3
"""Complete all-pairing mirror census from the origin to the new endpoint."""
import hashlib,json,time
from collections import Counter
from fresh_h9_new_host_20261013_edge import *
from fresh_h9_quartic_first_endpoint_20261013 import line_key,finite_orbit,root_test

EDGE=OUT
RELATIVE=ROOT/'results/fresh_h9_new_host_20261013_relative.json'
OUT=ROOT/'results/fresh_h9_new_host_20261013_first.json'

def run():
    started=time.monotonic();data=json.loads(EDGE.read_text());rel=json.loads(RELATIVE.read_text());fd=json.loads(FACE.read_text())
    g=Host(1).g;walls=[tuple(map(K.from_json,r))for r in fd['all_walls']];edge=[walls[i]for i in data['edge_indices']]
    eg=tuple(tuple(pair(r,g,s)for s in edge)for r in edge);ei=inverse(eg)
    start=(O,*(Z,)*5);end=tuple(map(K.from_json,rel['nearest_candidates'][0]['endpoint']));direction=tuple(map(K.from_json,data['direction']))
    H=ALPHA/(-pair(end,g,end))-1;assert H>0 and pair(end,g,end)<0
    span={};folded={};records=[];bounds=[];profiles=[];stats=Counter()
    for q in NORMS:
        alph=[balanced_alphabet((q*pair(e,g,e),)*4)for e in edge]
        qstart=len(profiles)
        for fixed in product(*alph):
            co=mv(ei,tuple(-x for x in fixed));r0=tuple(sum((c*e[j]for c,e in zip(co,edge)),Z)for j in range(6))
            R=q-pair(r0,g,r0)
            if R==0:
                if root_test(r0,q,g):span[line_key(r0)]=r0
                continue
            if any(R.sign(e)<=0 for e in range(4)):continue
            caps=(H*q/ALPHA,*(-R/ALPHA for _ in range(3)))
            positive=balanced_alphabet(caps);times=sorted({x for v in positive for x in(v,-v)})
            idx=len(profiles);profiles.append({'norm':q.to_json(),'fixed_magnitudes':enc(fixed),'fixed_component':enc(r0),
                    'residual':R.to_json(),'time_caps':enc(caps),'time_coordinates':enc(times)})
            for x in times:
                stats['quadratics']+=1;t=square_root((R+ALPHA*x*x)/6)
                if t is None:continue
                stats['field_squares']+=1
                for ss in {t,-t}:
                    r=tuple(a+ss*b for a,b in zip(r0,direction));r=(x,*r[1:]);assert pair(r,g,r)==q
                    if not root_test(r,q,g):continue
                    stats['primitive_integral']+=1;a,b=pair(start,g,r),pair(end,g,r)
                    if a*b>0:continue
                    assert a or b
                    lam=-a/(b-a);assert 0<=lam<=1
                    key=line_key(r)
                    if key in folded:continue
                    folded[key]=r;records.append({'profile':idx,'vector':enc(r),'norm':q.to_json(),'parameter':lam.to_json()})
        bounds.append({'norm':q.to_json(),'pairing_alphabets':[enc(v)for v in alph]})
        print({'q':str(q),'profiles':len(profiles)-qstart,'quadratics':stats['quadratics'],'segment_mirrors':len(records),'seconds':round(time.monotonic()-started,2)},flush=True)
    egroots=finite_orbit(edge,edge,g);spanall=finite_orbit(list(span.values()),edge,g);assert set(spanall)==set(egroots)
    expanded=finite_orbit(list(folded.values()),edge,g)
    initial={key:r for key,r in expanded.items()if pair(start,g,r)==0};initial.update(egroots)
    initial_expected=finite_orbit(walls[:5],walls[:5],g);assert set(initial)==set(initial_expected)
    positive=[r for r in records if K.from_json(r['parameter'])>0]
    nearest=min(K.from_json(r['parameter'])for r in positive);assert nearest==1
    endpoint={key:r for key,r in expanded.items()if pair(end,g,r)==0};endpoint.update(egroots)
    choices=[]
    for item in rel['nearest_candidates']:
        r9=tuple(map(K.from_json,item['vector']));r=(r9[0],*r9[1:5],r9[9])
        group=finite_orbit(edge+[r],edge+[r],g)
        choices.append({'relative_profile':item['profile'],'vector':enc(r),'norm':item['norm'],
                        'generated_mirrors':len(group),'complete_endpoint_stabilizer':set(group)==set(endpoint)})
    complete=[c for c in choices if c['complete_endpoint_stabilizer']];assert complete
    selected=min(complete,key=lambda c:ALPHA*K.from_json(c['vector'][0])**2/K.from_json(c['norm']))
    selected_r=tuple(map(K.from_json,selected['vector']));selected_gram=tuple(tuple(pair(r,g,s)for s in edge+[selected_r])for r in edge+[selected_r])
    lifting=[]
    for q in NORMS:
        bs,_=block_vectors(q)
        assert all(en==0 or en==q for en,v in bs)
        lifting.append({'norm':q.to_json(),'dominant_F4_vectors':[{'norm':en.to_json(),'vector':enc(v)}for en,v in bs]})
    out={'status':'COMPLETE_FULL_H9_FIRST_EDGE_MIRROR','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(EDGE,RELATIVE,FACE)},
         'gram':menc(g),'alpha':ALPHA.to_json(),'edge_indices':data['edge_indices'],'edge_roots':menc(edge),
         'start':enc(start),'end':enc(end),'end_norm':pair(end,g,end).to_json(),'height_cap':H.to_json(),
         'direction':enc(direction),'norm_bounds':bounds,'profiles':profiles,'span_roots':menc(list(span.values())),
         'folded_segment_roots':records,'expanded_segment_mirrors':len(expanded),'edge_mirrors':len(egroots),
         'complete_initial_mirrors':len(initial),'nearest_parameter':nearest.to_json(),'complete_endpoint_mirrors':len(endpoint),
         'endpoint_choices':choices,'selected_choice':selected,'selected_incident_gram':menc(selected_gram),
         'H9_other_F4_census':lifting,'H9_edge_mirrors':len(egroots)+24,'H9_initial_mirrors':len(initial)+24,
         'H9_endpoint_mirrors':len(endpoint)+24,'stats':dict(stats),'seconds':round(time.monotonic()-started,3),
         'scope':'Complete finite interval root census independent of the PARI orbit classification: every necessary primitive norm, all bounded integral pairings, all residual classes, all signed time coordinates and both square roots. Exact finite Weyl unfolding checks full initial and endpoint stabilizers. The complete bounded second-F4 census excludes mixed crossing mirrors and lifts the first-edge claim to H9. This proves one full-reflection-chamber edge, not finite volume or host reflectivity.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print({k:out[k]for k in('status','edge_mirrors','complete_initial_mirrors','complete_endpoint_mirrors','H9_endpoint_mirrors','stats','seconds')},flush=True)

if __name__=='__main__':run()
