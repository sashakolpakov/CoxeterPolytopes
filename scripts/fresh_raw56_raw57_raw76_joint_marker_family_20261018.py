#!/usr/bin/env python3
"""Complete three-family marker envelope and coupled charges.

No old marker or ordinary outcome is changed, and no optimizer is called.
"""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_raw56_raw57_extra_vertex_family_20261018 as T

ROOT=T.ROOT;M=T.M;U=T.U;D=U.D;A=T.A;R=T.R
TAG='fresh_raw56_raw57_raw76_joint_marker_family_20261018'


def relative(old,new,raw):
    p=dict(old['active_embedding']);q=dict(new['active_embedding'])
    if raw==57:return()
    assert raw==56
    names={p[0]:0,p[4]:1}
    for i in(q[0],q[5]):
        if i not in names:names[i]=len(names)
    return(names[q[0]],names[q[5]])


def cross_theorem(old,new):
    tests=[]
    for raw,forest,count in [(56,['A1']*8+['H3','I2(10)'],7),
                             (57,['A1']*6+['A2','H3','I2(10)'],1)]:
        row=next(r for r in old['forest_records']if Counter(r['forest'])==Counter(forest))
        other=next(r for r in new['forest_records']if r['vertex']==row['vertex'])
        target={int(i):{int(j):m for j,m in r.items()}for i,r in other['finite_frame'].items()}
        assert row['finite_frame']==other['finite_frame']
        J,active,components,tens=T.projection_data(target)
        ms=next(f['markers']for f in other['marker_families']if f['raw_family']==raw)
        anchor=row['markers'][0];reps={}
        for m in ms:reps.setdefault(relative(anchor,m,raw),m)
        assert len(reps)==count
        records=[]
        for pattern,m in sorted(reps.items()):
            result=T.pair_test(J,tens,anchor['native'],m['native'])
            assert result['allowed']==(raw==56 and pattern==(2,3))
            records.append(dict(pattern=list(pattern),new_marker=m,result=result))
            print('CROSS_MARKER',raw,pattern,result['allowed'],[r['code']for r in result['expressions']],flush=True)
        tests.append(dict(new_raw=raw,host_vertex=row['vertex'],host_forest=forest,
            finite_frame=other['finite_frame'],raw76_anchor=anchor,records=records))
    return tests


def compatible(ra,a,rb,b):
    if ra>rb:return compatible(rb,b,ra,a)
    if rb!=76:return T.allowed_pattern(ra,rb,T.relative(ra,a,rb,b))
    if ra==76:return U.marker_rule(a,b)
    return ra==56 and relative(b,a,ra)==(2,3)


def envelope(forest,families):
    """Complete count upper boundary, with actual graph-clique witnesses."""
    c=Counter(forest);a,h,t=c['A1'],c['H3'],c['I2(10)']
    cap76=min(a-1,h)if 76 in families else 0
    cap56=(a//2)*t if 56 in families else 0
    cap57=t if 57 in families else 0
    counts=[]
    if cap56:counts.append((cap56,0,0))
    for k in range(1,cap76+1):
        q=((a-k-1)//2)*t if cap56 else 0
        assert q>=0;counts.append((q,0,k))
    if cap57:counts.append((0,cap57,0))
    assert counts and len(counts)==len(set(counts))
    # Every clique lies coordinatewise below one row. Coefficients are
    # nonnegative, so these rows suffice for the entire down-closed domain.
    if cap56 and cap76:
        constraints=[(F(q),F(k),F(1))for q,_,k in counts if q or k]+[(F(1),F(),F()),(F(),F(1),F())]
        vertices=set()
        for (aa,bb,ee),(cc,dd,ff)in combinations(constraints,2):
            det=aa*dd-bb*cc
            if not det:continue
            x,y=(ee*dd-bb*ff)/det,(aa*ff-ee*cc)/det
            if x>=0 and y>=0 and all(q*x+k*y<=1 for q,_,k in counts):vertices.add((x,y))
        maximal=sorted(p for p in vertices if not any(p!=q and all(x<=y for x,y in zip(p,q))for q in vertices))
        weights=[(x,F(1,cap57)if cap57 else F(),y)for x,y in maximal]
    else:weights=[(F(1,cap56)if cap56 else F(),F(1,cap57)if cap57 else F(),F(1,cap76)if cap76 else F())]
    assert all(all(x>=0 for x in w)and all(sum(F(n)*x for n,x in zip(p,w))<=1 for p in counts)for w in weights)
    ones=sorted({dict(m['active_embedding'])[i]for raw in(56,76)if raw in families for m in families[raw]for i in((0,5)if raw==56 else(0,4))})
    if cap56 or cap76:assert len(ones)==a,(forest,sorted(families),ones,a)
    witnesses=[]
    for q,r,k in counts:
        chosen=[]
        if k:
            blocks=sorted({tuple(sorted(dict(m['active_embedding'])[j]for j in(1,3,6)))for m in families[76]})
            assert len(blocks)==h
            for i in range(k):
                m=next(m for m in families[76]if dict(m['active_embedding'])[0]==ones[0]
                    and dict(m['active_embedding'])[4]==ones[i+1]
                    and tuple(sorted(dict(m['active_embedding'])[j]for j in(1,3,6)))==blocks[i])
                chosen.append((76,m['index']))
        if q:
            avail=ones[k+1:]if k else ones
            pairs=set(zip(avail[::2],avail[1::2]))
            picked=[m for m in families[56]if(dict(m['active_embedding'])[0],dict(m['active_embedding'])[5])in pairs
                and dict(m['active_embedding'])[1]<dict(m['active_embedding'])[3]]
            assert len(picked)==q;chosen.extend((56,m['index'])for m in picked)
        if r:
            anchor=dict(families[57][0]['active_embedding'])
            picked=[m for m in families[57]if(dict(m['active_embedding'])[0],dict(m['active_embedding'])[5])==(anchor[0],anchor[5])
                and dict(m['active_embedding'])[1]<dict(m['active_embedding'])[3]]
            assert len(picked)==r;chosen.extend((57,m['index'])for m in picked)
        assert Counter(raw for raw,i in chosen)==Counter({raw:n for raw,n in zip((56,57,76),(q,r,k))if n})
        assert all(compatible(ra,families[ra][i],rb,families[rb][j])for(ra,i),(rb,j)in combinations(chosen,2))
        witnesses.append(dict(counts=[q,r,k],marker_clique=chosen))
    return dict(role_counts=dict(A1=a,H3=h,I210=t),raw_order=[56,57,76],
        boundary_counts=counts,coefficient_vertices=weights,attaining_graph_cliques=witnesses)


def combined(old,new,model):
    older={r['vertex']:r for r in old['forest_records']};newer={r['vertex']:r for r in new['forest_records']}
    assert len(older)==235 and len(newer)==1212 and len(older.keys()&newer.keys())==53
    dual=next(c for c in model['certificates']if c['kind']=='EXACT_PARTIAL_MODEL_DUAL')
    alpha={f:F(x)for f,x in dual['face_multipliers']};mu={v:F(x)for v,x in dual['vertex_multipliers']}
    baseline={v:F(x)for v,x in new['exported']['repair']['vertex_multipliers']}
    assert new['exported']['repair']['face_multipliers']==dual['face_multipliers']
    assert all(baseline.get(v,F())>=x for v,x in mu.items())
    candidates=defaultdict(list);vectors={};domains,represented=D.L.current_domains()
    vectors.update({v:d['domain']['maximal_count_vectors']for v,d in domains.items()})
    for tag in(D.L.TAG,D.B.TAG,D.THIRD):
        for r in M.read(tag)['forest_records']:
            candidates[r['vertex']].append((tag,r['certificate']['weights']))
            if tag!=D.L.TAG:vectors[r['vertex']]=r['maximal_count_vectors']
    columns=[];records=[];newmu=dict(baseline);overlap_pairs=retained_pairs=0
    for v in sorted(older.keys()|newer.keys()):
        ro=older.get(v);rn=newer.get(v);r=rn or ro
        if ro and rn:
            for key in('forest','finite_frame','old_nodes','old_groups','old_graph_edges','prior_history'):assert ro[key]==rn[key],(v,key)
        groups=r['old_groups'];owner={i:z['index']for z in groups for i in z['literal_stars']}
        fam={z['raw_family']:z['markers']for z in rn['marker_families']}if rn else{}
        faces={z['raw_family']:z['face']for z in rn['marker_families']}if rn else{}
        edges=[(raw,i,k)for z in rn['marker_families']for i,k in z['marker_ordinary_edges']for raw in[z['raw_family']]]if rn else[]
        if ro:
            fam[76]=ro['markers'];faces[76]=ro['face'];edges.extend((76,i,k)for i,k in ro['marker_ordinary_edges'])
        # Whole native covectors distinguish ordinary stars and all three
        # marker types. No identification or multiplicity is inferred from
        # coincident diagram names alone.
        keys=[D.H.frozen(m['native']['contacts'])for ms in fam.values()for m in ms]
        oldkeys={D.H.frozen(s['contacts'])for s in r['old_nodes']}
        assert len(keys)==len(set(keys))and not oldkeys&set(keys)
        env=envelope(r['forest'],fam);omitted=sorted({owner[i]for raw,i,k in edges})
        crosscount=retained=0
        if ro and rn:
            for raw,ms in fam.items():
                if raw==76:continue
                for x in fam[76]:
                    for y in ms:
                        crosscount+=1;retained+=compatible(76,x,raw,y)
        overlap_pairs+=crosscount;retained_pairs+=retained
        vv=vectors[v]if groups else[[]];cc=candidates[v]if groups else[('marker_only',[])];assert cc
        unique={}
        for tag,weights in cc:
            original=list(map(F,weights))+[F()]*(len(groups)-len(weights))
            assert len(original)==len(groups)and all(x>=0 for x in original)
            assert all(len(row)==len(original)and sum(w*n for w,n in zip(original,row))<=1 for row in vv)
            w=[F()if i in omitted else z for i,z in enumerate(original)]
            assert all(w[owner[i]]==0 for raw,i,k in edges)
            for mw in env['coefficient_vertices']:
                cf=Counter()
                for weight,z in zip(w,groups):
                    for face,n in z['positive_face_coefficients']:cf[face]+=weight*n
                for raw,weight in zip((56,57,76),mw):
                    if weight:assert raw in faces;cf[faces[raw]]+=weight
                key=tuple((f,x)for f,x in sorted(cf.items())if x)
                load=sum((alpha.get(f,F())*x for f,x in key),F())
                rec=dict(vertex=v,source_charge=tag,other_charge_sources=[],old_orbit_weights=list(map(str,original)),
                    omitted_orbits=omitted,retained_orbit_weights=list(map(str,w)),
                    marker_terms=[[raw,faces[raw],str(weight)]for raw,weight in zip((56,57,76),mw)if weight],
                    positive_face_coefficients=[[f,str(x)]for f,x in key],vertex_cost='1',load=str(load),
                    old_model_cost=str(mu.get(v,F())),baseline_cost=str(baseline.get(v,F())),
                    violation=str(load-baseline.get(v,F())),old_vector_count=len(vv),old_vectors_sha256=D.L.jhash(vv))
                if key in unique:unique[key]['other_charge_sources'].append(tag)
                else:unique[key]=rec
        cols=list(unique.values());columns.extend(cols)
        newmu[v]=max([baseline.get(v,F())]+[F(c['load'])for c in cols])
        records.append(dict(vertex=v,forest=r['forest'],raw_families=sorted(fam),source_faces=sorted(faces.items()),
            source_records=[[U.TAG if raw==76 else T.TAG,v,raw]for raw in sorted(fam)],
            old_nodes_sha256=D.L.jhash(r['old_nodes']),old_groups_sha256=D.L.jhash(groups),old_graph_sha256=D.L.jhash(r['old_graph_edges']),
            ordinary_omitted_orbits=omitted,marker_ordinary_unknowns=len(edges),cross_marker_positions=crosscount,
            retained_cross_marker_pairs=retained,envelope=env,baseline_cost=str(baseline.get(v,F())),
            new_cost=str(newmu[v]),increase=str(newmu[v]-baseline.get(v,F()))))
    preserved=[]
    for tag,d in((U.TAG,old),(T.TAG,new)):
        for c in d['exported']['columns']:
            load=sum((alpha.get(f,F())*F(x)for f,x in c['positive_face_coefficients']),F());v=c['vertex']
            assert c['vertex_cost']=='1'and load<=baseline.get(v,F())<=newmu.get(v,F())
            preserved.append(dict(source=tag,vertex=v,load=str(load),cost=str(newmu.get(v,F()))))
    multi=[]
    for c in model['prepared']['inherited_multi_vertex_columns']:
        load=sum((alpha.get(f,F())*n for f,n in c['integer_face_coefficients']),F())
        oldcost=sum((mu.get(v,F())*n for v,n in c['integer_vertex_capacities']),F())
        cost=sum((newmu.get(v,F())*n for v,n in c['integer_vertex_capacities']),F());assert load<=oldcost<=cost
        multi.append(dict(variable=c['variable'],all_vertex_costs=c['integer_vertex_capacities'],load=str(load),old_cost=str(oldcost),new_cost=str(cost)))
    assert len(records)==1394 and len(multi)==14 and len(preserved)==1454
    assert all(newmu.get(v,F())>=x for v,x in baseline.items())
    num=F(dual['objective_numerator']);den=sum(newmu.values(),F())
    assert num==F(new['exported']['repair']['objective_numerator'])
    return dict(forest_records=records,exported=dict(columns=columns,preserved_columns=preserved,preserved_multi_owner_costs=multi,
        repair=dict(face_multipliers=dual['face_multipliers'],vertex_multipliers=[[v,str(x)]for v,x in sorted(newmu.items())if x],
            objective_numerator=str(num),normalizer=str(den),bound=str(num/den),baseline_bound=new['exported']['repair']['bound'],
            scope='Least monotone fixed-alpha repair for prior plus the exported three-family columns only; not a full enlarged-cone optimum or global optimization.')),
        counts=dict(targets=len(records),overlap_targets=53,family_classes=1643,marker_orientations=24788,
            cross_marker_positions=overlap_pairs,retained_cross_marker_pairs=retained_pairs,
            exported_columns=len(columns),strictly_violated_columns=sum(F(c['violation'])>0 for c in columns),
            raised_allowances=sum(F(r['increase'])>0 for r in records)))


def run():
    old=M.read(U.TAG);new=M.read(T.TAG);model=M.read(T.MODEL)
    assert old['status']=='PASS_RAW76_COMPLETE_FAMILY_COUPLED_CHARGES'
    assert new['status']=='PASS_RAW56_RAW57_COMPLETE_EXTRA_VERTEX_FAMILY'
    for data in(old,new,model['prepared']):
        for p,h in data['sources'].items():assert M.digest(p)==h,p
    assert len(old['family'])+len(new['family'])==1643
    assert sum(len(r['markers'])for r in old['forest_records'])+sum(len(f['markers'])for r in new['forest_records']for f in r['marker_families'])==24788
    tests=cross_theorem(old,new);out=combined(old,new,model)
    support=M.read('fresh_raw76_complete_family_support_20261018')
    assert support['objective_numerator']==out['exported']['repair']['objective_numerator']
    remaining=[r for r in support['ranked_unambiguous_missing_families']if r['raw_family']not in(56,57)and r['missing_classes']>0]
    assert sum(r['missing_classes']for r in remaining)==1213 and len(remaining)==211
    repair=out['exported']['repair'];gap=F(repair['objective_numerator'])/12-F(repair['normalizer']);assert gap>0
    out['remaining_support_diagnostic']=dict(unrepresented_unambiguous_classes=1213,
        unrepresented_unambiguous_families=211,highest_saved_demand_families=remaining[:10],
        added_normalizer_needed_to_reach12=str(gap),
        scope='Saved fixed-mass provenance inventory after recognizing both new marker families. Demand ranks are selection cues, not new marker incidences or capacity promises; merged classes remain separate.')
    paths={str(Path(m.__file__).relative_to(ROOT))for m in(T,U,D,A,R,D.L)}|{str(Path(__file__).relative_to(ROOT))}
    for tag in(T.TAG,U.TAG,T.MODEL,D.L.TAG,D.B.TAG,D.THIRD):paths.update(f'{f}/{tag}.{e}'for f,e in(('scripts','py'),('results','json'),('results','md')))
    paths.update(f'{f}/fresh_raw76_complete_family_support_20261018.{e}'for f,e in(('scripts','py'),('results','json')))
    return T.normal(M.enc(dict(status='PASS_RAW56_RAW57_RAW76_JOINT_COMPLETE_FAMILY',
        sources={p:M.digest(p)for p in sorted(paths)},relative_tests=tests,
        relative_cases=sum(len(t['records'])for t in tests),
        retained_patterns=sum(r['result']['allowed']for t in tests for r in t['records']),**out,
        optimizer_invocations=0,current_graphs_changed=False,
        scope='Complete1643-class three-family marker count envelopes and mass-independent one-vertex charges. All old ordinary and marker unknowns/columns preserved. Eight new relative physical tests; prior geometry and ordinary outcomes inherited. No global solve.')))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(path.read_text())==out
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['relative_cases'],out['retained_patterns'],out['counts'],out['exported']['repair']['bound'],flush=True)
