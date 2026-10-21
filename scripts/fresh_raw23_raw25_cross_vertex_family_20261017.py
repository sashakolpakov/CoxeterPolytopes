#!/usr/bin/env python3
"""Transport the actual nine-root incidence and separate its weighted family.

Generation uses tiny local LPs to propose rational certificates. Verification
reconstructs the native maps and complete clique subsets, with no optimizer.
The stronger zero-event model is explicitly hypothetical, never a pair edit.
"""
import argparse,json,time
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import networkx as nx
import fresh_forest1325_cross_vertex_incidence_20261017 as P

A,R,N,W,V,K,M,ROOT=P.A,P.R,P.N,P.W,P.V,P.K,P.M,P.ROOT
TAG='fresh_raw23_raw25_cross_vertex_family_20261017'

def normal(x):return json.loads(json.dumps(x,sort_keys=True))
def kout(x):return x.out()
def kread(x):return K(*map(F,x))


def check_coupled(items,demand,quota,c):
    w=list(map(F,c['weights']));z=F(c['event_penalty'])
    assert len(w)==len(demand)and z>=0 and all(x>=0 for x in w)
    assert all(sum(a*b for a,b in zip(w,q['vector']))<=1+z*q['events']for q in items)
    mix=[(j,F(x))for j,x in c['mixture']]
    assert len({j for j,x in mix})==len(mix)and all(0<=j<len(items)and x>0 for j,x in mix)
    covered=[sum((x*items[j]['vector'][i]for j,x in mix),F())for i in range(len(demand))]
    assert all(a>=b for a,b in zip(covered,demand))
    events=sum((x*items[j]['events']for j,x in mix),F());assert events<=quota
    cost=sum((x for j,x in mix),F());load=sum((a*b for a,b in zip(w,demand)),F())-z*quota
    assert cost==load==F(c['optimal_mass'])and events==F(c['event_mass'])
    return cost


def coupled(items,demand,quota,saved=None):
    if saved is not None:check_coupled(items,demand,quota,saved);return saved
    import numpy as np
    from scipy.optimize import linprog
    scale=max(demand);mat=np.array([q['vector']for q in items],float).T
    ev=np.array([q['events']for q in items],float)
    b=np.array([float(x/scale)for x in demand])
    result=linprog(np.ones(len(items)),A_ub=np.vstack([-mat,ev]),b_ub=np.r_[-b,float(quota/scale)],
        bounds=(0,None),method='highs',options={'time_limit':20,'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
    assert result.success,result.message
    w=[F(float(x)).limit_denominator(1000000)for x in -result.ineqlin.marginals[:-1]]
    z=F(float(-result.ineqlin.marginals[-1])).limit_denominator(1000000)
    ids=[j for j,x in enumerate(result.x)if x>1e-10]
    covered=mat@result.x
    tight=[i for i,x in enumerate(covered)if abs(x-b[i])<1e-8*max(1e-4,abs(b[i]))]
    rows=[[items[j]['vector'][i]for j in ids]for i in tight];rhs=[demand[i]for i in tight]
    if z or abs(float(ev@result.x)-float(quota/scale))<1e-9:
        rows.append([items[j]['events']for j in ids]);rhs.append(quota)
    xs=W.solve_exact(rows,rhs,len(ids))
    c=dict(weights=list(map(str,w)),event_penalty=str(z),mixture=[[j,str(x)]for j,x in zip(ids,xs)if x],
        event_mass=str(sum((x*items[j]['events']for j,x in zip(ids,xs)),F())),optimal_mass=str(sum(xs,F())))
    check_coupled(items,demand,quota,c);return c


def native_map(a,b,nodes,target,J,reference):
    """Identify actual walls A,B,C afresh in one common complete vertex frame."""
    ca=N.native(nodes[a]);cb=N.native(nodes[b]);shared=set(ca)&set(cb);assert len(shared)==1
    C=next(iter(shared));AA=next(q for q in ca if q!=C);BB=next(q for q in cb if q!=C)
    ta=K(F(7,2),F(1,2));tb=K(6,2)
    assert len(AA)==2 and len(BB)==1 and len(C)==2
    qB=BB[0][0];assert V.squared(BB[0][1])==tb
    qA=next(i for i,e in AA if V.squared(e)==ta)
    w=next(i for i,e in AA if V.squared(e)==1);assert dict(C)[w][0]==3
    h3=next(i for i,e in C if i!=w);assert dict(C)[h3][0]==3
    assert len({qA,qB,w})==3 and all(not target[i]for i in(qA,qB,w))
    assert len(target[h3])==1
    center=next(iter(target[h3]));assert target[h3][center]==3 and len(target[center])==2
    h5=next(i for i in target[center]if i!=h3)
    assert target[center][h5]==5 and len(target[h5])==1
    active=[qB,qA,w,h5,center,h3];assert len(set(active))==6
    passive=sorted(set(target)-set(active));assert len(passive)==7
    assert all(j in active for i in active for j in target[i])
    ordered=[BB,AA,C];heights=[];selfchecks=[]
    for q in ordered:
        field=N.pair_field((q,));p=R.projection(q,q,J,field);assert not any(p[1:])
        h=p[0]-2;assert h.sign()>0;heights.append(h)
    for index,q1,q2 in((a,AA,C),(b,BB,C)):
        field=N.pair_field((q1,q2),(tuple(nodes[index]['native_omitted_product']),))
        pp=R.projection(q1,q2,J,field);native=field.negative_contact(tuple(nodes[index]['native_omitted_product']))
        hs=[heights[ordered.index(q)]for q in(q1,q2)]
        assert field.signs_difference(tuple(x-y for x,y in zip(pp,native)),hs[0]*hs[1])==0
        assert sorted(tuple(q)for q in nodes[index]['squared_timelike_heights'])==sorted(tuple(q.out())for q in hs)
        selfchecks.append(dict(star=index,projection=R.qout(pp),native=R.qout(native)))
    expected={(0,1):-K(F(1,2),F(1,2)),(0,2):K(-1),(1,2):K()};cross=[]
    for i,j in combinations(range(3),2):
        field=N.pair_field((ordered[i],ordered[j]));pp=R.projection(ordered[i],ordered[j],J,field)
        e=(expected[i,j],K(),K(),K())
        assert field.signs_difference(tuple(x-y for x,y in zip(pp,e)),heights[i]*heights[j])==0
        cross.append(dict(walls=[i,j],projection=R.qout(pp),height_product=(heights[i]*heights[j]).out(),product=expected[i,j].out()))
    # Rebuild the identical scaled Gram, rather than assuming a diagram label.
    scale=[K(1),ta]+[K(1)]*7;G=[[2*scale[i]if i==j else K()for j in range(9)]for i in range(9)]
    squares={}
    for i,j in combinations(range(6),2):
        if active[j]in target[active[i]]:squares[i,j]=V.squared(V.edge(target[active[i]][active[j]]))
    for j,q in enumerate(ordered,6):
        for i,e in q:squares[active.index(i),j]=V.squared(e)
    for(i,j),e in expected.items():squares[i+6,j+6]=e*e
    for(i,j),s in squares.items():
        edge=M.ksqrt(s*scale[i]*scale[j])if s else K();assert edge is not None
        G[i][j]=G[j][i]=-edge
    assert P.out(G)==reference['active_gram']
    assert set(target)-set(sum(nodes[b]['occupied_pairs'],[]))==set(passive+[qA,w])
    assert set(target)-set(sum(nodes[a]['occupied_pairs'],[]))==set(passive+[qB])
    return dict(stars=[a,b],active_base_roots=active,passive_roots=passive,
        external_contacts=normal(ordered),squared_timelike_heights=[h.out()for h in heights],self_checks=selfchecks,cross_checks=cross,
        triangle=dict(white_base=sorted(passive+[qB,qA,h5,center]),source_ports=[w,h3],opposite_wall='C'),
        target_support=dict(base=sorted(passive+[w,h5,center,h3]),external=['B','A']),
        target_H3_roles=dict(center='A',leaf5='B',leaf3=w))


def enumerate_items(g,owner,ng,pairs):
    found={};cliques=sorted(sorted(c)for c in nx.find_cliques(g))
    for c in cliques:
        for k in range(1,len(c)+1):
            for sub in combinations(c,k):
                vec=tuple(Counter(owner[i]for i in sub).get(j,0)for j in range(ng))
                x=sum(a in sub and b in sub for a,b in pairs);key=vec,x
                if key not in found or sub<found[key]:found[key]=sub
    return cliques,[dict(vector=v,events=x,clique=found[v,x])for v,x in sorted(found)]


def repair(base,alpha,columns,mode):
    mu={v:F(x)for v,x in base['vertex_multipliers']};new=dict(mu);loads=[]
    for c in columns:
        load=sum((F(q)*alpha.get(f,F())for f,q in c['positive_face_coefficients']),F())
        v,y=c['source_vertex'],c['target_vertex'];z=F(c['target_vertex_cost'])
        new[v]=max(new.get(v,F()),load-z*mu.get(y,F()))
        loads.append(str(load))
    assert all(new[v]+F(c['target_vertex_cost'])*new.get(c['target_vertex'],F())>=F(load)
        for c,load in zip(columns,loads)for v in[c['source_vertex']])
    num=F(base['objective_numerator']);den=sum(new.values(),F())
    return dict(domain=base['domain'],mode=mode,objective_numerator=str(num),normalizer=str(den),bound=str(num/den),loads=loads,
        vertex_updates=[[v,str(mu.get(v,F())),str(x)]for v,x in sorted(new.items())if x>mu.get(v,F())],
        vertex_multipliers=[[v,str(x)]for v,x in sorted(new.items())if x])


def run(verify=False):
    started=time.monotonic();p,changed,oldrows,patches,old=A.prepare();base=W.read(P.TAG)
    saved=W.read(TAG)if verify else None;sr={r['vertex']:r for r in saved['forest_records']}if saved else {}
    alphas=[{i:F(x)for i,x in base[k]}for k in('face_multipliers','restricted_face_multipliers')]
    mus=[{i:F(x)for i,x in r['vertex_multipliers']}for r in base['repairs']]
    vlookup={tuple(v['forest']):v['index']for v in W.read(R.P.FORESTS)['vertices']}
    models=dict(R.C.models());names={R.C.code(q):name for name,q in models.items()}
    records=[];zero=[];columns=[];hypo=[];counts=Counter();proposal_keys=set()
    selected=[r for r in W.read(P.INC)['constraints']if len(r['actual_profiles'])==2];assert len(selected)==13
    for inc in selected:
        af,bf=[q['face']for q in inc['actual_profiles']];v=vlookup[tuple(sorted(inc['spectator_forest']+['A1','H3']))]
        y,h=inc['vertex_capacity'][0]
        if not all(alphas[1].get(f,F())for f in(af,bf)):
            assert all(not all(alpha.get(f,F())for f in(af,bf))for alpha in alphas)
            # Zero-demand selected faces may be omitted from every local cover.
            # Verify this on both existing exact certificates, with all other
            # orbit demands retained. This is stronger than a mass-only claim.
            old,g,owner,unowned,dd,groups,forest,applied=A.current(v,p,changed,oldrows,patches)
            projected=[]
            for di,alpha in enumerate(alphas):
                demand=[sum((alpha.get(f,F())*m for f,m in q['positive_face_coefficients']),F())for q in groups]
                c=old['certificates'][di];W.check_certificate(old['maximal_count_vectors'],demand,c)
                terms=[]
                for vi,amount in c['mixture']:
                    cl=[i for i in old['literal_representatives'][vi]if demand[owner[i]]]
                    assert not(any(groups[owner[i]]['face']==af for i in cl)and any(groups[owner[i]]['face']==bf for i in cl))
                    terms.append(dict(stars=cl,amount=amount,vector=[Counter(owner[i]for i in cl).get(j,0)for j in range(len(groups))]))
                assert all(sum((F(t['amount'])*t['vector'][i]for t in terms),F())>=demand[i]for i in range(len(demand)))
                assert sum((F(t['amount'])for t in terms),F())==F(c['optimal_mass'])<=mus[di].get(v,F())
                projected.append(dict(demands=list(map(str,demand)),terms=terms,cover_mass=c['optimal_mass']))
            zero.append(dict(incidence_row=inc['index'],vertex=v,target_vertex=y,faces=[af,bf],corner_orbits=groups,projected_covers=projected));continue
        old,g,owner,unowned,dd,groups,forest,applied=A.current(v,p,changed,oldrows,patches)
        nodes=A.T.current(v,p,changed)[5];target=R.P.forest_graph(forest,models);J,active,components=R.projection_data(target)
        aa=[i for i,s in enumerate(nodes)if s['origin']['face']==af and s['origin']['omitted']==[2,5]]
        bb=[i for i,s in enumerate(nodes)if s['origin']['face']==bf and s['origin']['omitted']==[2,3]]
        assert aa and bb
        pairs=sorted((a,b)for a in aa for b in bb if g.has_edge(a,b));maps=[]
        for a,b in pairs:
            q=native_map(a,b,nodes,target,J,base)
            ps=set(q['passive_roots']);pg={i:{j:m for j,m in target[i].items()if j in ps}for i in ps}
            passive=sorted(names[R.C.code(z)]for z in R.C.pieces(pg))
            target_forest=sorted(passive+['H3','H3']);assert vlookup[tuple(target_forest)]==y and target_forest.count('H3')==h
            q['target_forest']=target_forest;maps.append(q)
        assert pairs
        cliques,items=enumerate_items(g,owner,len(groups),pairs)
        demands=[[sum((alpha.get(f,F())*m for f,m in q['positive_face_coefficients']),F())for q in groups]for alpha in alphas]
        certs=[];hc=[];cache={};zitems=[q for q in items if q['events']==0];zvectors=[q['vector']for q in zitems]
        for di,demand in enumerate(demands):
            quota=h*mus[di].get(y,F());key=tuple(demand),quota
            if saved:c=coupled(items,demand,quota,sr[v]['certificates'][di])
            elif key in cache:c=cache[key]
            else:c=coupled(items,demand,quota);cache[key]=c;proposal_keys.add((v,'coupled',key))
            certs.append(c)
            if saved:zc=sr[v]['zero_event_certificates'][di];W.check_certificate(zvectors,demand,zc)
            elif ('zero',tuple(demand))in cache:zc=cache['zero',tuple(demand)]
            else:zc=W.exact_lp(zvectors,demand);cache['zero',tuple(demand)]=zc;proposal_keys.add((v,'zero',tuple(demand)))
            W.check_certificate(zvectors,demand,zc);hc.append(zc)
            coefficients=Counter()
            for weight,q in zip(map(F,c['weights']),groups):
                for f,m in q['positive_face_coefficients']:coefficients[f]+=weight*m
            ccol=dict(source_vertex=v,target_vertex=y,target_vertex_cost=str(F(c['event_penalty'])*h),
                positive_face_coefficients=[[f,str(x)]for f,x in sorted(coefficients.items())if x])
            if ccol not in columns:columns.append(ccol)
            coefficients=Counter()
            for weight,q in zip(map(F,zc['weights']),groups):
                for f,m in q['positive_face_coefficients']:coefficients[f]+=weight*m
            hcol=dict(source_vertex=v,target_vertex=y,target_vertex_cost='0',positive_face_coefficients=[[f,str(x)]for f,x in sorted(coefficients.items())if x])
            if hcol not in hypo:hypo.append(hcol)
        rec=dict(incidence_row=inc['index'],vertex=v,forest=forest,target_vertex=y,pointed_capacity=h,faces=[af,bf],
            corner_orbits=groups,literal_stars=nodes,graph_edges=sorted(tuple(sorted(e))for e in g.edges()),compatible_selected_pairs=pairs,
            native_maps=maps,maximal_cliques=cliques,items=items,demands=[list(map(str,d))for d in demands],
            quotas=[str(h*mu.get(y,F()))for mu in mus],certificates=certs,zero_event_items=zitems,zero_event_certificates=hc)
        records.append(rec);counts.update(positive_contexts=1,compatible_native_pairs=len(pairs),maximal_cliques=len(cliques),event_count_vectors=len(items))
        print('CONTEXT',v,'pairs',len(pairs),'vectors',len(items),'actual',[c['optimal_mass']for c in certs],
            'zero_event',[c['optimal_mass']for c in hc],flush=True)
    assert len(records)==7 and len(zero)==6 and counts['compatible_native_pairs']==78
    assert len({r['target_vertex']for r in records})==7
    assert not({r['vertex']for r in records}&{r['target_vertex']for r in records})
    repairs=[repair(r,a,columns,'actual_family')for r,a in zip(base['repairs'],alphas)]
    ceilings=[repair(r,a,hypo,'hypothetical_all_selected_events_zero')for r,a in zip(base['repairs'],alphas)]
    for di in range(2):
        actual_mu={v:F(x)for v,x in repairs[di]['vertex_multipliers']}
        hypothetical_mu={v:F(x)for v,x in ceilings[di]['vertex_multipliers']}
        for r in records:
            v=r['vertex']
            assert actual_mu[v]==max(mus[di].get(v,F()),F(r['certificates'][di]['optimal_mass']))
            assert hypothetical_mu[v]==max(mus[di].get(v,F()),F(r['zero_event_certificates'][di]['optimal_mass']))
    # Every actual coupled cut is weaker than this zero-event comparison.
    for di,(c,alpha)in enumerate(zip(ceilings,alphas)):
        mu={v:F(x)for v,x in c['vertex_multipliers']}
        for col in columns:
            assert sum((F(x)*alpha.get(f,F())for f,x in col['positive_face_coefficients']),F())<=mu[col['source_vertex']]+F(col['target_vertex_cost'])*mu.get(col['target_vertex'],F())
    assert F(ceilings[1]['bound'])>12 and F(repairs[1]['bound'])>12
    sources=A.bindings()
    for tag in(P.TAG,P.INC):
        doc=W.read(tag)
        for path,sha in doc['sources'].items():assert W.digest(path)==sha
        for folder,ext in(('scripts','py'),('results','json'),('results','md')):
            path=f'{folder}/{tag}.{ext}';sources[path]=W.digest(path)
    sources[f'scripts/{TAG}.py']=W.digest(f'scripts/{TAG}.py')
    data=normal(dict(status='PASS_RAW23_RAW25_CROSS_VERTEX_FAMILY',sources=sources,baseline_source=P.TAG,
        face_multipliers=base['face_multipliers'],restricted_face_multipliers=base['restricted_face_multipliers'],
        forest_records=records,zero_endpoint_contexts=zero,columns=columns,repairs=repairs,
        hypothetical_zero_event_columns=hypo,hypothetical_zero_event_repairs=ceilings,counts=dict(counts,paired_contexts=13,zero_endpoint_contexts=6,columns=len(columns)),
        local_proposal_count=(saved['local_proposal_count']if saved else len(proposal_keys)),global_optimizer_invocations=0,new_pair_outcomes=0,
        scope='Compact arithmetic H13 over fixed K=Q(sqrt5). Seven positive paired raw23/raw25 contexts and six zero-endpoint contexts from the frozen52-family. Native maps transport the already proved actual finite-support/triangle/quad/injection geometry. Complete local count/event vectors supply weighted two-vertex columns. The zero-event comparison is hypothetical and leaves pair ownership unchanged; its above12 restricted witness is a quantitative obstruction to this entire chosen family, not to other incidence geometries or global realizability.'))
    if saved:assert data==saved
    else:(ROOT/f'results/{TAG}.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    print(data['status'],data['counts'],'actual',[float(F(q['bound']))for q in repairs],
        'zero_event',[float(F(q['bound']))for q in ceilings],'seconds',round(time.monotonic()-started,2),flush=True)


if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');args=a.parse_args();run(args.verify)
