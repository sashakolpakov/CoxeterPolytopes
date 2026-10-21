#!/usr/bin/env python3
"""Partial rank-13 whole-forest corner-orbit model with exact certificate replay."""
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017'
PARAMS=['fresh_sqrt5_rank13_irreducible_quad_parameters_20261017','fresh_sqrt5_rank13_small_d_quad_parameters_20261017']
QUADS=['fresh_sqrt5_rank13_irreducible_joint_quads_20261017','fresh_sqrt5_rank13_small_d_joint_quads_20261017']
OLDTRI='fresh_sqrt5_rank13_triangle_third_costs_20261017'
NEWTRI='fresh_sqrt5_rank13_compensating_endpoints_triangles_20261017'

def dump(x):return json.dumps(x,separators=(',',':'))
def read(tag):
    d=json.loads((ROOT/f'results/{tag}.json').read_text())
    for p,h in d['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    return d

def domain():
    models=[(n,g)for n,g in C.models()if len(g)<=13]
    graphs=dict(models);finite={C.code(g):n for n,g in models};assert len(finite)==len(models)
    raw=[];all_forests=set();counts=Counter()
    def corner(g,marks):
        pieces=C.pieces(g);forest=tuple(sorted(finite[C.code(p)]for p in pieces))
        assert sum(len(p)for p in pieces)==13 and len(set(marks)&set(g))==2
        code=dump(sorted(C.code(p,marks)for p in pieces));all_forests.add(forest)
        return forest,code
    old=read(OLDTRI)
    for i,p in enumerate(old['triangle_patterns']):
        cs=p['selected_corners']+[p['third_corner']];keys=[]
        for c in cs:
            f=tuple(c['forest']);all_forests.add(f);keys.append((f,dump(c['marked_forest_code'])))
        raw.append((2,keys,dict(kind='old_triangle',source=OLDTRI,index=i)))
    counts['old_marked_triangle_patterns']=len(old['triangle_patterns'])
    new=read(NEWTRI)
    for ti,t in enumerate(new['targets']):
        base={i:{}for i in range(13)}
        for i,j,m in t['edges']:base[i][j]=base[j][i]=m
        for oi,o in enumerate(t['triangle_occurrences']):
            port=o['port'];other=o['other'];attach=t['ports'][port]['attachments'][o['attachment']]
            full={i:dict(row)for i,row in base.items()};full[13]={}
            for i,m in attach['contacts']:full[i][13]=full[13][i]=m
            if o['label']!=2:full[port][13]=full[13][port]=o['label']
            keys=[]
            for removed,pair in [(13,(port,other)),(port,(other,13)),(other,(port,13))]:
                g={i:{j:m for j,m in row.items()if j!=removed}for i,row in full.items()if i!=removed}
                keys.append(corner(g,pair))
            assert keys[0][0]==tuple(t['forest']) and keys[1][0]==tuple(o['destination'])
            assert keys[2]==(tuple(o['third_forest']),dump(o['third_marked_forest_code']))
            raw.append((2,keys,dict(kind='new_triangle',source=NEWTRI,target=ti,occurrence=oi)))
    counts['new_literal_triangle_occurrences']=sum(len(t['triangle_occurrences'])for t in new['targets'])
    for ptag,qtag in zip(PARAMS,QUADS):
        pars=read(ptag);qs=read(qtag)
        for p in pars['records']:
            if not p['accepted']:continue
            q=qs['records'][p['index']];assert q['index']==p['index'];row=q['row']
            labels=dict(zip(combinations(range(row['nodes']),2),row['labels']));keys=[]
            for pair in [(0,1),(1,2),(2,3),(0,3)]:
                nodes=list(row['white'])+list(pair)
                g=C.graph(nodes,labels)
                assert all(m!=0 for a in g.values()for m in a.values())
                keys.append(corner(g,pair))
            raw.append((1,keys,dict(kind='quad',source=qtag,parameters=ptag,record=p['index'])))
            counts['retained_colored_quadrilaterals']+=1
    vertices=[];variables=[];lookup={}
    for forest in sorted(all_forests):
        g={}
        for name in forest:
            h=graphs[name];offset=len(g)
            for i,row in h.items():g[offset+i]={offset+j:m for j,m in row.items()}
        assert len(g)==13
        by_orbit=defaultdict(list)
        for pair in combinations(range(13),2):
            f,code=corner(g,pair);assert f==forest;by_orbit[code].append(list(pair))
        vid=len(vertices);vvars=[]
        for code,pairs in sorted(by_orbit.items()):
            index=len(variables);lookup[forest,code]=index;vvars.append([index,len(pairs)])
            variables.append(dict(index=index,vertex=vid,forest=list(forest),marked_forest_code=json.loads(code),
                multiplicity=len(pairs),literal_pairs=pairs))
        assert sum(n for _,n in vvars)==78
        vertices.append(dict(index=vid,forest=list(forest),variables=vvars))
    faces=[];face_lookup={}
    for need,keys,provenance in raw:
        coefficients=tuple(sorted(Counter(lookup[k]for k in keys).items()))
        k=(need,coefficients)
        if k not in face_lookup:
            face_lookup[k]=len(faces);faces.append(dict(index=len(faces),required=need,
                coefficients=[list(x)for x in coefficients],provenance=[]))
        faces[face_lookup[k]]['provenance'].append(provenance)
    counts.update(vertices=len(vertices),variables=len(variables),literal_face_rows=len(raw),distinct_face_constraints=len(faces))
    paths=[f'scripts/{TAG}.py',f'scripts/{C.TAG}.py']+[f'results/{t}.json'for t in PARAMS+QUADS+[OLDTRI,NEWTRI]]
    return dict(sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},counts=dict(counts),
        variables=variables,vertices=vertices,faces=faces)

def dual_bound(d,alphas):
    loads=defaultdict(F)
    for i,a in alphas.items():
        assert a>=0
        for j,c in d['faces'][i]['coefficients']:loads[j]+=a*c
    beta={};maximizers={}
    for v in d['vertices']:
        maximum=max((loads[j]/m for j,m in v['variables']),default=F(0))
        if maximum:
            beta[v['index']]=maximum
            maximizers[v['index']]=[j for j,m in v['variables']if loads[j]/m==maximum]
    numerator=sum((d['faces'][i]['required']*a for i,a in alphas.items()),F(0))
    denominator=sum(beta.values(),F(0));assert denominator>0
    for var in d['variables']:
        assert loads[var['index']]<=beta.get(var['vertex'],F(0))*var['multiplicity']
    return numerator/denominator,beta,maximizers,numerator,denominator

def recover_dual(d,dual):
    candidates=[]
    for denominator in(10**6,10**8,10**10):
        alpha={i:F(round(-float(x)*denominator),denominator)for i,x in enumerate(dual)if x<0}
        alpha={i:a for i,a in alpha.items()if a>0}
        if alpha:candidates.append((dual_bound(d,alpha)[0],alpha))
    bound,alpha=max(candidates,key=lambda z:z[0])
    if bound>12:
        # Try a small literal subcore before retaining the full dual support.
        order=sorted(alpha,key=lambda i:(-alpha[i]*d['faces'][i]['required'],i))
        for n in [1,2,4,8,16,32,64,128,256,512,1024,len(order)]:
            if n>=len(order):break
            trial={i:alpha[i]for i in order[:n]}
            if dual_bound(d,trial)[0]>12:alpha=trial;break
        # A deterministic deletion pass reduces the retained exact support.
        for i in sorted(alpha,key=lambda j:(alpha[j]*d['faces'][j]['required'],j)):
            trial={j:a for j,a in alpha.items()if j!=i}
            if trial and dual_bound(d,trial)[0]>12:alpha=trial
    bound,beta,maximizers,num,den=dual_bound(d,alpha)
    return dict(kind='EXACT_RESTRICTED_MODEL_DUAL',bound=str(bound),objective_numerator=str(num),normalizer=str(den),
        face_multipliers=[[i,str(a)]for i,a in sorted(alpha.items())],
        vertex_multipliers=[[i,str(b)]for i,b in sorted(beta.items())],
        tight_orbits={str(i):v for i,v in maximizers.items()},excludes_target12=bound>12)

def recover_primal(d,x):
    best=None
    for denominator in(10**6,10**8,10**10):
        weights=[F(max(0,round(float(v)*denominator)),denominator)for v in x]
        scale=F(1)
        for f in d['faces']:
            total=sum((c*weights[j]for j,c in f['coefficients']),F(0))
            if not total:break
            scale=max(scale,F(f['required'])/total)
        else:
            weights=[w*scale for w in weights]
            budgets=[sum((m*weights[j]for j,m in v['variables']),F(0))for v in d['vertices']]
            cert=dict(kind='EXACT_PARTIAL_MODEL_PRIMAL',maximum=str(max(budgets)),scale=str(scale),
                weights=[[i,str(w)]for i,w in enumerate(weights)if w],vertex_budgets=list(map(str,budgets)),
                reaches_target12=max(budgets)<=12)
            if best is None or F(cert['maximum'])<F(best['maximum']):best=cert
    return best

def verify_certificate(d,c):
    if c['kind']=='EXACT_RESTRICTED_MODEL_DUAL':
        alpha={i:F(a)for i,a in c['face_multipliers']}
        bound,beta,maximizers,num,den=dual_bound(d,alpha)
        assert c['bound']==str(bound)and c['objective_numerator']==str(num)and c['normalizer']==str(den)
        assert c['vertex_multipliers']==[[i,str(b)]for i,b in sorted(beta.items())]
        assert c['tight_orbits']=={str(i):v for i,v in maximizers.items()}
        assert c['excludes_target12']==(bound>12)
    elif c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL':
        w={i:F(a)for i,a in c['weights']};assert all(x>=0 for x in w.values())
        for f in d['faces']:assert sum((a*w.get(i,F(0))for i,a in f['coefficients']),F(0))>=f['required']
        budgets=[sum((m*w.get(j,F(0))for j,m in v['variables']),F(0))for v in d['vertices']]
        assert c['vertex_budgets']==list(map(str,budgets))and c['maximum']==str(max(budgets))
        assert c['reaches_target12']==(max(budgets)<=12)
    else:raise AssertionError(c['kind'])

def optimize(d):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    N=len(d['variables']);rows=[];cols=[];values=[];rhs=[]
    for f in d['faces']:
        r=len(rhs);rhs.append(-f['required'])
        for j,c in f['coefficients']:rows.append(r);cols.append(j);values.append(-c)
    for v in d['vertices']:
        r=len(rhs);rhs.append(0)
        for j,m in v['variables']:rows.append(r);cols.append(j);values.append(m)
        rows.append(r);cols.append(N);values.append(-1)
    matrix=coo_matrix((values,(rows,cols)),shape=(len(rhs),N+1)).tocsr()
    objective=np.zeros(N+1);objective[-1]=1
    print('PARTIAL_CORNER_LP_START',d['counts'],flush=True)
    result=linprog(objective,A_ub=matrix,b_ub=np.array(rhs,dtype=float),bounds=(0,None),
        method='highs',options={'time_limit':300})
    solver=dict(status=int(result.status),message=str(result.message),success=bool(result.success),
        objective=None if result.fun is None else float(result.fun))
    certificates=[]
    if result.success:
        primal=recover_primal(d,result.x[:-1])
        if primal is not None:certificates.append(primal)
        certificates.append(recover_dual(d,result.ineqlin.marginals[:len(d['faces'])]))
    for c in certificates:verify_certificate(d,c)
    return solver,certificates

def run(verify):
    d=domain();out=ROOT/f'results/{TAG}.json'
    if verify:
        saved=json.loads(out.read_text())
        for k,v in d.items():assert saved[k]==v,k
        for c in saved['certificates']:verify_certificate(d,c)
        return saved
    solver,certificates=optimize(d)
    feasible=any(c.get('reaches_target12',False)for c in certificates)
    excluded=any(c.get('excludes_target12',False)for c in certificates)
    assert not(feasible and excluded)
    status='PASS_EXACT_PARTIAL_CORNER_RULE_AT12'if feasible else 'PASS_EXACT_PARTIAL_CORNER_MODEL_OBSTRUCTION'if excluded else'UNRESOLVED_PARTIAL_CORNER_MODEL'
    d.update(status=status,dimension=13,target='12',solver=solver,certificates=certificates,
        definition='An independent nonnegative payment is allowed for every unordered marked-root-pair orbit of each entire finite vertex forest appearing in the saved local constraints. All78 root pairs are counted at each such vertex. No baseline or old coefficient is fixed. The common maximum vertex budget is minimized.',
        coverage='The model includes every retained colored quad in the two complete ADE finite-anchor packages; every3558 old marked ADE-endpoint triangle constraint; and all1373 newly enumerated non-ADE triangle occurrences. These are the entire stated input catalogues, not complete rank13 face coverage.',
        scope='A partial local weighting model only. A feasible budget12 rule still requires all remaining triangles, mixed/non-ADE quadrilaterals, and faces at absent vertex types. An exact dual excludes only rules satisfying this overcomplete necessary local-geometry model; it is not a dimension bound, a compact realization claim, or a proof that stronger global incidence cannot work.')
    out.write_text(json.dumps(d,indent=2)+'\n');return d

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run(a.verify)
    print(d['status'],d['solver'],[(c['kind'],c.get('maximum',c.get('bound')))for c in d['certificates']],flush=True)
