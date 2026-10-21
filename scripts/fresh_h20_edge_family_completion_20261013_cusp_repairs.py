#!/usr/bin/env python3
"""Two exact affine completions of the new H20 cusps; the 96-wall cone still escapes.

Replay is deterministic and uses no numerical proposals.  The old loader imports
NumPy/SciPy, but no numerical result is used by any assertion below.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from fresh_h20_edge_family_completion_20261013_four_volume import setup,ROOT,read,sha,check
from fresh_dimension_bridge_20261013_verify import E,Z,O,parse,vec,dot,mv
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension,fdot,fmv,same,label,rank
from fresh_h20_edge_family_completion_20261013_exact import positive_solve,extra_pair

PREFIX='fresh_h20_edge_family_completion_20261013'
BASE='data/enumerations/fresh_dimension_bridge_20261013.json'
P93='data/enumerations/'+PREFIX+'_exact.json'
P94='data/enumerations/'+PREFIX+'_four_exact.json'
V94='data/enumerations/'+PREFIX+'_four_volume.json'
OUT='data/enumerations/'+PREFIX+'_cusp_repairs.json'
D18_IDS=[0,8,17,5,23,14,24,10,1,9,18,3,15,38,22,7,16,30]
D18_COEFFS=[2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,1,1]

def sumvec(terms):return tuple(sum((a*r[k]for a,r in terms),Z)for k in range(21))
def base_root(f,r):
    if f is B:return r
    roots={E(F(4)):E(F(2)),E(F(16)):E(F(4)),O:O}
    assert f.d in roots
    return tuple(a+roots[f.d]*b for a,b in r)

def evaluate_base(r,Q,roots):
    rc=mv(Q,r);q=dot(r,rc);rows=[]
    for j,(f,s)in enumerate(roots):
        M=Q if f is B else tuple(tuple(map(f.lift,row))for row in Q)
        p=fdot(f,s,rc if f is B else tuple(map(f.lift,rc)))
        n=fdot(f,s,fmv(f,M,s));k=label(f,p,q if f is B else f.lift(q),n)
        rows.append({'old':j,'pairing':f.encode(p),'kind':k})
    return q,rows

def link(delta,Q,hc,roots,require_bounded=False):
    dc=mv(Q,delta);assert dot(delta,dc)==Z and dot(delta,hc).sign()<0
    pairs=[];ids=[];rs=[]
    for i,(f,r)in enumerate(roots):
        p=fdot(f,r,dc if f is B else tuple(map(f.lift,dc)));assert f.sign(p)<=0
        pairs.append(f.encode(p))
        if f.sign(p)==0:ids.append(i);rs.append(base_root(f,r))
    gram=[[dot(a,mv(Q,b))for b in rs]for a in rs]
    assert rank(gram)==19
    unseen=set(range(len(ids)));components=[]
    while unseen:
        todo=[min(unseen)];unseen.remove(todo[0]);indices=[]
        while todo:
            i=todo.pop();indices.append(i)
            neighbors=[j for j in unseen if gram[i][j].sign()]
            for j in neighbors:unseen.remove(j);todo.append(j)
        indices.sort();M=[[gram[i][j]for j in indices]for i in indices];n=len(M)
        assert all(M[i][j].sign()<=0 for i in range(n)for j in range(n)if i!=j)
        full=positive_solve(M,[Z]*n);component={'walls':[ids[i]for i in indices],'gram':[[x.enc()for x in row]for row in M]}
        if full is not None:
            assert not require_bounded
            component.update({'kind':'spherical','rank':n,'positive_pivots':[x.enc()for x in full[1]]})
        else:
            assert n>=2
            result=positive_solve([row[:-1]for row in M[:-1]],[-row[-1]for row in M[:-1]])
            assert result is not None
            a,piv=result;kernel=a+[O]
            assert all(x.sign()>0 for x in kernel)and all(sum((x*y for x,y in zip(row,kernel)),Z)==Z for row in M)
            total=sumvec([(a,rs[j])for a,j in zip(kernel,indices,strict=True)])
            k=next(i for i,x in enumerate(delta)if x.sign());mu=total[k]/delta[k]
            assert mu.sign()>0 and total==tuple(mu*x for x in delta)
            component.update({'kind':'affine','rank':n-1,'positive_pivots':[x.enc()for x in piv],'positive_kernel':[x.enc()for x in kernel],'kernel_sum_multiple':mu.enc()})
        components.append(component)
    assert sum(c['rank']for c in components)==19
    bounded=all(c['kind']=='affine'for c in components)
    if require_bounded:assert bounded
    return {'ray':[x.enc()for x in delta],'norm':['0','0'],'center_pairing':dot(delta,hc).enc(),'all_wall_pairings':pairs,'incident_walls':ids,'quotient_rank':19,'components':components,'bounded_link':bounded}

def final_vertices(Q,roots,cap_rows):
    out=[]
    for row in cap_rows:
        v=row['actual_finite_vertex']
        if v is None:continue
        f=Extension(B,parse(row['field']['radicand']));x=tuple(map(f.parse,v['vector']));Qf=tuple(tuple(map(f.lift,rr))for rr in Q);xc=fmv(f,Qf,x)
        incident=[];pairs=[];survives=True
        for j,(g,r)in enumerate(roots):
            if g is B:p=fdot(f,xc,tuple(map(f.lift,r)));h=f
            else:
                h=Extension(f,f.lift(g.d));r2=tuple((f.lift(a),f.lift(b))for a,b in r)
                p=fdot(h,tuple(map(h.lift,xc)),r2)
            sg=h.sign(p);pairs.append(h.encode(p))
            if sg>0:survives=False
            if sg==0:incident.append(j)
        if survives:assert incident==sorted(v['incident_walls'])
        out.append({'created_at_wall':v['incident_walls'][-1],'incident_walls':incident,'survives_final_96':survives,'all_wall_pairings':pairs})
    return out

def run():
    _,_,_,Q,h,hc,roots=setup();assert len(roots)==94
    base=read(BASE);p93,p94=read(P93),read(P94)
    r91=base_root(*roots[91]);r93=base_root(*roots[93]);r41=base_root(*roots[41]);r27=base_root(*roots[27]);r9=base_root(*roots[9])
    delta1=sumvec([(O,r91),(E(F(1,3)),r41)]);delta2=sumvec([(O,r93),(E(F(1,3)),r27)])
    before=[link(d,Q,hc,roots)for d in(delta1,delta2)]
    assert [[c['rank']for c in L['components']if c['kind']=='spherical']for L in before]==[[1],[18]]
    beta=[sumvec([(O if i<24 else E(F(1,3)),base_root(*roots[i]))])for i in D18_IDS]
    M=[[dot(a,mv(Q,b))for b in beta]for a in beta]
    assert all(M[i][i]==E(F(2))for i in range(18))
    assert all(M[i][j]in(Z,-O)for i in range(18)for j in range(i))
    assert positive_solve(M,[Z]*18)is not None
    theta=sumvec([(E(F(a)),r)for a,r in zip(D18_COEFFS,beta,strict=True)])
    tc=mv(Q,theta);assert dot(theta,tc)==E(F(2))
    theta_pairs=[dot(b,tc)for b in beta]
    assert Counter(theta_pairs)==Counter({Z:17,O:1})
    new=[sumvec([(-O,r9),(E(F(2)),delta1)]),sumvec([(-O,theta),(E(F(5)),delta2)])]
    records=[]
    for i,r in enumerate(new,94):
        q,pairs=evaluate_base(r,Q,roots);assert q==E(F(2))
        center=dot(r,hc);assert center.sign()<0
        cusps=[dot(r,mv(Q,vec(c['ray'])))for c in base['cusp_rows']]
        assert len(cusps)==107 and all(x.sign()<0 for x in cusps)
        records.append({'wall':i,'normal':[x.enc()for x in r],'norm':q.enc(),'strict_center_pairing':center.enc(),'old_pair_checks':pairs,'pair_counts':dict(sorted(Counter(p['kind']for p in pairs).items())),'inherited_cusp_strict_pairings':[x.enc()for x in cusps]})
        roots.append((B,r))
    after=[link(d,Q,hc,roots,True)for d in(delta1,delta2)]
    c=read(V94)['certificate'];escape=check(vec(c['vector']),Q,hc,roots)
    assert all(escape[k]==c[k]for k in('vector','norm','center_pairing'))
    vertices=final_vertices(Q,roots,p93['new_caps']+p94['new_caps'])
    return {'status':'EXACT_96_WALL_H20_TWO_BOUNDED_CUSP_LINKS_WITH_STRICT_ESCAPE','source_sha256':{p:sha(p)for p in(BASE,P93,P94,V94)},'retained_walls':94,'wall_count':96,
            'formulas':{'delta1':'r91+r41/3','delta2':'r93+r27/3','wall94':'-r9+2*delta1','wall95':'-theta+5*delta2','theta_walls':D18_IDS,'theta_coefficients':D18_COEFFS,'theta_long_root_normalization':'divide roots with index>=24 by3','theta_norm':dot(theta,tc).enc(),'theta_pairings':[x.enc()for x in theta_pairs]},
            'new_caps':records,'cusps_before':before,'cusps_after':after,'stage_vertex_survival':vertices,'strict_escape':escape,
            'scope':'Two exact cusp-link affine completions retain all94 walls and all107 inherited cusps. Every new pair has Coxeter angle or is parallel/dotted; both new cusp links have rank19 and only positive affine components, hence are bounded. The same exact strict spacelike escape survives all96 walls, so this is an infinite-volume chamber, not a new finite-volume example. No claim that these are all cusps or vertices.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
    if a.verify:assert out==read(OUT)
    else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status']);print([c['pair_counts']for c in out['new_caps']]);print('stage vertices',[(v['created_at_wall'],v['survives_final_96'])for v in out['stage_vertex_survival']])
