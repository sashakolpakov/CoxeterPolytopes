#!/usr/bin/env python3
"""Subdivide the dotted Q--J edge; solve connector high labels by quadratics.

Numerical necessary-domain search on the15 frozen continuous event/J seeds.
No exact normal or ambient volume certificate is asserted.
"""
import argparse,json
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
import fresh_split30_continuous_u_20261018 as C

N,S,ROOT,M=C.N,C.S,C.ROOT,C.M
TAG='fresh_qj_connector_20261018'

def topology():
    old=S.topology();removed=[{30,69,71},{69,70,71}]
    fs=[t for t in old['triangles']if set(t)not in removed]
    assert len(fs)==10
    fs += [[30,69,72],[30,71,72],[69,70,72],[70,71,72]]
    bc=old['boundary'];es=Counter(tuple(sorted(e))for t in fs for e in combinations(t,2));vs=set(sum(fs,[]))
    assert (len(vs),len(es),len(fs))==(10,23,14)and (69,71)not in es
    assert set(es.values())=={1,2}
    assert {e for e,n in es.items()if n==1}=={tuple(sorted(e))for e in zip(bc,bc[1:]+bc[:1])}
    links=[]
    for v in sorted(vs):
        ee=[tuple(i for i in t if i!=v)for t in fs if v in t];ds=Counter(i for e in ee for i in e);seen={next(iter(ds))}
        while True:
            nxt=seen|{i for e in ee if set(e)&seen for i in e}
            if nxt==seen:break
            seen=nxt
        assert seen==set(ds)
        assert sorted(ds.values())==([1,1]+[2]*(len(ds)-2)if v in bc else[2]*len(ds))
        links.append(dict(vertex=v,edges=ee))
    return dict(boundary=bc,triangles=fs,vertices=10,edges=23,faces=14,
                removed_edge=[69,71],new_vertex=72,K_link=[30,69,70,71,30],links=links)

@lru_cache(None)
def domain(q30,rj):
    # Label6 represents an unbounded finite label slot. In any finite
    # norm2 triple such a bond isolates its third node. The existing RJ3/4
    # bounds KJ/KR by5; at most one of K30/KQ can be high.
    out=[]
    for pi in range(7):
      for ms in product(range(2,7),repeat=4):
        m30,mQ,mR,mJ=ms
        cs=[N.P.corner3(0,0,pi,q30,m30,mQ),N.P.corner3(0,0,pi,2,m30,mJ),
            N.P.corner3(0,0,pi,2,mQ,mR),N.P.corner3(0,0,pi,rj,mR,mJ)]
        if any(c is None for c in cs):continue
        high=[i for i,m in enumerate(ms)if m==6]
        assert len(high)<=1 and all(i in(0,1)for i in high)
        if high:assert pi==0
        out.append(dict(profile=pi,labels=list(ms),high_position=high[0]if high else None,corner_types=cs))
    return out

def rebuild(f,e,j):
    s=f.scout;nd=dict(f.base);ty={};u=e['u']
    k,qs,pi,m7,mP,rs,q30,q52=e['key']
    recipes=[(67,f.recipes[67],f.recipes[67]['sheet']),
             (68,f.recipes[68],f.recipes[68]['sheet']),
             (69,dict(frame=[30,52,67],profile=0,labels=[q30,q52,k]),qs),
             (70,dict(frame=[7,67,68],profile=pi,labels=[m7,2,mP]),rs)]
    for node,rec,sign in recipes:
        roots,status=f.roots(rec,u,nd,ty);assert roots and not status.startswith('unresolved')
        v=next(v for sg,v in roots if sg==sign);nd[node]=v;ty[node]=rec['profile']
        saved=next(r['normal']for r in e['native']if r['node']==node)
        assert max(abs(v-np.array(saved)))<2e-5
    rhs=[-N.B.theta(m)-(s.profiles[j['profile']]['offsets'][s.index[i]]if i in s.index else s.off(j['profile'],ty[i]))
         for i,m in zip(j['frame'],j['labels'])]
    roots,info=N.N.norm_roots([nd[i]for i in j['frame']],rhs,s.profiles[j['profile']]['square'])
    assert roots and not info['status'].startswith('unresolved')
    v=next(v for sign,v in roots if sign==j['branch']);assert max(abs(v-np.array(j['normal'])))<2e-5
    nd[71]=v;ty[71]=j['profile']
    for i in range(67,72):
        vals,bad,cp=s.native(nd[i],ty[i],{a:v for a,v in nd.items()if a!=i},ty)
        assert not bad and cp< -1e-8,(e['index'],j['index'],i,bad,cp)
    return nd,ty

def run():
    prior=json.loads((ROOT/f'results/{C.TAG}.json').read_text())
    for p,h in prior['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    f=C.Family();s=f.scout;disk=topology();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[]
    s.control=s.score(s.old,{},s.cur['actual_sides'],'control63');deletion=s.score(f.base,{},f.ids,'delete47')
    seeds=[];all_domains={};records=[];candidates=[];unresolved=[];counts=Counter();improvements=[]
    for e in prior['events']:
      for j in e['J_bank']['caps']:
        nd,ty=rebuild(f,e,j);ids=f.ids+list(range(67,72))
        sid=s.score(nd,ty,ids,'reconstructed_seed');assert sid is not None
        assert s.ledgers[sid]['counts']==dict(finite=34,exterior=6)and s.ledgers[sid]['open_ends']==10
        frame=[30,69,70,71];mat=np.array([nd[i]@M for i in frame]);cond=float(np.linalg.cond(mat))
        seed=dict(index=len(seeds),event=e['index'],J=j['index'],u=e['u'],key=e['key'],J_labels=j['labels'],
                  score=sid,condition=cond,frame=frame);seeds.append(seed)
        if not np.isfinite(cond)or cond>1e10:
            unresolved.append(dict(stage='four_frame_condition',seed=seed['index'],condition=cond));continue
        q30=N.G.relation(s.dot(nd[30],0,nd[69],0));rj=N.G.relation(s.dot(nd[70],0,nd[71],0))
        assert q30 in(2,3)and rj in(3,4)
        dom=domain(q30,rj);all_domains[f'{q30},{rj}']=dom;seed['domain_key']=f'{q30},{rj}'
        for di,rec in enumerate(dom):
            pi=rec['profile'];labels=rec['labels'];high=rec['high_position'];counts['recipes']+=1
            offsets=np.array([s.profiles[pi]['offsets'][s.index[i]]if i in s.index else s.off(pi,ty[i])for i in frame])
            rhs=-np.array([0 if high==i else N.B.theta(m)for i,m in enumerate(labels)])-offsets
            v0=np.linalg.solve(mat,rhs);q=s.profiles[pi]['square']
            base=dict(seed=seed['index'],domain_record=di,profile=pi,labels=labels,high_position=high)
            values=[]
            if high is None:
                residual=float(v0@M@v0)-q
                base.update(norm_residual=residual);counts['low_recipes']+=1
                if abs(residual)<2e-7:values.append((None,v0))
            else:
                counts['high_families']+=1;unit=np.zeros(4);unit[high]=-1;v1=np.linalg.solve(mat,unit)
                a=float(v1@M@v1);b=2*float(v0@M@v1);c=float(v0@M@v0)-q
                base['quadratic']=[a,b,c];scale=max(1.,abs(a),abs(b),abs(c))
                if abs(a)<1e-12*scale:
                    if abs(b)<1e-12*scale:
                        unresolved.append(dict(stage='high_degenerate_quadratic',**base));records.append(base);continue
                    roots=[-c/b]
                else:
                    disc=b*b-4*a*c
                    if abs(disc)<1e-12*scale*scale:
                        unresolved.append(dict(stage='high_tangent_quadratic',**base,discriminant=disc));records.append(base);continue
                    roots=[]if disc<0 else[(-b-np.sqrt(disc))/(2*a),(-b+np.sqrt(disc))/(2*a)]
                base['real_roots']=list(map(float,roots));base['roots_in_high_domain']=[]
                for theta in roots:
                    if np.sqrt(3)-1e-10<=theta<2:
                        base['roots_in_high_domain'].append(float(theta));values.append((float(theta),v0+theta*v1))
            records.append(base)
            for theta,v in values:
                counts['norm_solutions']+=1
                ms=list(labels)
                if high is not None:
                    m=N.G.relation(-theta)
                    if not isinstance(m,int)or m<6:
                        counts['high_non_native']+=1;base.setdefault('non_native_high_roots',[]).append(theta);continue
                    ms[high]=m
                rhs2=-np.array([N.B.theta(m)for m in ms])-offsets
                error=max(abs(float(v@M@v)-q),float(np.max(abs(mat@v-rhs2))))
                if error>2e-7:
                    unresolved.append(dict(stage='connector_residual',**base,theta=theta,error=error));continue
                vals,bad,cp=s.native(v,pi,nd,ty)
                if bad or cp>=-1e-8:
                    base.setdefault('native_failures',[]).append(dict(theta=theta,blockers=[[i,vals[i]]for i in bad],center_product=cp))
                    counts['native_or_center_rejected']+=1;continue
                full={**nd,72:v};ft={**ty,72:pi};sid2=s.score(full,ft,ids+[72],'connectorK')
                if sid2 is None:continue
                corners=[s.corner(t,full,ft)for t in disk['triangles']];faces=[]
                for node in range(67,73):
                    rays=[r for r in s.ledgers[sid2]['rays']if node in r['active']]
                    expected={tuple(sorted(t))for t in disk['triangles']if node in t}
                    faces.append(dict(node=node,finite_exactly_expected=all(r['kind']=='finite'for r in rays)and
                        {tuple(sorted(r['active']))for r in rays}==expected))
                done=all(c is not None for c in corners)and all(r['finite_exactly_expected']for r in faces)
                cap=dict(seed=seed['index'],domain_record=di,profile=pi,labels=ms,theta=theta,normal=v.tolist(),
                    norm_residual=error,native_products=sorted(vals.items()),center_product=cp,score=sid2,
                    finite_corners=corners,whole_cap_faces=faces,completed_disk=done)
                candidates.append(cap)
                if s.ledgers[sid2]['open_ends']<10:improvements.append(len(candidates)-1)
        print('CONNECTOR_SEED',seed['index'],'event',e['index'],'J',j['index'],'condition',cond,
              'recipes',len(dom),'candidates',len(candidates),'improvements',len(improvements),flush=True)
    assert len(seeds)==15 and counts==dict(recipes=855,low_recipes=847,high_families=8)
    assert not candidates and not improvements and not unresolved and not s.unresolved
    low=[r for r in records if r['high_position']is None]
    closest=min(low,key=lambda r:abs(r['norm_residual']))
    diagnosis=dict(closest_low_record=closest,minimum_absolute_low_norm_residual=abs(closest['norm_residual']),
        all_high_roots_outside_domain=all(not r['roots_in_high_domain']for r in records if r['high_position']is not None),
        next_changed_duty='Change the fixed P contacts to7/12 rather than repeating the same K labels. P currently has prescribed[3,3,2] in frame[7,12,L]; changing a genuinely finite alternative changes R and its QR-closing u. This is a proposed next domain, not tested here.')
    sources=dict(prior['sources'])
    for folder,ext in [('scripts','py'),('results','json'),('results','md')]:
        p=ROOT/f'{folder}/{C.TAG}.{ext}';sources[str(p.relative_to(ROOT))]=sha256(p.read_bytes()).hexdigest()
    sources[f'scripts/{TAG}.py']=sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(status='NUMERICAL_QJ_CONNECTOR_NECESSARY_DOMAIN_SCOUT',sources=sources,disk=disk,
        deleted_walls=[47],retained_old_wall_count=62,exact_deletion_debit=prior['exact_deletion_debit'],
        control_score=s.control,deletion_score=deletion,seeds=seeds,joint_finite_domains=all_domains,
        frame_order=[30,69,70,71],records=records,candidates=candidates,counts=dict(counts),
        boundary_ledgers=s.ledgers,improvements=improvements,unresolved=unresolved+s.unresolved,diagnosis=diagnosis,
        scope='Numerical reconstruction of all15 frozen continuous-u event/J seeds. Exact disk/deletion combinatorics, inherited complete seven-profile finite-corner classifier. Four prescribed native contacts determine K linearly; low2..5 recipes plus every symbolic>=6 family solved by a norm quadratic overtheta in[sqrt3,2). This is floating numerical root coverage at fixed seeds, not exact interval/all-parameter certification. Every accepted candidate must pass all62old/native/new contacts, strict center, complete targetcone/independentPD support agreement and whole-cap ray checks. No ambient finitevolume claim or exactnewnormal.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],'seeds',len(out['seeds']),out['counts'],'candidates',len(out['candidates']),
          'improvements',len(out['improvements']),'unresolved',len(out['unresolved']),flush=True)
