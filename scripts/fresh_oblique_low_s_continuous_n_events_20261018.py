#!/usr/bin/env python3
"""Reverse fixed low-S/R duties, then solve the continuous high-N right event.

The algebraic event equation covers a closed continuous overdomain; numerical
root/native tests are a bounded discovery diagnostic, not an exact exclusion.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_oblique_high_n_right_r_join_20261018 as H
ROOT,T,C,Q,F=H.ROOT,H.T,H.C,H.Q,H.F
TAG='fresh_oblique_low_s_continuous_n_events_20261018'
CONTINUOUS='fresh_oblique_continuous_high_n_high_s_interval_20261018'
LOW=(2,3,4,5)


def n_family(base,a):
    frame=[10,67,69];ff=np.array([base[i]for i in frame]);g=ff@T.M@ff.T;inv=np.linalg.inv(g)
    r0=np.array([-float(a),0.,0.]);r1=np.array([0.,0.,-1.])
    u=r0@inv@ff;v=r1@inv@ff;cv=ff@T.M@base[64]
    rem=base[64]-cv@inv@ff;sig=float(rem@T.M@rem);assert sig>0
    b=rem/np.sqrt(sig)
    rho=[float(2-u@T.M@u),float(-2*u@T.M@v),float(-v@T.M@v)]
    return dict(frame=frame,Gram=g.tolist(),inverse=inv.tolist(),old64_residual_square=sig,
                a=a,u=u.tolist(),v=v.tolist(),b=b.tolist(),rho=rho)


def quadratic(c):
    # A zero/constant event, or the double root caused by C=0, is not
    # discarded: it is unresolved. No such case occurs in this saved domain.
    scale=max(1.,max(abs(x)for x in c));eps=1e-11*scale
    if abs(c[2])<=eps:
        if abs(c[1])<=eps:
            return [],dict(status='unresolved_constant_or_zero_polynomial',scale=scale)
        return [-c[0]/c[1]],dict(status='linear',scale=scale)
    disc=c[1]*c[1]-4*c[2]*c[0];tol=1e-10*max(1.,c[1]*c[1],abs(4*c[2]*c[0]))
    if abs(disc)<=tol:return [],dict(status='unresolved_near_double_root',discriminant=disc,tolerance=tol)
    if disc<0:return [],dict(status='no_real_event',discriminant=disc)
    # Stable quadratic formula, preserving both roots.
    q=-.5*(c[1]+np.copysign(np.sqrt(disc),c[1]))
    roots=sorted([q/c[2],c[0]/q])
    return [float(z)for z in roots],dict(status='two_real_events',discriminant=disc)


def events(nf,rv):
    u,v,b=(np.array(nf[k])for k in('u','v','b'));rho=nf['rho']
    ab=[float(rv@T.M@z)for z in(u,v,b)];aa,bb,cc=ab
    poly=[aa*aa-cc*cc*rho[0],2*aa*bb-cc*cc*rho[1],bb*bb-cc*cc*rho[2]]
    roots,info=quadratic(poly);rows=[]
    for theta in roots:
        inside=np.sqrt(3)-1e-9<=theta<=2+1e-9
        row=dict(theta=theta,in_closed_high_domain=bool(inside))
        if inside:
            rem=sum(q*theta**k for k,q in enumerate(rho));row['norm_remainder']=rem
            if rem<=0:row['status']='unresolved_nonpositive_remainder'
            else:
                nv=u+theta*v+np.sqrt(rem)*b
                val=aa+bb*theta+cc*np.sqrt(rem);row.update(status='unsquared_match'if abs(val)<1e-6 else'squared_extraneous',
                    unsquared_residual=float(val),direct_RN=float(rv@T.M@nv),N_normal=nv.tolist(),
                    N_norm_residual=abs(float(nv@T.M@nv)-2),NX_relation=T.G.relation(-theta))
                if theta<2:
                    label=float(np.pi/np.arccos(theta/2));lo=int(np.floor(label))
                    row.update(effective_NX_label=label,adjacent_integer_labels=[lo,lo+1],
                        adjacent_theta_gaps=[float(theta-T.B.theta(lo)),float(T.B.theta(lo+1)-theta)])
        else:row['status']='outside_closed_high_domain'
        rows.append(row)
    return dict(coefficients=ab,polynomial=poly,root_analysis=info,roots=rows)


def r_recipes(s,nd,tt):
    frame=[69,71,72];sm=T.G.relation(s.dot(nd[71],s.ps(71,tt),nd[72],s.ps(72,tt)))
    assert isinstance(sm,int)and sm in LOW
    assert T.G.relation(s.dot(nd[69],s.ps(69,tt),nd[71],s.ps(71,tt)))==2
    recipes=[(pi,[2,m,r])for pi,m,r in product(range(7),LOW,LOW)
        if T.P.corner3(s.ps(69,tt),s.ps(71,tt),pi,2,2,m)is not None
        and T.P.corner3(s.ps(71,tt),s.ps(72,tt),pi,sm,m,r)is not None]
    return frame,sm,recipes


def run():
    names=(H.TAG,H.I.TAG,H.F.TAG,CONTINUOUS)
    docs={n:json.loads((ROOT/f'results/{n}.json').read_text())for n in names}
    for d in docs.values():
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    interface=H.H.old_frame_interface(s);assert interface==docs[H.TAG]['old_frame_interface']
    debit=T.H.deletion(s.primary,s.cur,C.DELETED);assert json.loads(json.dumps(debit))==docs[H.TAG]['exact_deletion_debit']
    pp=H.parents(s,docs[H.I.TAG],docs[H.F.TAG]);top=F.topology();records=[];accepted=[];all_events=[]
    for pre in pp:
        base,bt=pre['normals'],pre['types'];assert set(bt)=={67,68,69,71}and all(i not in base for i in(70,72,73))
        # Old64 is an explicit recipe witness, never a retained facet.
        anchors={**base,64:s.old[64]}
        nf=[n_family(anchors,a)for a in(0,1)]
        sb=H.D.bank(s,base,bt,'S')
        row=dict(index=pre['index'],origins=pre['origins'],instruction_key=pre['instruction_key'],reconstruction=pre['reconstruction'],
                 prefix_normals={str(i):base[i].tolist()for i in bt},prefix_profiles=bt,N_families=nf,S_bank=sb,R_banks=[])
        records.append(row)
        for cap in sb['caps']:
            nd={**base,72:np.array(cap['normal'])};tt={**bt,72:cap['profile']}
            frame,sm,rr=r_recipes(s,nd,tt)
            bank=dict(S=cap['index'],frame=frame,SM=sm,recipes=[],native_caps=[]);row['R_banks'].append(bank)
            for pi,ms in rr:
                rhs=[-T.B.theta(m)-s.off(pi,tt[i])for i,m in zip(frame,ms)]
                roots,info=T.N.norm_roots([nd[i]for i in frame],rhs,s.profiles[pi]['square'])
                recipe=dict(index=len(bank['recipes']),profile=pi,labels=ms,norm=info,branches=[]);bank['recipes'].append(recipe)
                if info['status'].startswith('unresolved'):s.unresolved.append(dict(parent=pre['index'],S=cap['index'],recipe=recipe['index'],norm=info))
                for branch,rv in roots:
                    vals,bad,cp=s.native(rv,pi,nd,tt)
                    error=max(abs(float(rv@T.M@rv)-s.profiles[pi]['square']),max(abs(float(rv@T.M@nd[i])-q)for i,q in zip(frame,rhs)))
                    br=dict(branch=branch,normal=rv.tolist(),native_products=sorted(vals.items()),blockers=[[i,vals[i]]for i in bad],center_product=cp,norm_residual=error)
                    recipe['branches'].append(br)
                    if error>2e-5:s.unresolved.append(dict(parent=pre['index'],S=cap['index'],recipe=recipe['index'],branch=branch,error=error));continue
                    if bad or cp>=-1e-8:continue
                    full={**nd,70:rv};ft={**tt,70:pi}
                    assert s.corner([69,71,70],full,ft)is not None and s.corner([71,70,72],full,ft)is not None
                    out=dict(recipe=recipe['index'],branch=branch,profile=pi,labels=ms,events=[]);bank['native_caps'].append(out)
                    for family in nf:
                        ev=events(family,rv);ev['a']=family['a'];out['events'].append(ev)
                        if ev['root_analysis']['status'].startswith('unresolved'):s.unresolved.append(dict(parent=pre['index'],S=cap['index'],recipe=recipe['index'],branch=branch,a=family['a'],event=ev))
                        for root in ev['roots']:
                            all_events.append(root)
                            if root['status']!='unsquared_match':continue
                            nv=np.array(root['N_normal']);nvals,nbad,ncp=s.native(nv,0,full,ft)
                            root.update(native_products=sorted(nvals.items()),blockers=[[i,nvals[i]]for i in nbad],center_product=ncp)
                            if nbad or ncp>=-1e-8:continue
                            chamber={**full,73:nv};ct={**ft,73:0}
                            # Every new wall is checked against all remaining walls again.
                            gates=[]
                            for i in ct:
                                vv,bb,cc=s.native(chamber[i],ct[i],{j:w for j,w in chamber.items()if j!=i},ct)
                                gates.append(dict(wall=i,blockers=[[j,vv[j]]for j in bb],center_product=cc))
                            root['all_new_wall_gates']=gates
                            if any(g['blockers']or g['center_product']>=-1e-8 for g in gates):continue
                            corners=[s.corner(f,chamber,ct)for f in top['triangles']]
                            root['corner_types']=corners
                            if not all(q is not None for q in corners):continue
                            score=Q.P.B.B.score(s,chamber,ct,[i for i in s.cur['actual_sides']if i not in C.DELETED]+list(ct),'continuous_N_event',[i for i in F.WALLS if i in chamber])
                            faces={f['wall']:f for f in score['whole_faces']}
                            exact={str(w):faces[w]['all_rays_finite']and sorted(sorted(q['active'])for q in faces[w]['rays'])==sorted(sorted(f)for f in top['triangles']if w in f)for w in(67,68,69,71)}
                            accepted.append(dict(parent=pre['index'],S=cap['index'],recipe=recipe['index'],branch=branch,a=family['a'],theta=root['theta'],
                                corner_types=corners,exact_whole_interior_caps=exact,whole_compact_faces={str(w):faces[w]['all_rays_finite']for w in faces},**score))
            print('PARENT',pre['index'],'S',cap['index'],'recipes',len(bank['recipes']),'native_R',len(bank['native_caps']),flush=True)
    banks=[b for r in records for b in r['R_banks']];recipes=[q for b in banks for q in b['recipes']];branches=[x for q in recipes for x in q['branches']]
    native=[x for b in banks for x in b['native_caps']]
    summary=dict(parents=len(records),parent_origins=sum(len(r['origins'])for r in records),S_recipes=sum(len(r['S_bank']['recipes'])for r in records),
        S_branches=sum(len(q['branches'])for r in records for q in r['S_bank']['recipes']),native_S=len(banks),R_recipes=len(recipes),R_branches=len(branches),
        native_R=len(native),event_equations=2*len(native),event_roots=len(all_events),event_status=dict(Counter(x['status']for x in all_events)),
        R_blocker_histogram=dict(Counter(str(tuple(i for i,v in b['blockers']))for b in branches)),accepted=len(accepted),new_cones=len(accepted),unresolved=len(s.unresolved))
    paths={Path(__file__),Path(H.__file__)};paths.update(ROOT/p for d in docs.values()for p in d['sources'])
    paths.update(ROOT/f'{folder}/{n}.{ext}'for n in names for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return dict(status='BOUNDED_LOW_S_CONTINUOUS_HIGH_N_EVENT_SEARCH',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        deleted_walls=list(C.DELETED),exact_deletion_debit=debit,old_frame_interface=interface,topology=top,
        numerical_snapshot=docs[H.TAG]['numerical_snapshot'],
        wall_count_breakdown=dict(projected_originals=52,fixed_original_whites=6,fixed_new_before_S=4,
            fixed_new_before_R=5,fixed_new_before_N=6),
        parameter_domain=[float(np.sqrt(3)),2.],N_families=[[2,2,'m>=6'],[3,2,'m>=6']],records=records,accepted=accepted,summary=summary,
        boundary_ledgers=s.ledgers,unresolved=s.unresolved,
        scope='Two exact-instruction parent classes/four origins. Complete low-S seven-profile cube and low-R RXright RM/RS2..5 duties, before continuous N event equations. Numerical normals, quadratic roots and native tests; not a certified interval/all-label exclusion. High RS/RM separate. Exact original support debit, no ambient volume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    assert args.verify or not path.exists();out=json.loads(json.dumps(run()))
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
