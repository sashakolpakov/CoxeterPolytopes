#!/usr/bin/env python3
"""Twenty high-N geometries with R-X and R-N right, then all low S joins.

Two exact-instruction parent classes retain all four selected origins.
This is a bounded numerical search, not an all-high-label exclusion.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_oblique_simultaneous_r_s_join_20261018 as D
ROOT=D.ROOT;I=D.I;F=D.F;A=D.A;C=D.C;T=D.T;Q=D.Q;H=D.H
TAG='fresh_oblique_high_n_right_r_join_20261018'
HIGH=(6,7,8,10,12);LOW=(2,3,4,5)


def parents(s,old,forward):
    out=[];seen={}
    for row in old['records']:
        nd,tt,steps,removed=D.reconstruct(s,row,forward)
        removed['73']=dict(normal=nd.pop(73).tolist(),profile=tt.pop(73))
        chain=row['reconstruction'][:5]
        key=json.dumps([[z['wall'],z['frame'],z['profile'],z['labels'],z['branch']]for z in chain],separators=(',',':'))
        origin=dict(index=row['index'],forward_prefix=row['selection']['prefix'])
        if key in seen:
            prev=out[seen[key]]
            assert tt==prev['types']and all(np.max(abs(nd[i]-prev['normals'][i]))<1e-9 for i in tt)
            prev['origins'].append(origin);continue
        seen[key]=len(out)
        out.append(dict(index=len(out),origins=[origin],instruction_key=key,reconstruction=steps,
            removed=removed,normals=nd,types=tt))
    assert len(out)==2 and [len(p['origins'])for p in out]==[3,1]
    return out


def constrained_recipes(s,nd,tt):
    frame=[69,71,73];duties=[(0,1),(0,2)];bonds=[2]*9
    for a,b in duties:
        m=T.G.relation(s.dot(nd[frame[a]],s.ps(frame[a],tt),nd[frame[b]],s.ps(frame[b],tt)))
        assert isinstance(m,int);bonds[3*a+b]=m
    assert bonds[1]==2 and bonds[2]in HIGH
    rr=[(pi,[2,m,2])for pi in range(7)for m in LOW
        if all(T.P.corner3(s.ps(frame[a],tt),s.ps(frame[b],tt),pi,bonds[3*a+b],ms[a],ms[b])is not None
            for ms in([2,m,2],)for a,b in duties)]
    full=[(pi,list(ms))for pi in range(7)for ms in product(LOW,repeat=3)
        if all(T.P.corner3(s.ps(frame[a],tt),s.ps(frame[b],tt),pi,bonds[3*a+b],ms[a],ms[b])is not None for a,b in duties)]
    assert sorted(rr)==sorted(full)
    return frame,duties,bonds,rr


def rbank(s,nd,tt):
    frame,duties,bonds,rr=constrained_recipes(s,nd,tt)
    row=dict(kind='R',frame=frame,duties=duties,duty_bonds=bonds,frame_profiles=[s.ps(i,tt)for i in frame],recipes=[],caps=[])
    seen=set();digest=sha256();hist=Counter()
    for pi,ms in rr:
        rhs=[-T.B.theta(m)-s.off(pi,tt[i])for i,m in zip(frame,ms)]
        roots,info=T.N.norm_roots([nd[i]for i in frame],rhs,s.profiles[pi]['square'])
        rec=dict(index=len(row['recipes']),profile=pi,labels=ms,norm=info,branches=[]);row['recipes'].append(rec)
        if info['status'].startswith('unresolved'):
            s.unresolved.append(dict(kind='R',recipe=rec['index'],norm=info));continue
        for branch,v in roots:
            vals,bad,cp=s.native(v,pi,nd,tt)
            err=max(abs(float(v@T.M@v)-s.profiles[pi]['square']),max(abs(float(v@T.M@nd[i])-q)for i,q in zip(frame,rhs)))
            b=dict(branch=branch,normal=v.tolist(),native_products=sorted(vals.items()),blockers=[[i,vals[i]]for i in bad],center_product=cp,norm_residual=err)
            rec['branches'].append(b);hist[str(tuple(bad))]+=1
            digest.update(json.dumps([rec['index'],branch,pi,ms,b['native_products'],cp],separators=(',',':')).encode())
            if err>2e-5:s.unresolved.append(dict(kind='R',recipe=rec['index'],branch=branch,error=err));continue
            if bad or cp>=-1e-8:continue
            key=(pi,*np.round(v,8))
            if key in seen:b['duplicate']=True;continue
            seen.add(key);cap=dict(index=len(row['caps']),recipe=rec['index'],profile=pi,labels=ms,**b)
            b['cap']=cap['index'];row['caps'].append(cap)
    row['branch_histogram']=dict(hist);row['complete_branch_digest']=digest.hexdigest();return row


def run():
    names=(D.TAG,I.TAG,F.TAG);docs={n:json.loads((ROOT/f'results/{n}.json').read_text())for n in names}
    for doc in docs.values():
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert docs[D.TAG]['summary']['pairs']==2844 and docs[D.TAG]['summary']['finite_MRS']==0
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    interface=H.old_frame_interface(s);assert interface==docs[D.TAG]['old_frame_interface']
    debit=T.H.deletion(s.primary,s.cur,C.DELETED);assert json.loads(json.dumps(debit))==docs[D.TAG]['exact_deletion_debit']
    top=F.topology();oldids=[i for i in s.cur['actual_sides']if i not in C.DELETED]
    score=lambda nd,tt,stage:Q.P.B.B.score(s,nd,tt,oldids+list(tt),stage,[i for i in F.WALLS if i in nd])
    control=Q.P.B.B.score(s,s.old,{},s.cur['actual_sides'],'original63_control',F.WALLS[:8])
    assert(control['finite_rays'],control['open_ends'])==(26,10)
    pp=parents(s,docs[I.TAG],docs[F.TAG]);records=[];accepted=[];relations=Counter()
    for pre in pp:
        base,bt=pre['normals'],pre['types'];frame=[10,67,69]
        _,duties,bonds,sy,_=F.recipes(s,base,bt)
        high=[list(z)for z in sy if z[2]]
        assert len(high)==2 and all(z[0]==0 and list(z[2])==[2]and list(z[1])in([2,2,6],[3,2,6])for z in high)
        G=np.array([[float(base[i]@T.M@base[j])for j in frame]for i in frame]);inv=np.linalg.inv(G)
        eig=np.linalg.eigvalsh(G);assert sum(eig>1e-8)==2 and sum(eig< -1e-8)==1
        row=dict(index=pre['index'],origins=pre['origins'],instruction_key=pre['instruction_key'],reconstruction=pre['reconstruction'],removed=pre['removed'],
            prefix_normals={str(i):base[i].tolist()for i in bt},prefix_profiles=bt,before=score(base,bt,'N_R_S_removed'),
            frame=frame,duties=duties,duty_bonds=bonds,frame_profiles=[s.ps(i,bt)for i in frame],symbolic_high_families=high,
            projected_Gram=G.tolist(),projected_inverse=inv.tolist(),eigenvalues=eig.tolist(),recipes=[],caps=[])
        records.append(row);digest=sha256()
        for n10,m in product((3,2),HIGH):
            labels=[n10,2,m];rhs=np.array([-T.B.theta(z)-(s.profiles[0]['offsets'][s.index[i]]if i in s.index else s.off(0,bt[i]))for i,z in zip(frame,labels)])
            b=rhs.copy();b[2]=0;v=np.array([0.,0.,-1.]);coeff=[float(2-b@inv@b),float(-2*b@inv@v),float(-v@inv@v)]
            theta=T.B.theta(m);res=coeff[0]+coeff[1]*theta+coeff[2]*theta*theta
            assert abs(res-(2-rhs@inv@rhs))<1e-8
            roots,info=T.N.norm_roots([base[i]for i in frame],rhs,s.profiles[0]['square'])
            rec=dict(index=len(row['recipes']),profile=0,labels=labels,theta=theta,residual_square_polynomial=coeff,residual_square=res,norm=info,branches=[]);row['recipes'].append(rec)
            if info['status'].startswith('unresolved'):s.unresolved.append(dict(kind='N',parent=pre['index'],recipe=rec['index'],norm=info));continue
            for branch,nv in roots:
                vals,bad,cp=s.native(nv,0,base,bt)
                err=max(abs(float(nv@T.M@nv)-2),max(abs(float(nv@T.M@base[i])-q)for i,q in zip(frame,rhs)))
                br=dict(branch=branch,normal=nv.tolist(),native_products=sorted(vals.items()),blockers=[[i,vals[i]]for i in bad],center_product=cp,norm_residual=err)
                rec['branches'].append(br);digest.update(json.dumps([rec['index'],branch,0,labels,br['native_products'],cp],separators=(',',':')).encode())
                if err>2e-5:s.unresolved.append(dict(kind='N',parent=pre['index'],recipe=rec['index'],branch=branch,error=err));continue
                if bad or cp>=-1e-8:continue
                nd={**base,73:nv};tt={**bt,73:0}
                assert all(s.corner([frame[a],frame[b],73],nd,tt)is not None for a,b in duties)
                cap=dict(index=len(row['caps']),recipe=rec['index'],profile=0,labels=labels,**br,before=score(nd,tt,'high_N_R_S_removed'),
                    R_bank=rbank(s,nd,tt),S_bank=D.bank(s,nd,tt,'S'),pairs=[],accepted=[])
                br['cap']=cap['index'];row['caps'].append(cap);pdigest=sha256()
                def full_pair(a,b,diagnostic=False):
                    full={**nd,70:np.array(a['normal']),72:np.array(b['normal'])};ft={**tt,70:a['profile'],72:b['profile']}
                    for i in ft:
                        _,fail,center=s.native(full[i],ft[i],{j:w for j,w in full.items()if j!=i},ft);assert not fail and center<-1e-8
                    sc=score(full,ft,'high_N_closest_pair'if diagnostic else'high_N_finite_join');faces={f['wall']:f for f in sc['whole_faces']}
                    cs=[s.corner(f,full,ft)for f in top['triangles']]
                    exact={str(w):faces[w]['all_rays_finite']and sorted(sorted(q['active'])for q in faces[w]['rays'])==sorted(sorted(f)for f in top['triangles']if w in f)for w in(67,68,69,71)}
                    return dict(R=a['index'],S=b['index'],diagnostic_only=diagnostic,corner_types=cs,all_disk_corners=all(z is not None for z in cs),
                        exact_whole_interior_caps=exact,whole_compact_faces={str(w):faces[w]['all_rays_finite']for w in faces},**sc)
                for a in cap['R_bank']['caps']:
                    for b in cap['S_bank']['caps']:
                        value=s.dot(np.array(a['normal']),a['profile'],np.array(b['normal']),b['profile']);rel=T.G.relation(value);cs=None
                        if isinstance(rel,int):cs=s.corner([71,70,72],{**nd,70:np.array(a['normal']),72:np.array(b['normal'])},{**tt,70:a['profile'],72:b['profile']})
                        pair=dict(R=a['index'],S=b['index'],product=value,relation=rel,finite_MRS=cs);cap['pairs'].append(pair);relations[str(rel)]+=1
                        pdigest.update(json.dumps(pair,separators=(',',':')).encode())
                        if cs is not None:
                            done=full_pair(a,b);cap['accepted'].append(done);accepted.append(dict(parent=row['index'],N=cap['index'],**done))
                cap['complete_pair_digest']=pdigest.hexdigest();cap['closest_native_diagnostic']=None
                native=[p for p in cap['pairs']if p['relation']=='dotted'or isinstance(p['relation'],int)]
                if native:
                    best=max(native,key=lambda p:(p['product'],-p['R'],-p['S']))
                    cap['closest_native_diagnostic']=full_pair(cap['R_bank']['caps'][best['R']],cap['S_bank']['caps'][best['S']],True)
                print('PARENT',row['index'],'N',labels,'banks',len(cap['R_bank']['caps']),len(cap['S_bank']['caps']),'pairs',len(cap['pairs']),'joins',len(cap['accepted']),flush=True)
        row['complete_branch_digest']=digest.hexdigest()
    caps=[c for r in records for c in r['caps']]
    summary=dict(parent_origins=4,distinct_parents=2,N_recipes=sum(len(r['recipes'])for r in records),N_branches=sum(len(q['branches'])for r in records for q in r['recipes']),native_N=len(caps),
        R_recipes=sum(len(c['R_bank']['recipes'])for c in caps),S_recipes=sum(len(c['S_bank']['recipes'])for c in caps),
        R_branches=sum(len(q['branches'])for c in caps for q in c['R_bank']['recipes']),S_branches=sum(len(q['branches'])for c in caps for q in c['S_bank']['recipes']),
        R_caps=sum(len(c['R_bank']['caps'])for c in caps),S_caps=sum(len(c['S_bank']['caps'])for c in caps),pairs=sum(len(c['pairs'])for c in caps),pair_relations=dict(relations),
        finite_MRS=len(accepted),completed_disks=sum(c['score']is not None and c['all_disk_corners']and all(c['exact_whole_interior_caps'].values())for c in accepted),
        diagnostic_cones=sum(c['closest_native_diagnostic']is not None for c in caps),unresolved=len(s.unresolved))
    paths={Path(__file__),Path(D.__file__)};paths.update(ROOT/p for p in docs[D.TAG]['sources'])
    paths.update(ROOT/f'{folder}/{n}.{ext}'for n in names for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return dict(status='BOUNDED_HIGH_N_RIGHT_R_JOIN',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        deleted_walls=list(C.DELETED),exact_deletion_debit=debit,old_frame_interface=interface,topology=top,control=control,records=records,accepted=accepted,summary=summary,boundary_ledgers=s.ledgers,unresolved=s.unresolved,
        high_labels=list(HIGH),necessary_R_rights=dict(NX_squared_lower_bound=3,other_nonright_squared_lower_bound=1,determinant_upper_bound=0,RX=2,RN=2),
        wall_count_breakdown=dict(projected_originals=52,fixed_original_whites=6,fixed_new_before_N=4,fixed_new_before_R_S=5,explicit_N_gates=56,explicit_R_or_S_gates=57),
        numerical_snapshot=docs[D.TAG]['numerical_snapshot'],
        scope='Two complete-instruction parent classes with allfour origins; profile0 N families[3,2,m]/[2,2,m] at m6/7/8/10/12, both norm branches. Allretained gates before complete constrained seven-profile low-R and independent low-S banks and Cartesian joins. Exact necessary right-contact lemma conditional on native finite corner; allnew normals/cones numerical. No continuous/all-high-label or ambientvolume claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    assert a.verify or not path.exists();out=json.loads(json.dumps(run()))
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
