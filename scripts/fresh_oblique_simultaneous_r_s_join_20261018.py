#!/usr/bin/env python3
"""Independent R/S low banks on four frozen N geometries, then all joins."""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import numpy as np
import fresh_oblique_joint_n_r_pilot_20261018 as I
ROOT=I.ROOT;F=I.F;R=I.R;A=I.A;O=I.O;C=I.C;T=I.T;Q=I.Q;H=I.H
TAG='fresh_oblique_simultaneous_r_s_join_20261018'


def reconstruct(s,row,forward):
    prior=forward['records'][row['selection']['prefix']]
    nd={i:v for i,v in s.old.items()if i not in C.DELETED};tt={};checks=[]
    assert len(row['reconstruction'])==8
    for j,ch in enumerate(row['reconstruction']):
        if j==6:nd.pop(70);tt.pop(70)
        i=ch['wall'];saved=prior['old_R_normal']if j==3 else row['former_R']if j==6 else row['prefix_normals'][str(i)]
        v,proof=C.solve_saved(s,ch['frame'],ch['profile'],ch['labels'],ch['branch'],s.old if j<2 else nd,{}if j<2 else tt,saved)
        nd[i]=v;tt[i]=ch['profile'];checks.append(dict(wall=i,**proof))
    removed={str(i):dict(normal=nd.pop(i).tolist(),profile=tt.pop(i))for i in(70,72)}
    assert set(tt)=={67,68,69,71,73}
    for i in tt:
        _,bad,cp=s.native(nd[i],tt[i],{j:v for j,v in nd.items()if j!=i},tt);assert not bad and cp<-1e-8
    return nd,tt,checks,removed


def bank(s,nd,tt,kind):
    frame,duties,bonds,sy,rr=(I.recipes if kind=='R'else A.recipes)(s,nd,tt)
    row=dict(kind=kind,frame=frame,duties=duties,duty_bonds=bonds,frame_profiles=[s.ps(i,tt)for i in frame],
        symbolic_high_families=[list(z)for z in sy if z[2]],recipes=[],caps=[])
    seen=set();digest=sha256();hist=Counter()
    for pi,ms in rr:
        rhs=[-T.B.theta(m)-(s.profiles[pi]['offsets'][s.index[i]]if i in s.index else s.off(pi,tt[i]))for i,m in zip(frame,ms)]
        roots,info=T.N.norm_roots([nd[i]for i in frame],rhs,s.profiles[pi]['square'])
        rec=dict(index=len(row['recipes']),profile=pi,labels=ms,norm=info,branches=[]);row['recipes'].append(rec)
        if info['status'].startswith('unresolved'):
            s.unresolved.append(dict(kind=kind,recipe=rec['index'],norm=info));continue
        for branch,v in roots:
            vals,bad,cp=s.native(v,pi,nd,tt)
            err=max(abs(float(v@T.M@v)-s.profiles[pi]['square']),max(abs(float(v@T.M@nd[i])-q)for i,q in zip(frame,rhs)))
            b=dict(branch=branch,normal=v.tolist(),native_products=sorted(vals.items()),blockers=[[i,vals[i]]for i in bad],center_product=cp,norm_residual=err)
            rec['branches'].append(b);hist[str(tuple(bad))]+=1
            digest.update(json.dumps([rec['index'],branch,pi,ms,b['native_products'],cp],separators=(',',':')).encode())
            if err>2e-5:s.unresolved.append(dict(kind=kind,recipe=rec['index'],branch=branch,error=err));continue
            if bad or cp>=-1e-8:continue
            key=(pi,*np.round(v,8))
            if key in seen:b['duplicate']=True;continue
            seen.add(key);cap=dict(index=len(row['caps']),recipe=rec['index'],profile=pi,labels=ms,**b)
            b['cap']=cap['index'];row['caps'].append(cap)
    row['branch_histogram']=dict(hist);row['complete_branch_digest']=digest.hexdigest()
    return row


def run():
    names=(I.TAG,F.TAG);docs={n:json.loads((ROOT/f'results/{n}.json').read_text())for n in names}
    for d in docs.values():
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    old=docs[I.TAG];forward=docs[F.TAG];assert old['summary']['native_R']==80 and len(old['records'])==4
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    interface=H.old_frame_interface(s);assert interface==old['old_frame_interface']
    debit=T.H.deletion(s.primary,s.cur,C.DELETED);assert json.loads(json.dumps(debit))==old['exact_deletion_debit']
    top=F.topology();ids=[i for i in s.cur['actual_sides']if i not in C.DELETED]+[67,68,69,71,73]
    score=lambda nd,tt,ii,stage:Q.P.B.B.score(s,nd,tt,ii,stage,[i for i in F.WALLS if i in nd])
    control=score(s.old,{},s.cur['actual_sides'],'original63_control');assert(control['finite_rays'],control['open_ends'])==(26,10)
    records=[];accepted=[];relations=Counter()
    for src in old['records']:
        nd,tt,steps,removed=reconstruct(s,src,forward)
        row=dict(index=len(records),selection=src['selection'],reconstruction=steps,removed=removed,
            prefix_normals={str(i):nd[i].tolist()for i in tt},prefix_profiles=tt,
            before=score(nd,tt,ids,'both_R_S_removed'),R_bank=bank(s,nd,tt,'R'),S_bank=bank(s,nd,tt,'S'),pairs=[],accepted=[])
        records.append(row);digest=sha256()
        def pair_full(a,b,diagnostic=False):
            full={**nd,70:np.array(a['normal']),72:np.array(b['normal'])};ft={**tt,70:a['profile'],72:b['profile']}
            for i in ft:
                _,bad,cp=s.native(full[i],ft[i],{j:v for j,v in full.items()if j!=i},ft);assert not bad and cp<-1e-8
            sc=score(full,ft,ids+[70,72],'diagnostic_pair'if diagnostic else'finite_MRS_join')
            faces={f['wall']:f for f in sc['whole_faces']};cs=[s.corner(f,full,ft)for f in top['triangles']]
            exact={str(w):faces[w]['all_rays_finite']and sorted(sorted(q['active'])for q in faces[w]['rays'])==sorted(sorted(f)for f in top['triangles']if w in f)for w in(67,68,69,71)}
            return dict(R=a['index'],S=b['index'],diagnostic_only=diagnostic,corner_types=cs,all_disk_corners=all(z is not None for z in cs),
                exact_whole_interior_caps=exact,whole_compact_faces={str(w):faces[w]['all_rays_finite']for w in faces},**sc)
        for a in row['R_bank']['caps']:
            for b in row['S_bank']['caps']:
                value=s.dot(np.array(a['normal']),a['profile'],np.array(b['normal']),b['profile']);rel=T.G.relation(value)
                cs=None
                if isinstance(rel,int):cs=s.corner([71,70,72],{**nd,70:np.array(a['normal']),72:np.array(b['normal'])},{**tt,70:a['profile'],72:b['profile']})
                rec=dict(R=a['index'],S=b['index'],product=value,relation=rel,finite_MRS=cs)
                row['pairs'].append(rec);relations[str(rel)]+=1;digest.update(json.dumps(rec,separators=(',',':')).encode())
                if cs is not None:
                    out=pair_full(a,b);row['accepted'].append(out);accepted.append(dict(prefix=row['index'],**out))
        row['complete_pair_digest']=digest.hexdigest();native=[p for p in row['pairs']if p['relation']=='dotted'or isinstance(p['relation'],int)]
        row['closest_native_diagnostic']=None
        if native:
            best=max(native,key=lambda p:(p['product'],-p['R'],-p['S']))
            row['closest_native_diagnostic']=pair_full(row['R_bank']['caps'][best['R']],row['S_bank']['caps'][best['S']],True)
        print('PREFIX',row['index'],'banks',len(row['R_bank']['caps']),len(row['S_bank']['caps']),'pairs',len(row['pairs']),'finiteMRS',len(row['accepted']),flush=True)
    summary=dict(prefixes=4,R_recipes=sum(len(r['R_bank']['recipes'])for r in records),S_recipes=sum(len(r['S_bank']['recipes'])for r in records),
        branches=sum(len(q['branches'])for r in records for k in('R_bank','S_bank')for q in r[k]['recipes']),
        R_caps=sum(len(r['R_bank']['caps'])for r in records),S_caps=sum(len(r['S_bank']['caps'])for r in records),pairs=sum(len(r['pairs'])for r in records),
        pair_relations=dict(relations),finite_MRS=len(accepted),completed_disks=sum(c['score']is not None and c['all_disk_corners']and all(c['exact_whole_interior_caps'].values())for c in accepted),
        diagnostic_cones=sum(r['closest_native_diagnostic']is not None for r in records),unresolved=len(s.unresolved))
    paths={Path(__file__),Path(I.__file__)};paths.update(ROOT/p for p in old['sources'])
    paths.update(ROOT/f'{folder}/{n}.{ext}'for n in names for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    return dict(status='BOUNDED_SIMULTANEOUS_R_S_LOW_JOIN',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        deleted_walls=list(C.DELETED),retained_old_walls=58,wall_count_breakdown=dict(projected_originals=52,fixed_original_whites=6,
            fixed_new_walls=5,per_candidate_explicit_gates=57,total_fixed_facets=63),
        exact_deletion_debit=debit,old_frame_interface=interface,topology=top,control=control,records=records,accepted=accepted,
        summary=summary,boundary_ledgers=s.ledgers,unresolved=s.unresolved,numerical_snapshot=old['numerical_snapshot'],
        scope='Complete independent seven-profile label2..5 R and S banks and all pair joins on four frozen N geometries only. Both former R/S removed before candidate gates. Exact old19/12/1 debit; new normals and target cones numerical. Three closest-native diagnostic cones are not closure candidates. High families/other prefixes unsearched; no ambientH9 volume certificate.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    assert a.verify or not path.exists();out=json.loads(json.dumps(run()))
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['summary'],flush=True)
