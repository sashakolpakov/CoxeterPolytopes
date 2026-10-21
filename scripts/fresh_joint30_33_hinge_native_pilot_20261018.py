#!/usr/bin/env python3
"""Move the shared30 hinge jointly with33 while rebuilding compact10 duties.

A bounded two-stage numerical family; all old/new native and cone gates remain.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
from itertools import product
from itertools import combinations
from fractions import Fraction as F
import numpy as np
import fresh_wall66_seam_subdivision_20261018 as C

T=C.T;ROOT=C.ROOT;TAG='fresh_joint30_33_hinge_native_pilot_20261018'


def old_frame_interface(s):
    """Recover all actual-side A2 profiles; compare with saved native enclosures.

    The inherited scout only listed profiles for its previous boundary. Its
    native gate already used every original wall offset; this local extension
    is needed for newly chosen old-frame walls such as10 in corner classifiers.
    Native normals below are reconstructed from the frozen identity
    facet_point_i = h - (h.n_i)n_i/2, without numerical normal input.
    """
    p=s.primary
    dec=lambda q:tuple(map(F,q))
    add=lambda a,b:(a[0]+b[0],a[1]+b[1])
    sub=lambda a,b:(a[0]-b[1],a[1]-b[0])
    def mul(a,b):
        z=[x*y for x in a for y in b];return min(z),max(z)
    def div(a,b):
        assert b[0]*b[1]>0
        return mul(a,(1/b[1],1/b[0]))
    enc=lambda q:list(map(str,q))
    h=list(map(dec,p['common_center']));g=[[dec(q)for q in row]for row in p['computational_gram']]
    normals={r['wall']:[div(mul((F(2),F(2)),sub(a,dec(b))),dec(p['center_products'][str(r['wall'])]))for a,b in zip(h,r['point'])]for r in p['facet_witnesses']}
    gates={r['wall']:r for r in p['native_gates']}
    def native_interval(i,j):
        if i==j:return(F(2),F(2)),'original_norm2'
        if 66 in(i,j):
            row=gates[j if i==66 else i]
            if row['label']==2:return(F(0),F(0)),'original66_prescribed_right'
            return dec(row['product']),'original66_native_gate'
        a,b=normals[i],normals[j];z=(F(0),F(0))
        for k,x in enumerate(a):
            for l,y in enumerate(b):z=add(z,mul(mul(x,g[k][l]),y))
        return z,'original_facet_center_identity'
    def compare(x,q):
        assert float(q[0])-2e-5<=x<=float(q[1])+2e-5,(x,enc(q))
    rows=[]
    for i in s.cur['actual_sides']:
        k=s.index[i];a=s.profiles[1]['offsets'][k];b=s.profiles[4]['offsets'][k]
        contacts=[-(2*a-b),-(2*b-a)]
        matches=[j for j,r in enumerate(s.profiles)if max(abs(x-y)for x,y in zip(contacts,r['contacts']))<2e-5]
        assert len(matches)==1,(i,contacts,matches)
        pi=matches[0];s.types[i]=pi
        square=float(s.old[i]@T.M@s.old[i]);assert abs(square-s.profiles[pi]['square'])<2e-5
        checks=[]
        for j in[0,1,2,3,5,6]:
            q,origin=native_interval(i,j);x=0. if j<4 else contacts[j-5];compare(x,q)
            checks.append(dict(white=j,numerical_product=x,exact_native_enclosure=enc(q),origin=origin))
        rows.append(dict(wall=i,profile=pi,stored_A2_contacts=contacts,intrinsic_square=square,original_native_checks=checks))
    assert s.types[7]==4 and s.types[10]==1 and s.types[60]==6
    pairs=[]
    for i,j in combinations([7,10,55,66],2):
        x=s.dot(s.old[i],s.ps(i,{}),s.old[j],s.ps(j,{}));q,origin=native_interval(i,j);compare(x,q)
        pairs.append(dict(walls=[i,j],profiles=[s.ps(i,{}),s.ps(j,{})],native_product=x,
            native_relation=T.G.relation(x),exact_native_enclosure=enc(q),origin=origin))
    assert next(r for r in pairs if r['walls']==[7,10])['native_relation']==2
    return dict(actual_side_profiles=rows,all_old_frame_pairs=pairs,
        scope='Exact rational interval calculations from frozen original native data compared with numerical stored-profile/normal products. The profile choice and numerical comparison use tolerance2e-5; no new physical enclosure or exact new-normal certificate.')


def face(bd,wall):
    rr=[r for r in bd['rays']if wall in r['active']]
    return dict(wall=wall,rays=rr,all_rays_finite=bool(rr)and all(r['kind']=='finite'for r in rr))


def first_bank(s,nd,ids):
    frame=[7,10,66];duties=[(0,1),(0,2)];bonds=[2]*9
    for a,b in duties:
        m=T.G.relation(s.dot(nd[frame[a]],s.ps(frame[a],{}),nd[frame[b]],s.ps(frame[b],{})))
        assert isinstance(m,int);bonds[3*a+b]=m
    symbolic=T.E.recipes(tuple(s.ps(i,{})for i in frame),tuple(bonds),tuple(duties))
    recipes=[(pi,list(ms))for pi,ms,high in symbolic if not high and all(2<=m<=5 for m in ms)]
    direct=[(pi,list(ms))for pi in range(7)for ms in product(range(2,6),repeat=3)
        if all(T.P.corner3(s.ps(frame[a],{}),s.ps(frame[b],{}),pi,bonds[3*a+b],ms[a],ms[b])is not None for a,b in duties)]
    assert sorted(recipes)==sorted(direct)
    tests=[];caps=[];seen=set()
    for pi,ms in recipes:
        rhs=[-T.B.theta(m)-s.profiles[pi]['offsets'][s.index[i]]for i,m in zip(frame,ms)]
        roots,info=T.N.norm_roots([nd[i]for i in frame],rhs,s.profiles[pi]['square'])
        row=dict(index=len(tests),profile=pi,labels=ms,norm=info,branches=[]);tests.append(row);s.counts['norm_equations']+=1
        if info['status'].startswith('unresolved'):
            s.unresolved.append(dict(stage='first_norm',recipe=row['index'],info=info));continue
        for sign,v in roots:
            s.counts['norm_branches']+=1
            err=max(abs(float(v@T.M@v)-s.profiles[pi]['square']),max(abs(float(v@T.M@nd[i])-x)for i,x in zip(frame,rhs)))
            vals,bad,cp=s.native(v,pi,nd,{})
            br=dict(branch=sign,norm_residual=err,center_product=cp,blockers=[[i,vals[i]]for i in bad]);row['branches'].append(br)
            if err>2e-5:
                s.unresolved.append(dict(stage='first_residual',recipe=row['index'],branch=sign,error=err));continue
            if bad or cp>=-1e-8:continue
            key=(pi,)+tuple(np.round(v,8))
            if key in seen:br['duplicate']=True;continue
            seen.add(key);nn={**nd,67:v};tt={67:pi}
            assert all(s.corner([frame[a],frame[b],67],nn,tt)is not None for a,b in duties)
            before=len(s.unresolved);sid=s.score(nn,tt,ids+[67],'replace33')
            diagnostic=None if sid is not None else s.unresolved[before:]
            # A native intermediate can remain provisional until new30 resolves
            # the deletion's nonsimple exterior ray. Final acceptance is strict.
            br['cap']=len(caps);caps.append(dict(index=len(caps),recipe=row['index'],profile=pi,labels=ms,
                branch=sign,normal=v.tolist(),native_products=sorted(vals.items()),center_product=cp,norm_residual=err,
                score=sid,cone_validated=sid is not None,cone_diagnostic=diagnostic))
    return dict(frame=frame,duties=duties,complete_low_recipe_count=len(recipes),recipes=tests,caps=caps,
        symbolic_high_families_not_tested=sum(bool(h)for _,_,h in symbolic))


def rank_score(s,c):
    if c['score']is not None:
        b=s.ledgers[c['score']];return b['open_ends'],-b['counts'].get('finite',0)
    ds=c['cone_diagnostic'];assert len(ds)==1
    pd=ds[0]['PD_incidence']
    return pd.get('open_ends',100000),-pd.get('finite',0)


def run():
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    interface=old_frame_interface(s)
    debit=T.H.deletion(s.primary,s.cur,[30,33]);assert debit['lost_finite_count']==15 and debit['compact_faces_by_dimension']=={2:10,3:1}
    assert all(8 not in v for v in debit['lost_finite_supports'])
    nd={i:v for i,v in s.old.items()if i not in(30,33)};ids=[i for i in s.cur['actual_sides']if i not in(30,33)]
    control=s.score(s.old,{},s.cur['actual_sides'],'control63');deletion=s.score(nd,{},ids,'delete30_33')
    assert control is not None
    deletion_diagnostic=list(s.unresolved)if deletion is None else []
    first=first_bank(s,nd,ids)
    for c in first['caps']:
        c['restores_old33_numerically']=c['profile']==0 and float(np.max(abs(np.array(c['normal'])-s.old[33])))<2e-6
    # All first-bank caps are scored. At most eight genuinely changed normals
    # enter the second bank; every old33 control also enters and is flagged.
    ranked=sorted([c for c in first['caps']if not c['restores_old33_numerically']],key=lambda c:(
        *rank_score(s,c),c['profile'],c['labels'],c['branch']))
    selected=[c['index']for c in ranked[:8]]+[c['index']for c in first['caps']if c['restores_old33_numerically']]
    records=[];accepted=[];corner_failures=[]
    print('HINGE_FIRST',len(first['caps']),'selected',selected,'unresolved',len(s.unresolved),flush=True)
    for ci in selected:
        a=first['caps'][ci];nn={**nd,67:np.array(a['normal'])};tt={67:a['profile']};actual=ids+[67]
        frame_checks=[]
        for x,y in combinations([10,55,67],2):
            val=s.dot(nn[x],s.ps(x,tt),nn[y],s.ps(y,tt))
            if y==67:
                saved=dict(a['native_products'])[x];origin='first_cap_full_native_gate'
            else:
                saved=next(r['native_product']for r in interface['all_old_frame_pairs']if r['walls']==[x,y]);origin='original_native_comparison_above'
            assert abs(val-saved)<2e-5 and T.G.relation(val)==T.G.relation(saved)
            frame_checks.append(dict(walls=[x,y],profiles=[s.ps(x,tt),s.ps(y,tt)],native_product=val,reference_product=saved,origin=origin))
        # Only10/67 is a mandatory pair in this frame. No10/55 or55/67
        # corner is invented;7/55/new30 is checked after all native gates.
        b=C.lowbank(s,[10,55,67],[(0,2)],nn,tt,actual,68,'replace30')
        ri=len(records);records.append(dict(first_cap=ci,frame_native_checks=frame_checks,bank=b))
        for c in b['caps']:
            n2={**nn,68:np.array(c['normal'])};t2={**tt,68:c['profile']}
            duties=[[10,67,68],[7,55,68]]
            failed=[f for f in duties if s.corner(f,n2,t2)is None]
            if failed:
                corner_failures.append(dict(record=ri,cap=c['index'],failed_corners=failed,score=c['score']));continue
            # Recheck both normals in the joint complete native list.
            for new in(67,68):
                vals,bad,cp=s.native(n2[new],t2[new],{i:v for i,v in n2.items()if i!=new},t2)
                assert not bad and cp<-1e-8
            bd=s.ledgers[c['score']]
            f10=face(bd,10);f55=face(bd,55)
            accepted.append(dict(record=ri,cap=c['index'],score=c['score'],old33_restored=a['restores_old33_numerically'],
                old30_restored=c['profile']==0 and float(np.max(abs(n2[68]-s.old[30])))<2e-6,
                open_ends=bd['open_ends'],finite_rays=bd['counts'].get('finite',0),
                whole10=f10,whole55=f55,new33_face=face(bd,67),new30_face=face(bd,68),
                exterior_supports=[r['active']for r in bd['rays']if r['kind']=='exterior'],
                open_edges=[e for e in bd['edges']if e['kind']in('one_open','two_open')]))
        print('HINGE_SECOND',ci,'recipes',b['complete_low_recipe_count'],'caps',len(b['caps']),'accepted',len(accepted),'unresolved',len(s.unresolved),flush=True)
    paths=[Path(__file__),Path(C.__file__),Path(T.__file__),Path(T.H.__file__),Path(T.E.__file__),Path(T.N.__file__),Path(T.P.__file__),Path(T.G.__file__),T.H.PARENT,T.H.BOUNDARY,
        Path(T.H.A.__file__),Path(T.H.A.C.__file__),T.H.A.INPUT,T.H.A.C.PARENT,Path(T.S.__file__),Path(T.S.C.__file__)]
    return dict(status='BOUNDED_NUMERICAL_JOINT30_33_HINGE_PILOT',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        deleted_walls=[30,33],retained_old_walls=61,exact_deletion_debit=debit,control_score=control,deletion_score=deletion,deletion_diagnostic=deletion_diagnostic,
        old_frame_interface=interface,low_labels=[2,3,4,5],first_bank=first,second_selection=selected,
        selection_rule='All first-bank native prefixes attempted complete cone/PD scoring; provisional exterior-incidence failures retained only as search inputs, ranked by saved PD diagnostic. Best eight genuinely changed33 normals ranked by openends,negativefinite,profile,labels,branch, plus all old33 numerical controls. Each selected prefix gets complete seven-profile2..5 cube for frame10/55/67 with only10/67 finite duty; final7/55/68 checked separately.',
        records=records,accepted=accepted,required_corner_failures=corner_failures,boundary_ledgers=s.ledgers,
        counts=dict(s.counts),histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n)for(a,b,c,d),n in sorted(s.hist.items())],
        blockers=[dict(stage=a,walls=list(b),count=n)for(a,b),n in s.blockers.most_common()],unresolved=s.unresolved,
        scope='Joint30/33 physical reconstruction changes the common boundary hinge and the formerly fixed33 frame; exact old15/10/1 support debit. Bounded low2..5 seven-profile banks with selected second-stage prefixes, both norm branches, all61 retained old and allmutual native/commoncenter gates. First-stage native prefixes with unresolved exterior incidence are provisional only; every admitted final pair requires clean full numericalS6 cone and independentPD incidence. Whole10/55 target-face records are explicit. No high-label/continuous completeness, exact new normals, ambient finite census, finite-volume or dimension claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'first',len(out['first_bank']['caps']),'accepted',len(out['accepted']),'unresolved',len(out['unresolved']),flush=True)
