#!/usr/bin/env python3
"""Change the old58/59 bond, with exact old deletion debit and numerical gates."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_release58_66_seven_boundary_20261017 as T
ROOT,M=T.ROOT,T.M;TAG='fresh_wall59_bond_replacement_20261018'
FRAME=[55,30,58]
LABEL30=[2,3,4,5,6,8,12]
LABEL58=[2,3,4,5,6,7,8,9,10,11,12,16,24]

def run():
    s=T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[]
    control=s.score(s.old,{},s.cur['actual_sides'],'control63')
    debits=[];bases={};stages=[]
    for release in([59],[59,60],[59,60,61]):
        debit=T.H.deletion(s.primary,s.cur,release)
        assert debit['lost_finite_count']==7 and debit['compact_faces_by_dimension']=={2:6,3:1}
        debit['lost_supports_containing8']=sum(8 in v for v in debit['lost_finite_supports'])
        assert debit['lost_supports_containing8']==0
        debits.append(debit);key=','.join(map(str,release))
        nd={i:v for i,v in s.old.items()if i not in release};ids=[i for i in s.cur['actual_sides']if i not in release]
        sid=s.score(nd,{},ids,'deletion'+key);assert sid is not None
        bases[key]=(nd,ids);stages.append(dict(deleted_walls=release,retained_old_walls=63-len(release),score=sid))
    assert debits[0]['lost_finite_supports']==debits[1]['lost_finite_supports']
    assert debits[0]['compact_faces_touched']==debits[1]['compact_faces_touched']
    assert debits[0]['lost_finite_supports']==debits[2]['lost_finite_supports']
    assert debits[0]['compact_faces_touched']==debits[2]['compact_faces_touched']
    frame=np.array([s.old[i]for i in FRAME]);records=[];accepted=[];counts=Counter();blockers=Counter();unresolved=[]
    for m30,m58 in product(LABEL30,LABEL58):
        labels=[2,m30,m58]
        # Profile0, right to55; both displayed boundary corners are finite.
        assert T.P.corner3(0,0,0,2,2,m30)is not None
        assert T.P.corner3(0,0,0,2,2,m58)is not None
        rhs=[-T.B.theta(m)-s.profiles[0]['offsets'][s.index[i]]for i,m in zip(FRAME,labels)]
        rr,info=T.N.norm_roots(frame,rhs,s.profiles[0]['square']);counts['norm_equations']+=1
        row=dict(index=len(records),profile=0,frame=FRAME,labels=labels,norm=info,branches=[]);records.append(row)
        if info['status'].startswith('unresolved'):
            unresolved.append(dict(stage='norm',record=row['index'],info=info));continue
        for sign,v in rr:
            counts['norm_branches']+=1
            err=max(abs(float(v@M@v)-2),max(abs(float(v@M@s.old[i])-r)for i,r in zip(FRAME,rhs)))
            if err>2e-5:
                unresolved.append(dict(stage='norm_residual',record=row['index'],branch=sign,error=err));continue
            vals,bad,cp=s.native(v,0,bases['59'][0],{})
            rec=dict(branch=sign,normal=v.tolist(),native_products=sorted(vals.items()),blockers=bad,center_product=cp,norm_error=err)
            row['branches'].append(rec);blockers[tuple(bad)]+=1
            if cp>=-1e-8:counts['wrong_common_center']+=1;continue
            for key,(nd,ids)in bases.items():
                bad2=[i for i in bad if i in nd]
                if bad2:continue
                nn={**nd,67:v};tt={67:0};sid=s.score(nn,tt,ids+[67],'replacement'+key)
                if sid is None:continue
                bd=s.ledgers[sid];rays55=[r for r in bd['rays']if 55 in r['active']]
                whole55=all(r['kind']=='finite'for r in rays55)and{tuple(sorted(r['active']))for r in rays55}=={(7,30,55),(7,55,58),(30,55,67),(55,58,67)}
                accepted.append(dict(index=len(accepted),record=row['index'],branch=sign,release=list(map(int,key.split(','))),profile=0,labels=labels,normal=v.tolist(),score=sid,
                    center_product=cp,norm_error=err,whole55_quadrilateral=whole55,old60_native=vals[60],
                    finite_old_contacts=[[i,T.G.relation(x)]for i,x in vals.items()if i in nd and isinstance(T.G.relation(x),int)]))
                counts['accepted_'+key]+=1
        print('REPLACEMENT_RECIPE',m30,m58,'accepted',len(accepted),'unresolved',len(unresolved)+len(s.unresolved),flush=True)
    baseline=[r for r in accepted if r['release']==[59]and r['labels']==[2,5,9]]
    assert len(baseline)==1
    assert np.max(abs(np.array(baseline[0]['normal'])-s.old[59]))<2e-6
    low=[r for r in accepted if r['labels'][2]<9]
    paths=[Path(__file__),Path(T.__file__),Path(T.H.__file__),Path(T.N.__file__),Path(T.P.__file__),Path(T.G.__file__),T.H.PARENT,T.H.BOUNDARY,
        ROOT/'results/fresh_wall66_single_edge_conic_20261018.json',ROOT/'results/fresh_nonarithmetic_alternative66_release59_20261017.json',
        ROOT/'scripts/fresh_nonarithmetic_pi8_quad_completion_20261017.py',ROOT/'results/fresh_nonarithmetic_pi8_quad_completion_20261017.json']
    return dict(status='BOUNDED_NUMERICAL_WALL59_BOND_REPLACEMENT',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        exact_deletion_debits=debits,control_score=control,deletion_stages=stages,profile=0,frame=FRAME,fixed_contact55=2,label30_samples=LABEL30,label58_samples=LABEL58,
        records=records,accepted=accepted,low58_label_accepted_indices=[r['index']for r in low],counts=dict(counts),blockers=[dict(walls=list(k),count=v)for k,v in blockers.most_common()],
        boundary_ledgers=s.ledgers,histograms=[dict(stage=a,finite=b,exterior=c,open_ends=d,count=n)for(a,b,c,d),n in sorted(s.hist.items())],
        unresolved=unresolved+s.unresolved,
        scope='Purposeful replacement of59 changes old58/59 label9. Profile0 right55, bounded native labels to30/58; release59,59/60,and59/60/61 are compared with identical exact old deletion debit. Every accepted replacement passes all retained native products/commoncenter and complete numericalS6 cone/independentPD supports. Lost old certificates and numerical open-end changes do not imply improved volume. No exact newnormal or finite-volume theorem.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'accepted',len(out['accepted']),'low58',len(out['low58_label_accepted_indices']),'unresolved',len(out['unresolved']),flush=True)
