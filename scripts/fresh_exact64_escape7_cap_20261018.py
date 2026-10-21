#!/usr/bin/env python3
"""Native cap of the old7/30/47 escape on the exact64-wall seed."""
import argparse,json
from collections import Counter
from itertools import combinations
from hashlib import sha256
from pathlib import Path
import numpy as np
import fresh_remote60_exact_original_extension_20261018 as A
import fresh_remote60_compact_face_bridge_20261018 as D
ROOT=A.ROOT;TAG='fresh_exact64_escape7_cap_20261018';FRAME=[7,30,47];NEW=76
V,add,sub,mul,div,neg,sqrt,dot=A.V,A.add,A.sub,A.mul,A.div,A.neg,A.sqrt,A.dot
enc=A.enc;WHITE=A.WHITE


def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())


def numerical():
    prior=read(D.TAG)
    for path,h in prior['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    s=D.T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    assert D.P.H.H.old_frame_interface(s)==prior['old_frame_interface']
    u=next(b for b in prior['remote_bank']['recipes'][1]['branches']if b['native_admissible'])
    nd={**s.old,75:np.array(u['normal'])};tt={75:0};rows=[]
    assert [s.ps(i,tt)for i in FRAME]==[4,0,0]
    assert [D.T.G.relation(s.dot(nd[7],4,nd[i],0))for i in(30,47)]==[5,3]
    for m in(2,3,4):
        duties=[D.T.P.corner3(4,0,0,b,2,x)for b,x in((5,2),(3,m))]
        assert all(x is not None for x in duties)
        roots,info=D.T.N.norm_roots([nd[i]for i in FRAME],[0,0,-D.T.B.theta(m)],2)
        assert len(roots)==2
        row=dict(label=m,duties=duties,norm=info,branches=[]);rows.append(row)
        for sign,v in roots:
            prods,bad,cp=s.native(v,0,nd,tt)
            row['branches'].append(dict(branch=sign,normal=v.tolist(),products=sorted(prods.items()),
                blockers=[[i,prods[i]]for i in bad],center_product=cp,native=not bad and cp<0))
    assert all([b['native']for b in r['branches']]==[True,False]for r in rows)
    cones=[]
    for row in [None]+rows:
        full=nd.copy()
        if row:full[NEW]=np.array(row['branches'][0]['normal'])
        cone=D.T.S.C.detailed_cone(full,list(full));assert not cone['warnings']
        section7=[r for r in cone['rays']if 7 in r['active']]
        section60=[r for r in cone['rays']if 60 in r['active']]
        assert len(section60)==3 and all(r['kind']=='finite'for r in section60)
        if row:assert len(section7)==11 and all(r['kind']=='finite'for r in section7)
        cones.append(dict(label=None if row is None else row['label'],cone=cone,section7=section7,section60=section60))
    return dict(frame=FRAME,profile=0,declared_labels=[2,3,4],recipes=rows,cones=cones,
        snapshot=dict(metric=D.T.M.tolist(),center=s.center.tolist(),old_normals={str(i):v.tolist()for i,v in nd.items()}))


def physical():
    d=read(A.TAG);p=read(A.PRIMARY)
    for x in(d,p):
        for path,h in x['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    z=d['physical_input'];metric=[list(map(V,r))for r in z['metric']];h=list(map(V,z['center']));h2=V(z['center_square'])
    roots={int(i):list(map(V,v))for i,v in z['roots'].items()};roots[75]=list(map(V,d['normal']));ids=sorted(roots)
    hc={int(i):V(v)for i,v in z['center_products'].items()};hc[75]=V(d['center_product'])
    wi=[list(map(V,r))for r in d['white_inverse']];wc={int(i):list(map(V,v))for i,v in d['white_covectors'].items()};wc[75]=[V(0)]*6
    hw=[hc[w]for w in WHITE];hp={i:sub(hc[i],dot(hw,wi,wc[i]))for i in ids};hp2=sub(h2,dot(hw,wi,hw))
    native={};identities=[]
    for i in FRAME:
        native[i]={}
        for j in ids:
            raw=dot(roots[i],metric,roots[j]);v,pr=A.refine(raw,i,j);native[i][j]=v
            if pr:identities.append(pr)
    assert next(r['identity']for r in identities if r['pair']==[7,30])==5
    assert next(r['identity']for r in identities if r['pair']==[7,47])==3
    assert native[30][47][1]<-2
    g={i:{j:sub(native[i][j],dot(wc[i],wi,wc[j]))for j in ids}for i in FRAME}
    gram=[[g[i][j]for j in FRAME]for i in FRAME];inv,det=A.P.inverse(gram)
    assert gram[0][0][0]>0 and A.P.determinant([r[:2]for r in gram[:2]])[0]>0 and det[1]<0
    hf=[hp[i]for i in FRAME];q=sub(hp2,dot(hf,inv,hf));assert q[0]>0
    # All-right member chosen; the other two label variants are numerical.
    rhs=[V(0)]*3;rem=V(2);fac=neg(sqrt(div(rem,q)))
    products={j:mul(fac,sub(hp[j],dot(hf,inv,[g[i][j]for i in FRAME])))for j in ids}
    for j in WHITE+FRAME:products[j]=V(0)
    cp=mul(fac,q);assert cp[1]<0
    pj={j:A.vector_sub(roots[j],A.linear(A.row_times(wc[j],wi),[roots[w]for w in WHITE]))for j in ids}
    ph=A.vector_sub(h,A.linear(A.row_times(hw,wi),[roots[w]for w in WHITE]))
    res=A.vector_sub(ph,A.linear(A.row_times(hf,inv),[pj[i]for i in FRAME]));v=[mul(fac,x)for x in res]
    direct={j:dot(v,metric,roots[j])for j in ids};schur=products.copy()
    products={j:(max(x[0],direct[j][0]),min(x[1],direct[j][1]))for j,x in products.items()}
    assert all(x[0]<=x[1]for x in products.values())
    dotted=set(ids)-set(WHITE+FRAME);assert len(dotted)==55 and all(products[i][1]<-2 for i in dotted),[(i,enc(products[i]))for i in dotted if products[i][1]>=-2]
    norm=dot(v,metric,v);assert norm[0]<=2<=norm[1]
    oldw={i:sub(cp,mul(div(hc[i],V(2)),products[i]))for i in ids}
    neww={i:sub(hc[i],mul(div(cp,V(2)),products[i]))for i in ids}
    nsq=sub(h2,div(mul(cp,cp),V(2)));assert nsq[1]<0 and all(x[1]<0 for x in oldw.values())and all(x[1]<0 for x in neww.values())
    for i in FRAME:g[i][NEW]=products[i]
    g[NEW]={**products,NEW:V(2)};hp[NEW]=cp;pj[NEW]=v
    vertices=[]
    for support in([7,30,NEW],[7,47,NEW]):
        gv=[[g[i][j]for j in support]for i in support];vi,piv=A.I.inverse(gv);assert all(x[0]>0 for x in piv)
        cv=[hp[i]for i in support];coeff=A.row_times(cv,vi);point=A.vector_sub(ph,A.linear(coeff,[pj[i]for i in support]))
        sq=sub(hp2,dot(cv,vi,cv));assert sq[1]<0
        pp={}
        for j in ids+[NEW]:
            val=sub(hp[j],dot(cv,vi,[g[i][j]for i in support]))
            if j in WHITE+support:val=V(0)
            else:assert val[1]<0,(support,j,enc(val))
            pp[j]=val
        vertices.append(dict(support=WHITE+support,section_support=support,Gram=gv,inverse=vi,pivots=piv,point=point,square=sq,products=pp))
    trials=[]
    for subset in combinations(sorted(WHITE+FRAME),8):
        support=sorted([NEW,*subset])
        if 30 in subset and 47 in subset:trials.append(dict(support=support,dotted_pair=[30,47]))
        else:
            assert support in[sorted(r['support'])for r in vertices]
            trials.append(dict(support=support,finite=True))
    assert len(trials)==9 and sum('finite'in r for r in trials)==2
    oldescape=p['strict_escape'];oldprod=dot(v,metric,list(map(V,oldescape['point'])))
    assert oldprod[0]>0
    # A fresh exact ray at the retained exterior30/59/75 region, followed by
    # a strict center perturbation. All gates below are physical inequalities.
    ef=[30,59,75];eg={}
    for i in ef:
        eg[i]={}
        for j in ids:
            raw=dot(roots[i],metric,roots[j]);val,proof=A.refine(raw,i,j)
            eg[i][j]=sub(val,dot(wc[i],wi,wc[j]))
    ei,edet=A.P.inverse([[eg[i][j]for j in ef]for i in ef]);assert edet[1]<0
    ec=[hp[i]for i in ef];er=A.vector_sub(ph,A.linear(A.row_times(ec,ei),[pj[i]for i in ef]));er=[neg(x)for x in er]
    escape=None
    for power in range(1,50):
        epsilon=V(A.F(1,2**power));point=[add(x,mul(epsilon,y))for x,y in zip(er,h)]
        sq=dot(point,metric,point);future=dot(point,metric,h)
        pp={i:dot(point,metric,r)for i,r in{**roots,NEW:v}.items()}
        if sq[0]>0 and future[1]<0 and all(x[1]<0 for x in pp.values()):
            escape=dict(frame=ef,power=power,point=point,square=sq,future=future,products=pp);break
    return dict(old_walls=ids,new_wall=NEW,deleted_original_walls=[],frame=FRAME,profile=0,labels=[2,2,2],
        identity_refinements=identities,frame_Gram=gram,frame_inverse=inv,frame_determinant=det,
        center_complement_square=q,factor=fac,normal=v,direct_norm=norm,center_product=cp,
        schur_products=schur,direct_products=direct,native_products=products,
        old_facet_witness_products=oldw,new_facet_witness_products=neww,new_facet_witness_square=nsq,
        new_finite_vertices=vertices,new_finite_support_trials=trials,previous_escape_new_product=oldprod,
        remaining_strict_escape=escape,
        summary=dict(total_walls=65,old_finite_vertices_retained=201,new_finite_vertices=2,complete_finite_vertices=203,
            strict_dotted_contacts=55,right_contacts=9,infinite_volume_proved=escape is not None))


def run():
    n=numerical();print('NUMERICAL_PASS3native/6roots;4completecones',flush=True)
    p=physical();print('PHYSICAL_PASS65walls203finite;escape',p['remaining_strict_escape']is not None,flush=True)
    sources={f'{folder}/{stem}.{ext}':sha256((ROOT/f'{folder}/{stem}.{ext}').read_bytes()).hexdigest()
        for stem,pairs in[(TAG,[('scripts','py')]),(A.TAG,[('scripts','py'),('results','json'),('results','md')]),
            (D.TAG,[('scripts','py'),('results','json'),('results','md')]),(A.PRIMARY,[('scripts','py'),('results','json'),('results','md')])]
        for folder,ext in pairs}
    return dict(status='PASS_EXACT65_ESCAPE7_CAP',sources=sources,numerical=n,physical=p,
        scope=dict(full_two_corner_duty_domain=True,all_right_member_exact=True,other_two_members_numerical=True,
            compact7_section_numerical=True,full_ambient_facet7_compact=False,finite_volume_proved=False))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    out=enc(run())
    if args.write:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    else:assert json.loads(path.read_text())==out
    print('PASS',TAG,flush=True)


if __name__=='__main__':main()
