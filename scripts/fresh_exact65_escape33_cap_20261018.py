#!/usr/bin/env python3
"""Changed-anchor cap at30/33/66, preserving the exact65-wall chamber."""
import argparse,json
from collections import Counter
from itertools import combinations,product
from hashlib import sha256
import numpy as np
import fresh_exact64_escape7_cap_20261018 as F
A,D=F.A,F.D
ROOT=A.ROOT;TAG='fresh_exact65_escape33_cap_20261018';FRAME=[30,33,66];NEW=77
V,add,sub,mul,div,neg,sqrt,dot=A.V,A.add,A.sub,A.mul,A.div,A.neg,A.sqrt,A.dot
enc=A.enc;WHITE=A.WHITE

def read(stem):return json.loads((ROOT/f'results/{stem}.json').read_text())

def numerical():
    prior=read(F.TAG)
    for path,h in prior['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    s=D.T.Scout();s.ledgers=[];s.ledgerkeys={};s.hist=Counter();s.unresolved=[];s.counts=Counter();s.blockers=Counter()
    D.P.H.H.old_frame_interface(s)
    nd={int(i):np.array(v)for i,v in prior['numerical']['snapshot']['old_normals'].items()}
    nd[76]=np.array(prior['numerical']['recipes'][0]['branches'][0]['normal']);tt={75:0,76:0}
    domains=[]
    for frame,bonds in [([30,59,75],[5,3]),(FRAME,[2,4])]:
        assert [s.ps(i,tt)for i in frame]==[0,0,0]
        rows=[];rejected=[]
        for pi,ms in product(range(7),product(range(2,6),repeat=3)):
            duties=[D.T.P.corner3(0,0,pi,b,*xy)for b,xy in zip(bonds,(ms[:2],ms[1:]))]
            if any(x is None for x in duties):
                rejected.append(dict(profile=pi,labels=ms,duties=duties));continue
            rhs=[-D.T.B.theta(m)for m in ms]
            roots,info=D.T.N.norm_roots([nd[i]for i in frame],rhs,s.profiles[pi]['square'])
            assert len(roots)==2
            row=dict(profile=pi,labels=ms,duties=duties,norm=info,branches=[]);rows.append(row)
            for sign,v in roots:
                prods,bad,cp=s.native(v,pi,nd,tt)
                row['branches'].append(dict(branch=sign,normal=v.tolist(),products=sorted(prods.items()),
                    blockers=[[i,prods[i]]for i in bad],center_product=cp,native=not bad and cp<0))
        domains.append(dict(frame=frame,bonds=bonds,labels=[2,3,4,5],profiles=list(range(7)),recipes=rows,rejected_finite_duties=rejected))
    assert [len(x['recipes'])for x in domains]==[17,38]
    assert not any(b['native']for r in domains[0]['recipes']for b in r['branches'])
    cones=[]
    for row in domains[1]['recipes']:
        for b in row['branches']:
            if not b['native']:continue
            full={**nd,NEW:np.array(b['normal'])}
            cone=D.T.S.C.detailed_cone(full,list(full));assert not cone['warnings']
            faces={str(i):[r for r in cone['rays']if i in r['active']]for i in[7,33,60,75,76,NEW]}
            assert all(r['kind']=='finite'for i in('7','33','60')for r in faces[i])
            cones.append(dict(profile=row['profile'],labels=row['labels'],branch=b['branch'],normal=b['normal'],cone=cone,faces=faces))
            print('CONE',len(cones),row['profile'],row['labels'],len(cone['rays']),cone['open_ends'],flush=True)
    assert len(cones)==28
    return dict(domains=domains,accepted=cones,snapshot=dict(metric=D.T.M.tolist(),center=s.center.tolist(),
        old_normals={str(i):v.tolist()for i,v in nd.items()},profiles=s.profiles),
        summary=dict(eligible_recipes=[17,38],norm_roots=[34,76],native=[0,28],deleted_original_walls=[]))

def physical():
    d=read(A.TAG);p=read(A.PRIMARY);prior=read(F.TAG)
    for x in(d,p,prior):
        for path,h in x['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    z=d['physical_input'];metric=[list(map(V,r))for r in z['metric']];h=list(map(V,z['center']));h2=V(z['center_square'])
    roots={int(i):list(map(V,v))for i,v in z['roots'].items()};roots[75]=list(map(V,d['normal']));roots[76]=list(map(V,prior['physical']['normal']));ids=sorted(roots)
    hc={int(i):V(v)for i,v in z['center_products'].items()};hc[75]=V(d['center_product']);hc[76]=V(prior['physical']['center_product'])
    wi=[list(map(V,r))for r in d['white_inverse']];wc={int(i):list(map(V,v))for i,v in d['white_covectors'].items()};wc[75]=[V(0)]*6;wc[76]=[V(0)]*6
    hw=[hc[w]for w in WHITE];hp={i:sub(hc[i],dot(hw,wi,wc[i]))for i in ids};hp2=sub(h2,dot(hw,wi,hw))
    native={};identities=[]
    for i in FRAME:
        native[i]={}
        for j in ids:
            raw=dot(roots[i],metric,roots[j]);v,pr=A.refine(raw,i,j);native[i][j]=v
            if pr:identities.append(pr)
    assert next(r['identity']for r in identities if r['pair']==[30,33])==2
    assert next(r['identity']for r in identities if r['pair']==[33,66])==4
    assert native[30][66][1]<-2
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
    dotted=set(ids)-set(WHITE+FRAME);assert len(dotted)==56 and all(products[i][1]<-2 for i in dotted),[(i,enc(products[i]))for i in dotted if products[i][1]>=-2]
    norm=dot(v,metric,v);assert norm[0]<=2<=norm[1]
    oldw={i:sub(cp,mul(div(hc[i],V(2)),products[i]))for i in ids}
    neww={i:sub(hc[i],mul(div(cp,V(2)),products[i]))for i in ids}
    nsq=sub(h2,div(mul(cp,cp),V(2)));assert nsq[1]<0 and all(x[1]<0 for x in oldw.values())and all(x[1]<0 for x in neww.values())
    for i in FRAME:g[i][NEW]=products[i]
    g[NEW]={**products,NEW:V(2)};hp[NEW]=cp;pj[NEW]=v
    vertices=[]
    for support in([30,33,NEW],[33,66,NEW]):
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
        if 30 in subset and 66 in subset:trials.append(dict(support=support,dotted_pair=[30,66]))
        else:
            assert support in[sorted(r['support'])for r in vertices]
            trials.append(dict(support=support,finite=True))
    assert len(trials)==9 and sum('finite'in r for r in trials)==2
    # The prior strict escape need not be the region this changed anchor cuts.
    oldescape=prior['physical']['remaining_strict_escape'];oldprod=dot(v,metric,list(map(V,oldescape['point'])))
    # A fresh exact ray at the retained exterior30/52/76 region, followed by
    # a strict center perturbation. All gates below are physical inequalities.
    ef=[30,52,76];eg={}
    for i in ef:
        eg[i]={}
        for j in ids:
            raw=dot(roots[i],metric,roots[j]);val,proof=A.refine(raw,i,j)
            eg[i][j]=sub(val,dot(wc[i],wi,wc[j]))
    ei,edet=A.P.inverse([[eg[i][j]for j in ef]for i in ef]);assert edet[1]<0
    ec=[hp[i]for i in ef];eq=sub(hp2,dot(ec,ei,ec));assert eq[0]>0
    er=A.vector_sub(ph,A.linear(A.row_times(ec,ei),[pj[i]for i in ef]));er=[neg(x)for x in er]
    escape=None
    for power in range(1,50):
        epsilon=V(A.F(1,2**power));point=[add(x,mul(epsilon,y))for x,y in zip(er,h)]
        direct_sq=dot(point,metric,point);direct_future=dot(point,metric,h)
        sq=sub(add(eq,mul(mul(epsilon,epsilon),h2)),mul(V(2),mul(epsilon,eq)))
        future=add(neg(eq),mul(epsilon,h2))
        direct_pp={i:dot(point,metric,r)for i,r in{**roots,NEW:v}.items()}
        pp={}
        for i in ids+[NEW]:
            cols=[eg[j][i]for j in ef]if i!=NEW else[products[j]for j in ef]
            val=add(sub(dot(ec,ei,cols),hp[i]),mul(epsilon,hc[i]if i!=NEW else cp))
            if i in WHITE+ef:val=mul(epsilon,hc[i])
            pp[i]=(max(val[0],direct_pp[i][0]),min(val[1],direct_pp[i][1]));assert pp[i][0]<=pp[i][1]
        sq=(max(sq[0],direct_sq[0]),min(sq[1],direct_sq[1]));assert sq[0]<=sq[1]
        future=(max(future[0],direct_future[0]),min(future[1],direct_future[1]));assert future[0]<=future[1]
        if sq[0]>0 and future[1]<0 and all(x[1]<0 for x in pp.values()):
            escape=dict(frame=ef,power=power,point=point,square=sq,future=future,products=pp,complement_square=eq,direct_square=direct_sq,direct_future=direct_future,direct_products=direct_pp);break
    return dict(old_walls=ids,new_wall=NEW,deleted_original_walls=[],frame=FRAME,profile=0,labels=[2,2,2],
        identity_refinements=identities,frame_Gram=gram,frame_inverse=inv,frame_determinant=det,
        center_complement_square=q,factor=fac,normal=v,direct_norm=norm,center_product=cp,
        schur_products=schur,direct_products=direct,native_products=products,
        old_facet_witness_products=oldw,new_facet_witness_products=neww,new_facet_witness_square=nsq,
        new_finite_vertices=vertices,new_finite_support_trials=trials,previous_escape_new_product=oldprod,
        remaining_strict_escape=escape,
        summary=dict(total_walls=66,old_finite_vertices_retained=203,new_finite_vertices=2,complete_finite_vertices=205,
            strict_dotted_contacts=56,right_contacts=9,infinite_volume_proved=escape is not None))


def run():
    n=numerical();p=physical();assert p['summary']['infinite_volume_proved']
    sources={f'{folder}/{stem}.{ext}':sha256((ROOT/f'{folder}/{stem}.{ext}').read_bytes()).hexdigest()
        for stem,pairs in[(TAG,[('scripts','py')]),(F.TAG,[('scripts','py'),('results','json'),('results','md')]),
            (A.TAG,[('scripts','py'),('results','json'),('results','md')]),(D.TAG,[('scripts','py'),('results','json'),('results','md')]),
            (A.PRIMARY,[('scripts','py'),('results','json'),('results','md')])]for folder,ext in pairs}
    return dict(status='PASS_EXACT66_ESCAPE33_CAP',sources=sources,numerical=n,physical=p,
        scope=dict(first_two_corner_duty_domain_complete=True,second_low_label_domain_complete=True,
            second_high_labels_unsearched=True,all_right_member_exact=True,other_members_numerical=True,
            compact33_section_numerical=True,full_ambient_facet33_compact=False,finite_volume_proved=False))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    out=enc(run())
    if args.write:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    else:assert json.loads(path.read_text())==out
    print('PASS',TAG,flush=True)

if __name__=='__main__':main()
