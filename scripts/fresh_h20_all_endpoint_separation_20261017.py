#!/usr/bin/env python3
"""All finite endpoint profiles: low singleton angles and continuous high tail.

Inherits the frozen fixed/fixed gate. No parabolic endpoint census is claimed.
"""
import argparse
import hashlib
import json
import math
from fractions import Fraction as F
import sympy as s
import fresh_h20_e8_exact_cluster_20261017 as P
import fresh_h20_e8_endpoint_belt_20261017_exact_gate as E

ROOT=P.ROOT
TAG='fresh_h20_all_endpoint_separation_20261017'
END='results/fresh_h20_e8_partner_endpoints_20261017.json'
PRIOR='results/fresh_h20_e8_endpoint_belt_20261017_exact_gate.json'
pt,add,sub,neg,mul,div,sqrt=P.pt,P.add,P.sub,P.neg,P.mul,P.div,P.sqrt
scale,plus,sx=E.scale,E.plus,E.sx

def run():
    C=json.loads((ROOT/P.CUSP).read_text())
    D=json.loads((ROOT/END).read_text())
    prior=json.loads((ROOT/PRIOR).read_text())
    assert prior['status']=='PASS_EXACT_ALL158LABEL_FIXED_ENDPOINT_SEPARATION'
    rejected={(r['edge'],r['profile'],r['label158']) for r in prior['rejected_discretes']}
    ids=C['core_ids'];G=s.Matrix(C['core_gram']);GI=G.inv()
    IG=[[P.alg(x) for x in row] for row in G.tolist()]
    def gp(v):return P.mv(IG,v)
    Z=json.loads((ROOT/'results/fresh_h20_e8_exact_cluster_20261017.json').read_text())
    H=json.loads((ROOT/P.GEOM).read_text())
    Q=[list(map(F,row)) for row in H['form']]
    nz=[(i,j,q) for i,row in enumerate(Q) for j,q in enumerate(row) if q]
    N=[[tuple(map(F,z)) for z in v] for v in H['normals']]+[[tuple(map(F,z)) for z in r['physical_normal']] for r in Z['steps']]
    qs=list(map(F,H['norms']))+[F(2)]*4
    thresholds=[sqrt(pt(2*q)) for q in qs]
    def physical(v,w):return P.total(mul(pt(q),mul(v[i],w[j])) for i,j,q in nz)
    BP=[[physical(N[i],n) for n in N] for i in ids]
    hp=[physical(N[i],[tuple(map(F,z)) for z in Z['strict_center']]) for i in ids]
    xlabels=[(2,s.Integer(0)),(3,s.Integer(1)),(4,s.sqrt(2)),(5,(1+s.sqrt(5))/2),(6,s.sqrt(3))]
    fixed=[];low=[];high=[];low_rejected=[];corners=[]
    for ci,cen in enumerate(D['endpoint_censuses']):
        assert 151 in cen['support']
        wi=[ids.index(i) for i in cen['support']]
        SB=s.eye(21)[:,wi]
        WJ=s.Matrix(cen['inverse'])
        ix=ids.index(158);ee=s.eye(21)[:,ix]
        u=(ee-SB*WJ*SB.T*G*ee).applyfunc(s.simplify)
        sigma=s.simplify((u.T*G*u)[0]);assert sigma.is_positive
        unit=sx((u/s.sqrt(sigma)).applyfunc(s.simplify))
        missing=next(i for i in (159,160,161) if i not in cen['edge'])
        mi=ids.index(missing)
        V=sx((-GI[:,mi]/s.sqrt(-GI[mi,mi])).applyfunc(s.simplify))
        assert P.inner(V,hp)[1]<0
        j151=cen['support'].index(151)
        assert WJ[j151,j151]==s.Rational(1,2)
        assert all(WJ[j151,j]==0 for j in range(19) if j!=j151)
        corners.append(dict(edge=ci,support=cen['support'],sigma=str(sigma),
            isolated151_inverse_diagonal='1/2',future_vertex_pair=P.save(P.inner(V,hp))))
        def make(contacts,name,threshold,kind):
            dd=s.zeros(19,1)
            for j,label in contacts:dd[cen['support'].index(j)]=2*s.cos(s.pi/label)
            pp=(-SB*WJ*dd).applyfunc(s.simplify)
            h=s.simplify(2-(pp.T*G*pp)[0]);assert h.is_positive
            off=s.simplify((pp.T*G*ee)[0]);p=sx(pp)
            a=scale(plus(V,scale(unit,pt(-1))),div(sqrt(P.alg(h)),pt(2)))
            b=scale(plus(scale(V,pt(-1)),scale(unit,pt(-1))),div(sqrt(P.alg(h)),pt(2)))
            t=div(P.alg(threshold+off),sqrt(P.alg(sigma*h)));assert t[0]>1
            k=add(t,sqrt(sub(mul(t,t),pt(1))))
            aa=scale(a,k);bb=scale(b,div(pt(1),k))
            nums=[plus(plus(p,aa),bb),plus(p,scale(aa,pt(2))),aa]
            fam=dict(ci=ci,name=name,kind=kind,contacts=contacts,nums=nums,gnums=[gp(v) for v in nums],discrete=[])
            if kind=='high':return fam
            for label,x in xlabels:
                rad=s.simplify((x+off)**2/sigma-h)
                if kind=='fixed' and (ci,name,label) in rejected:continue
                ri=P.alg(rad)
                if ri[1]<=0 or rad==0:
                    assert kind=='low'
                    low_rejected.append(dict(edge=ci,singleton_label=name,label158=label,reason='nonpositive norm radicand',radicand=str(rad)))
                    continue
                assert ri[0]>0
                r=plus(plus(p,scale(unit,neg(div(P.alg(x+off),sqrt(P.alg(sigma)))))),scale(V,sqrt(ri)))
                if kind=='low':
                    fail=None
                    for j in range(162):
                        if j in cen['support'] or j==158:continue
                        g=P.inner(r,[BP[z][j] for z in range(21)])
                        test=add(g,thresholds[j])
                        if test[0]>0:
                            fail=dict(wall=j,pair=P.save(g),dotted_test=P.save(test));break
                        assert test[1]<0,(ci,name,label,j,P.save(test))
                    if fail:
                        low_rejected.append(dict(edge=ci,singleton_label=name,label158=label,reason='old unprescribed contact not dotted',witness=fail));continue
                fam['discrete'].append(dict(label=label,r=r,gr=gp(r)))
            return fam
        for row in cen['profiles']:
            fixed.append(make(row['contacts'],row['id'],s.Rational(9,5),'fixed'))
        for m in (3,4,5):low.append(make([[151,m]],m,s.Rational(9,5),'low'))
        high.append(make([[151,6]],6,s.Integer(0),'high'))
        print('EDGE',ci,'fixed',len(cen['profiles']),'low3','high1',flush=True)
    assert [sum(len(f['discrete']) for f in fixed if f['ci']==i) for i in range(3)]==prior['discrete_counts']
    stats={};maxima={};exceptions=[];gaps=[];stream=hashlib.sha256()
    def coefficients(kind,vals,rec):
        upper=max(x[1] for x in vals)
        stats[kind]=stats.get(kind,0)+1
        maxima[kind]=max(maxima.get(kind,F(-10**99)),upper)
        stream.update(json.dumps([kind,rec,[P.save(z) for z in vals]],sort_keys=True).encode())
        if upper>=0:exceptions.append(dict(kind=kind,record=rec,coefficients=[P.save(z) for z in vals]))
    # High singleton labels: actual finite neighbors must be right to151 and
    # to the cap. The m6 conic with x158>=0 bounds all d>=sqrt3 after scaling.
    for h in high:
        for f in fixed:
            if h['ci']==f['ci'] or any(j==151 for j,label in f['contacts']):continue
            rec=[h['ci'],f['ci'],f['name']]
            coefficients('high_tail', [P.inner(v,w) for v in h['nums'] for w in f['gnums']],rec)
            for r in f['discrete']:
                coefficients('high_discrete',[P.inner(v,r['gr']) for v in h['nums']],rec+[r['label']])
        print('HIGH_EDGE',h['ci'],'exceptions',len(exceptions),flush=True)
    # Low singleton labels3,4,5: full pair separation, with every158 label.
    for ni,f in enumerate(low):
        for g in fixed+low[:ni]:
            if f['ci']==g['ci']:continue
            rec=[f['ci'],f['name'],g['ci'],g['kind'],g['name']]
            coefficients('low_tail_tail',[add(P.inner(v,w),pt(2 if a<2 and b<2 else 0)) for a,v in enumerate(f['nums']) for b,w in enumerate(g['gnums'])],rec)
            for r in f['discrete']:
                coefficients('low_discrete_tail',[add(P.inner(r['r'],w),pt(2 if b<2 else 0)) for b,w in enumerate(g['gnums'])],rec+[r['label']])
            for r in g['discrete']:
                coefficients('low_tail_discrete',[add(P.inner(v,r['gr']),pt(2 if a<2 else 0)) for a,v in enumerate(f['nums'])],rec+[r['label']])
            for r in f['discrete']:
                for z in g['discrete']:
                    dot=P.inner(r['r'],z['gr']);item=rec+[r['label'],z['label']]
                    stats['low_discrete_discrete']=stats.get('low_discrete_discrete',0)+1
                    if dot[1]<-2:continue
                    if dot[0]>0:
                        gaps.append(dict(record=item,type='obtuse',pair=P.save(dot)));continue
                    co=mul(pt(F(-1,2)),dot)
                    assert 0<co[0]<co[1]<1,(item,P.save(co))
                    lab=max(2,int(math.pi/math.acos(float(sum(co)/2))))
                    lo,hi=E.cosine(lab),E.cosine(lab+1)
                    assert lo[1]<co[0] and co[1]<hi[0],(item,lab,P.save(co))
                    gaps.append(dict(record=item,type='strict finite-label gap',labels=[lab,lab+1],cosine=P.save(co),lower=P.save(lo),upper=P.save(hi)))
        print('LOW',f['ci'],f['name'],'exceptions',len(exceptions),'gaps',len(gaps),flush=True)
    assert not exceptions,exceptions[:1]
    # A norm-one3-node Gram with one bond>=cos(pi/6) and any attached bond
    # >=1/2 has determinant<=0 (the omitted third bond only decreases it).
    assert 1-s.Rational(3,4)-s.Rational(1,4)==0
    assert E.cosine(7)[0]>F(9,10)
    paths=[f'scripts/{TAG}.py',END,PRIOR,P.CUSP,P.GEOM,
           'results/fresh_h20_e8_exact_cluster_20261017.json',
           'results/fresh_h20_e8_f4_cusps_20261017_height_gate.json',
           f'scripts/{P.TAG}.py',f'scripts/{E.TAG}.py']
    return dict(status='PASS_EXACT_ALL_FINITE_ENDPOINT_PROFILE_NEIGHBOR_EXCLUSION',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        inherited_fixed_profiles=[57,45,57],low_singleton_profiles=9,
        low_discrete_counts=[sum(len(f['discrete']) for f in low if f['ci']==i) for i in range(3)],
        low_discrete_rejections=low_rejected,
        high_right151_fixed_profile_counts=[sum(not any(j==151 for j,m in f['contacts']) for f in fixed if f['ci']==i) for i in range(3)],
        high_domain='Every integer151 label m>=6, by d>=sqrt3 and the continuous m6 conic x158>=0. Any actualT-neighbor must be right to both151 and the cap.',
        corner_checks=corners,case_counts=stats,
        strict_coefficient_maxima={k:str(v) for k,v in maxima.items()},
        coefficient_stream_sha256=stream.hexdigest(),low_discrete_gaps=gaps,exceptions=exceptions,
        scope='All spherical rank20 endpoint profiles W+r, with all unprescribed old162 contacts except158 strictly dotted, all158 Coxeter values, finite coherent batches retaining168. Distinct caps cannot supply actualT-neighbor sides. Combined with the stated actual-face argument this excludes finite-volume closure using only these finite endpoint caps. Parabolic/ideal endpoint profiles, additional finite old contacts, and other intermediate facets remain outside.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
