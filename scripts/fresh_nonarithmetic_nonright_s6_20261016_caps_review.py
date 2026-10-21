#!/usr/bin/env python3
"""Independent native-Gram reconstruction of the two nonright-S6 cutters.

No producer arithmetic imports. Strict signs use independent rational intervals;
projection definitions, not interval overlaps, establish exact equalities.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import fresh_nonarithmetic_escape52_20261016_cap_review as C

B,ROOT,I=C.B,C.ROOT,C.I
TAG='fresh_nonarithmetic_nonright_s6_20261016_caps_review'
PRIMARY='results/fresh_nonarithmetic_nonright_s6_20261016_caps.json'
OLD_CONE='results/fresh_nonarithmetic_escape52_20261016_cone.json'
S=[0,1,2,3,6,7,8,10]
PROFILES=[(34,{6:5,8:3,34:5}),(12,{3:3,10:10,12:4})]


def overlap(x,raw):
    lo,hi=map(F,raw);assert x.lo<=hi and lo<=x.hi


def positive_pivots(matrix):
    a=[row[:] for row in matrix];p=[]
    for k in range(len(a)):
        q=a[k][k];p.append(q.saved())
        if q.hi<0:return False,p
        assert q.lo>0,('unresolved pivot',k,q.saved())
        for i in range(k+1,len(a)):
            for j in range(i,len(a)):
                a[i][j]=a[j][i]=a[i][j]-a[i][k]*a[k][j]/q
    return True,p


def geometry():
    g=C.geometry();primary=B.read(PRIMARY);oldcone=B.read(OLD_CONE)
    cone_review=B.read('results/fresh_nonarithmetic_escape52_20261016_cone_review.json')
    assert cone_review['sources'][OLD_CONE]==B.digest(OLD_CONE)
    assert primary['status']=='PASS_EXACT_TWO_NONRIGHT_S6_ESCAPE_CUTTERS'
    norms=[B.iv(B.dot(n,c)) for n,c in zip(g['N'],g['cov'])]+[I(2)]*4
    assert len(norms)==len(g['NI'])==53
    assert B.dot(g['N'][12],g['cov'][12])==F(1,4)
    E=[B.iv(B.dec(x)) for x in oldcone['strict_escape']['point']]
    ec=B.imv(g['QI'],E);hh=B.idot(g['HI'],g['hI'])
    theta={2:I(),3:I(1),4:I(2).sqrt(),5:(I(1)+I(5).sqrt())/2,
           6:I(3).sqrt(),10:((I(5)+I(5).sqrt())/2).sqrt()}
    records=[];newnormals=[];censuses=[];affinechecks=[]
    for new,(anchor,profile) in enumerate(PROFILES,53):
        old=primary['caps'][new-53];frame=S+[anchor]
        assert old['index']==new and old['frame']==frame
        assert dict(old['finite_contacts'])=={i:profile.get(i,2) for i in frame}
        # Reverse elimination differs from the producer's forward inversion.
        native=[[B.iv(B.dot(g['N'][i],g['cov'][j])) for j in frame] for i in frame]
        inv,pivots=C.inverse_interval(native)
        contacts=[-(norms[i]/2).sqrt()*theta[profile.get(i,2)] for i in frame]
        coeff=B.imv(inv,contacts);particular=C.combine(coeff,[g['NI'][i] for i in frame])
        hs=[B.idot(g['HI'],g['covI'][i]) for i in frame]
        hcoeff=B.imv(inv,hs);hproj=C.combine(hcoeff,[g['NI'][i] for i in frame])
        kernel=[x-y for x,y in zip(g['HI'],hproj)]
        qk=hh-B.idot(hs,hcoeff);pp=B.idot(contacts,coeff)
        assert qk.lo>0
        delta=(2-pp)/qk;assert delta.lo>0
        branches=[];accepted=[]
        for sign in (-1,1):
            normal=[x+sign*delta.sqrt()*y for x,y in zip(particular,kernel)]
            hp=B.idot(normal,g['hI']);cut=B.idot(normal,ec)
            off=[B.idot(normal,g['covI'][j]) for j in range(53) if j not in frame]
            strict=all(x.hi<0 for x in off)
            branches.append(dict(sign=sign,center_pair=hp.saved(),escape_pair=cut.saved(),strict_off_frame_acute=strict))
            if hp.hi<0<cut.lo and strict:accepted.append((sign,normal,hp,cut))
        assert len(accepted)==1
        sign,normal,hp,cut=accepted[0];assert sign==old['selected_sign']==-1
        for x,raw in zip(normal,old['normal_interval']):overlap(x,raw)
        for x,raw in zip(contacts,old['native_contact_vector']):overlap(x,raw)
        if anchor==12:
            assert contacts[-1].lo<=F(-1,2)<=contacts[-1].hi
            assert old['finite_contacts'][-1]==[12,4]
        relations=[]
        for j,(cov,q) in enumerate(zip(g['covI'],norms)):
            value=B.idot(normal,cov)
            if j in frame:
                label=profile.get(j,2);kind='right' if label==2 else 'finite'
                expected=contacts[frame.index(j)]
                assert value.lo<=expected.hi and expected.lo<=value.hi
            else:
                kind='dotted';label=None
                assert value.hi<0 and value.hi*value.hi>2*q.hi,(new,j,value.saved())
            raw=old['relations'][j];assert raw['kind']==kind
            if label is not None:assert label==raw['label']
            overlap(value,raw['pairing_interval'])
            relations.append(dict(wall=j,kind=kind,label=label,bounds=value.saved()))
        g['NI'].append(normal);g['covI'].append(B.imv(g['QI'],normal));norms.append(I(2));newnormals.append(normal)
        # Every candidate finite vertex containing the new root needs eight
        # of its exactly nine finite neighbors. A dotted pair or a certified
        # negative Schur pivot excludes each candidate, without symbolic
        # determinants or copied primary inertia classifications.
        for support in combinations(frame,8):
            ids=list(support)+[new];witness=None
            for a,b in combinations(support,2):
                value=B.idot(g['NI'][a],g['covI'][b])
                if value.hi<0 and value.hi*value.hi>norms[a].hi*norms[b].hi:
                    witness=dict(dotted_pair=[a,b],bounds=value.saved());break
            if witness is None:
                gram=[[norms[a] if a==b else B.idot(g['NI'][a],g['covI'][b]) for b in ids] for a in ids]
                pd,lp=positive_pivots(gram);assert not pd
                witness=dict(negative_schur_pivot=lp[-1],preceding_positive_pivots=lp[:-1])
            censuses.append(dict(support=ids,**witness))
        # Enumerate every connected standard-label subset containing new.
        # An irreducible Euclidean Coxeter diagram uses only labels3,4,6
        # (or the parallel rank-two diagram, already excluded by strict gates).
        def standard_label(a,b):
            if a==b:return 2
            if new in (a,b):return profile.get(b if a==new else a,2)
            value=B.dot(g['N'][a],g['cov'][b])
            if value==0:return 2
            assert B.sg(value)<0
            square=value*value;product=B.dot(g['N'][a],g['cov'][a])*B.dot(g['N'][b],g['cov'][b])
            for label,factor in ((3,F(1,4)),(4,F(1,2)),(6,F(3,4))):
                if square==product*factor:return label
            return 0
        standard=[i for i in frame if profile.get(i,2) in (2,3,4,6)]
        checks=[]
        for size in range(len(standard)+1):
            for sub in combinations(standard,size):
                ids=list(sub)+[new];labels={(a,b):standard_label(a,b) for a,b in combinations(ids,2)}
                if any(m not in (2,3,4,6) for m in labels.values()):continue
                reached={new}
                while True:
                    nxt=reached|{b for (a,b),m in labels.items() if m!=2 and a in reached}|{a for (a,b),m in labels.items() if m!=2 and b in reached}
                    if nxt==reached:break
                    reached=nxt
                if len(reached)!=len(ids):continue
                gram=[[norms[a] if a==b else B.idot(g['NI'][a],g['covI'][b]) for b in ids] for a in ids]
                pd,ps=positive_pivots(gram);assert pd
                checks.append(dict(support=ids,positive_pivots=ps))
        affinechecks.append(dict(new_wall=new,all_connected_standard_subsets=checks))
        records.append(dict(index=new,frame=frame,native_contact_vector=[x.saved() for x in contacts],
            inverse_pivots=pivots,positive_kernel_norm=qk.saved(),particular_norm=pp.saved(),
            positive_radical_square=delta.saved(),branches=branches,selected_sign=sign,
            normal=[x.saved() for x in normal],relations=relations,center_pair=hp.saved(),escape_pair=cut.saved()))
    assert len(censuses)==18
    primary_supports={tuple(x['support']) for x in primary['finite_support_census']}
    assert {tuple(x['support']) for x in censuses}==primary_supports
    assert not primary['new_finite_vertices']
    mutual=B.idot(newnormals[0],g['covI'][54]);assert mutual.hi<-2
    facets=[];centerpairs=[]
    for n,q in zip(g['NI'],norms):
        hp=B.idot(n,g['hI']);assert hp.hi<0
        fq=hh-hp*hp/q;assert fq.hi<0
        centerpairs.append(hp.saved());facets.append(fq.saved())
    g.update(two_cap_records=records,two_cap_censuses=censuses,two_cap_affine=affinechecks,
        all55_center_pairs=centerpairs,all55_facet_squares=facets,mutual=mutual.saved())
    return g


def run():
    g=geometry()
    paths=[f'scripts/{TAG}.py',PRIMARY,OLD_CONE,f'scripts/{C.TAG}.py',f'results/{C.TAG}.json',
        'results/fresh_nonarithmetic_escape52_20261016_cone_review.json',
        'results/fresh_nonarithmetic_affine_filtered_20261016_affine_gate.json']
    return dict(status='PASS_INDEPENDENT_TWO_NONRIGHT_S6_CUTTERS_AND_ZERO_NEW_FINITE_VERTICES',
        sources={p:B.digest(p) for p in paths},caps=g['two_cap_records'],mutual_dotted_pair=g['mutual'],
        actual_facets=55,all55_center_pairs=g['all55_center_pairs'],facet_projection_squares=g['all55_facet_squares'],
        complete_new_finite_vertex_count=0,all18_rejected_supports=g['two_cap_censuses'],
        no_new_affine_checks=g['two_cap_affine'],
        proof='Independent53geometry. Two reverse native-Gram inverses define exact projection/radical roots. All107new contacts, strict common center, both escape cuts and mutual dottedness are checked. Native old12 norm1/4 gives contact-1/2. Every18possible finite support has a dotted pair or negative Schur pivot. All connected standard-label new subdiagrams are positive definite, so no new affine circuit occurs.',
        scope='Exact local55wall geometry and complete new-facet finite census only. Ambient cone and renewed escape are separate; no finite-volume or new-dimension claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
