#!/usr/bin/env python3
"""Independent53rd cap from reviewed52 geometry and literal nine-root Gram."""
import argparse
import json
from fractions import Fraction as F
from itertools import combinations
import fresh_nonarithmetic_affine_filtered_20261016_cluster_review as C

B,ROOT,I=C.B,C.ROOT,C.I
TAG='fresh_nonarithmetic_escape52_20261016_cap_review'
PRIMARY='results/fresh_nonarithmetic_escape52_20261016_cap.json'
PRIOR='results/fresh_nonarithmetic_affine_filtered_20261016_cone.json'
SUPPORT=[0,1,2,3,5,6,8,34]
FRAME=SUPPORT+[30]


def inverse_interval(matrix):
    n=len(matrix)
    a=[row[:]+[I(int(i==j)) for j in range(n)] for i,row in enumerate(matrix)]
    pivots=[]
    for k in reversed(range(n)):
        row=next(i for i in range(k+1) if a[i][k].lo>0 or a[i][k].hi<0)
        a[k],a[row]=a[row],a[k]
        q=a[k][k];pivots.append(dict(column=k,row=row,pivot=q.saved()))
        a[k]=[x/q for x in a[k]]
        for i in range(n):
            if i!=k:
                f=a[i][k];a[i]=[x-f*y for x,y in zip(a[i],a[k])]
    return [row[n:] for row in a],pivots


def combine(c,vs):return [B.idot(c,[v[j] for v in vs]) for j in range(10)]


def geometry():
    g=C.geometry();target=B.read(PRIMARY);prior=B.read(PRIOR)
    assert target['right_support']==SUPPORT and target['frame']==FRAME
    assert target['finite_anchor']==dict(wall=30,label=4)
    G=[[B.dot(g['N'][i],g['cov'][j]) for j in FRAME] for i in FRAME]
    assert all(G[i][i]==2 for i in range(9))
    GI=[[B.iv(x) for x in row] for row in G]
    J,pivots=inverse_interval(GI)
    c=[I()]*8+[-I(2).sqrt()]
    weights=B.imv(J,c);p=combine(weights,[g['NI'][i] for i in FRAME])
    hc=[B.idot(g['HI'],g['covI'][i]) for i in FRAME]
    hweights=B.imv(J,hc)
    hprojection=combine(hweights,[g['NI'][i] for i in FRAME])
    kernel=[x-y for x,y in zip(g['HI'],hprojection)]
    hh=B.idot(g['HI'],g['hI'])
    qk=hh-B.idot(hc,hweights);assert qk.lo>0
    E=[B.iv(B.dec(x)) for x in prior['strict_escape']['point']]
    ec=B.imv(g['QI'],E)
    orientation=B.idot(kernel,ec).strict_sign()
    if orientation<0:kernel=[-x for x in kernel]
    pp=B.idot(c,weights);delta=(2-pp)/qk;assert delta.lo>0
    normal=[x+delta.sqrt()*y for x,y in zip(p,kernel)]
    covector=B.imv(g['QI'],normal)
    hp,cut=B.idot(normal,g['hI']),B.idot(normal,ec)
    assert hp.hi<0<cut.lo
    oldnorms=[B.iv(B.dot(n,cov)) for n,cov in zip(g['N'],g['cov'])]+[I(2)]*3
    contacts=[]
    for i,(cov,norm) in enumerate(zip(g['covI'],oldnorms)):
        value=B.idot(normal,cov)
        if i in FRAME:
            label=4 if i==30 else 2
            expected=-I(2).sqrt() if i==30 else I()
            assert value.lo<=expected.hi and expected.lo<=value.hi
            kind='finite' if i==30 else 'right'
        else:
            label=None;kind='dotted'
            assert value.hi<0 and value.hi*value.hi>2*norm.hi,(i,value.saved())
        assert target['relations'][i]['kind']==kind
        if label is not None:assert target['relations'][i]['label']==label
        contacts.append(dict(wall=i,kind=kind,label=label,bounds=value.saved()))
    for value,raw in zip(normal,target['normal_interval']):
        lo,hi=map(F,raw);assert value.lo<=hi and lo<=value.hi
    # The8-root support is PD. Orthogonality of the new root adds A1.
    sgram=[row[:8] for row in GI[:8]]
    schur=[row[:] for row in sgram];pd=[]
    for k in range(8):
        pivot=schur[k][k];assert pivot.lo>0;pd.append(pivot.saved())
        for i in range(k+1,8):
            for j in range(i,8):
                schur[i][j]=schur[j][i]=schur[i][j]-schur[i][k]*schur[k][j]/pivot
    JS,_=inverse_interval(sgram)
    scoeff=B.imv(JS,hc[:8]);sproj=combine(scoeff,[g['NI'][i] for i in SUPPORT])
    vertex=[x-y-hp*z/2 for x,y,z in zip(g['HI'],sproj,normal)]
    vnorm=B.idot(vertex,B.imv(g['QI'],vertex));vtime=B.idot(vertex,g['hI'])
    assert vnorm.hi<0 and vtime.hi<0
    for i,cov in enumerate(g['covI']+[covector]):
        value=B.idot(vertex,cov)
        if i in SUPPORT+[52]:assert value.lo<=0<=value.hi
        else:assert value.hi<0,(i,value.saved())
    rejected=[]
    for support in combinations(FRAME,8):
        if list(support)==SUPPORT:continue
        witnesses=[]
        for i,j in combinations(support,2):
            value=B.dot(g['N'][i],g['cov'][j])
            if B.sg(value)<0 and B.sg(value*value-4)>0:witnesses.append([i,j])
        assert witnesses,support
        rejected.append(dict(support=sorted(list(support)+[52]),dotted_pair=witnesses[0]))
    assert len(rejected)==8
    assert [row['support'] for row in target['new_finite_vertices']]==[SUPPORT+[52]]
    full=g['NI']+[normal];norms=oldnorms+[I(2)]
    facets=[]
    for n,nq in zip(full,norms):
        h=B.idot(n,g['hI']);assert h.hi<0
        value=hh-h*h/nq;assert value.hi<0
        facets.append(value.saved())
    # No connected finite-label path continues from new52 beyond old30:
    # every other finite neighbor is right or dotted to30.
    for i in SUPPORT:
        value=B.dot(g['N'][30],g['cov'][i])
        assert value==0 or (B.sg(value)<0 and B.sg(value*value-4)>0)
    g['NI']=full;g['covI'].append(covector)
    g.update(new_cap=normal,new_contacts=contacts,new_facet_squares=facets,
             new_radical=delta,new_kernel_norm=qk,new_inverse_pivots=pivots,
             new_cut=cut,new_center_pair=hp,new_support_pivots=pd,
             new_vertex=dict(support=SUPPORT+[52],point=[x.saved() for x in vertex],
                             norm=vnorm.saved(),center_pair=vtime.saved()),
             rejected_finite_supports=rejected)
    return g


def run():
    g=geometry()
    paths=[f'scripts/{TAG}.py',PRIMARY,PRIOR,f'scripts/{C.TAG}.py',f'results/{C.TAG}.json',
           'results/fresh_nonarithmetic_affine_filtered_20261016_cone_review.json',C.PRIMARY,
           'results/fresh_nonarithmetic_affine_filtered_20261016_affine_gate.json']
    return dict(status='PASS_INDEPENDENT_53RD_CAP_AND_COMPLETE_NEW_FINITE_VERTEX_CENSUS',
                sources={p:B.digest(p) for p in paths},norm_exact=2,normal=[x.saved() for x in g['new_cap']],
                contacts=g['new_contacts'],actual_facets=53,facet_projection_squares=g['new_facet_squares'],
                inverse_pivots=g['new_inverse_pivots'],kernel_norm=g['new_kernel_norm'].saved(),
                positive_radical_square=g['new_radical'].saved(),escape_cut=g['new_cut'].saved(),
                center_pair=g['new_center_pair'].saved(),support_pivots=g['new_support_pivots'],
                complete_new_finite_vertices=[g['new_vertex']],rejected_finite_supports=g['rejected_finite_supports'],
                no_new_affine_circuit=True,
                proof='Independent52physical geometry; reverse interval inversion of literal9-root Gram defines exact projection and positive-radical normal. Formula identities prove norm2,8rights andpi4to30; strict intervals prove43dots. PD8support plusorthogonalnewroot givesonefinitecorner; allother8possible supports containa dottedpair. Commonstrictcenter proves53actualfacets. Newconnected standardlabel components areonlyA1orI2(4), so noaffine circuit isintroduced.',
                scope='Exact local53-wall reconstruction and complete finitevertex census onnewfacet. Ambientcone andvolume are separate. Old52affine absence is inherited; nofinitevolume/newdimension claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
