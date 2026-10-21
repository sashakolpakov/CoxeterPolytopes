#!/usr/bin/env python3
"""Exact norm2 extension closing S16+{151,156} to a compact quadrilateral.

The frozen independent predecessor is replayed in full through its geometric
assertions. Its local geometry is returned before packaging the old report;
the predecessor source is never modified. New assertions are below.
"""
import argparse
import ast
import hashlib
import json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_compact_quad_cap_20261016'
PREVIOUS='scripts/fresh_h20_release144_cap_20261016_review.py'

def geometry():
    path=ROOT/PREVIOUS
    tree=ast.parse(path.read_text())
    run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run')
    stop=max(i for i,n in enumerate(run.body) if isinstance(n,ast.Assign)
             and any(isinstance(t,ast.Name) and t.id=='paths' for t in n.targets))
    run.body=run.body[:stop]+ast.parse('def f():\n return locals()').body[0].body
    ns={'__name__':'_frozen_predecessor_geometry','__file__':str(path)}
    exec(compile(ast.fix_missing_locations(tree),str(path),'exec'),ns)
    return ns,ns['run']()

def run():
    ns,L=geometry()
    pt,plus,minus,neg,times,div,root=[ns[k] for k in ['pt','plus','minus','neg','times','div','root']]
    saved,vec,R=[ns[k] for k in ['saved','vec','R']]
    ipair,ilinear,ivadd,ivscale,matvec,inverse=[L[k] for k in ['ipair','ilinear','ivadd','ivscale','matvec','interval_inverse']]
    N,qs,H,Hres,S16=[L[k] for k in ['full157','norms157','H','Hres','S16']]
    ids=[95,145,151,156]
    a=ipair(N[95],N[145]);assert a[1]<-2
    gram=[[pt(2) if i==j else a if {i,j}=={95,145}
           else pt(-1) if {i,j}=={145,156} else pt(0) for j in ids] for i in ids]
    inv,pivots=inverse(gram)
    target=[neg(root(pt(3))),pt(-1),pt(0),pt(0)]
    p=ilinear([N[j] for j in ids],matvec(inv,target))
    k=ivadd(Hres,ivscale(ilinear([N[j] for j in ids],matvec(inv,[ipair(Hres,N[j]) for j in ids])),pt(-1)))
    p2=ns['inner'](target,matvec(inv,target));k2=ipair(k,k)
    assert k2[0]>0
    rad=div(minus(pt(2),p2),k2);assert rad[0]>0
    branches=[];accepted=[]
    rights=S16|{151,156}
    for sign in [-1,1]:
        r=ivadd(p,ivscale(k,times(pt(sign),root(rad))))
        vals=[ipair(r,n) for n in N]
        bad=[j for j,v in enumerate(vals) if j not in rights|{95,145}
             and not(v[1]<0 and v[1]**2>2*qs[j])]
        hc=ipair(r,H)
        branches.append(dict(sign=sign,failed_dotted_gates=bad,center_pair=saved(hc)))
        if not bad and hc[1]<0:accepted.append((sign,r,vals))
    assert len(accepted)==1
    sign,r,vals=accepted[0]
    assert R.overlap(ipair(r,r),pt(2))
    contacts=[]
    for j,value in enumerate(vals):
        if j in rights:kind='right';expected=pt(0)
        elif j==95:kind='pi/6';expected=neg(root(pt(3)))
        elif j==145:kind='pi/3';expected=pt(-1)
        else:kind='dotted';expected=None
        if expected is not None:assert R.overlap(value,expected)
        contacts.append(dict(wall=j,relation=kind,pair_bounds=saved(value)))
    counts=Counter(c['relation'] for c in contacts)
    assert counts==dict(right=18,**{'pi/6':1,'pi/3':1,'dotted':137})
    full=N+[r];norms=qs+[F(2)]
    hp=[ipair(H,n) for n in full];assert all(x[1]<0 for x in hp)
    facets=[minus(ipair(H,H),div(times(x,x),pt(q))) for x,q in zip(hp,norms)]
    assert all(x[1]<0 for x in facets)
    cusps=[ipair(r,v) for v in L['cusp_vectors']+[L['nu2']]]
    assert len(cusps)==5 and all(x[1]<0 for x in cusps)
    corners=[]
    for sidepair in [[95,142],[142,145],[145,157],[157,95]]:
        ci=[151,156,*sidepair]
        cg=[[pt(norms[i]) if i==j else ipair(full[i],full[j]) for j in ci] for i in ci]
        # Fill exact construction zeros to keep the proof identities explicit.
        for i in range(4):
            for j in range(4):
                pair={ci[i],ci[j]}
                if i!=j and (151 in pair or pair in [{95,156},{95,142},{142,145},{156,157}]):cg[i][j]=pt(0)
                if pair=={95,157}:cg[i][j]=neg(root(pt(3)))
                if pair in [{145,156},{145,157}]:cg[i][j]=pt(-1)
                if pair=={142,156}:cg[i][j]=neg(root(pt(qs[142]/2)))
        schur=[row[:] for row in cg];pd=[]
        for i in range(4):
            pivot=schur[i][i];assert pivot[0]>0;pd.append(saved(pivot))
            for j in range(i+1,4):
                for kindex in range(j,4):
                    schur[j][kindex]=schur[kindex][j]=minus(schur[j][kindex],div(times(schur[j][i],schur[kindex][i]),pivot))
        cinv,_=inverse(cg)
        v=ivadd(Hres,ivscale(ilinear([full[j] for j in ci],matvec(cinv,[ipair(H,full[j]) for j in ci])),pt(-1)))
        vn,vh=ipair(v,v),ipair(v,H);assert vn[1]<0 and vh[1]<0
        support=sorted(S16|set(ci));assert len(support)==20
        vp=[ipair(v,n) for n in full]
        for j,x in enumerate(vp):
            if j in support:assert R.overlap(x,pt(0))
            else:assert x[1]<0,(sidepair,j,saved(x))
        corners.append(dict(side_pair=sidepair,incident=support,point_bounds=vec(v),norm_bounds=saved(vn),
                            center_pair=saved(vh),all158_pairs=vec(vp),positive_pivots=pd))
    neighbors=sorted(rights|{95,145});assert len(neighbors)==20
    rejections=[]
    for omit in neighbors:
        if omit in [95,145]:continue
        assert {95,145}<=set(neighbors)-{omit}
        rejections.append(dict(omitted=omit,dotted_pair=[95,145]))
    assert len(rejections)==18
    # Four feasible finite vertices, each strict on the other two side walls,
    # give a convex quadrilateral in the projective face plane. Each side is
    # a supporting line of that convex hull, with a nonzero finite segment.
    # Their four inward halfplanes intersect in precisely this quadrilateral;
    # every remaining wall contains its vertices and hence its convex hull.
    sides=[95,142,145,157]
    side_pattern=[]
    for rec in corners:
        row=[]
        for j in sides:
            x=tuple(map(F,rec['all158_pairs'][j]))
            if j in rec['side_pair']:assert x[0]<=0<=x[1];row.append(0)
            else:assert x[1]<0;row.append(-1)
        side_pattern.append(row)
    assert side_pattern==[[0,0,-1,-1],[-1,0,0,-1],[-1,-1,0,0],[0,-1,-1,0]]
    E7=L['E7'];ep=[ipair(E7,n) for n in full]
    assert all(x[1]<0 for x in ep) and ipair(E7,H)[1]<0 and L['q7']>0
    # The old whole compact triangle survives: its three old finite vertices
    # pair strictly with the new wall, hence so does its convex hull.
    triangle=[]
    for rec in L['alltrianglevertices']:
        v=[tuple(map(F,x)) for x in rec['point']];val=ipair(v,r);assert val[1]<0;triangle.append(saved(val))
    files=[f'scripts/{TAG}.py',PREVIOUS,'results/fresh_h20_release144_cap_20261016_review.json',
           'results/fresh_h20_release144_cap_20261016.json',
           'results/fresh_h20_pentagon_closure_scout_20261016.json']
    return dict(status='PASS_EXACT_158_WALL_COMPACT_QUADRILATERAL_CAP',
       sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
       inherited_integrity=ns['integrity'](files),
       definition=dict(contact_frame=ids,target_products=['-sqrt3','-1','0','0'],
          kernel_seed='H projected off S16 and the four contact roots',gram_inverse_pivots=pivots,
          particular_square=saved(p2),kernel_square=saved(k2),radical_square=saved(rad),
          branches=branches,selected_sign=sign,physical_normal_bounds=vec(r),norm_exact=2),
       all157_contacts=contacts,relations=dict(counts),actual_facets=158,
       all158_center_pairs=vec(hp),facet_projection_squares=vec(facets),five_proper_cusp_pairs=vec(cusps),
       new_finite_vertices=2,known_finite_vertices_at_least=184,finite_neighbors=neighbors,rejected_supports=rejections,
       compact_quadrilateral=dict(support=sorted(S16|{151,156}),sides=sides,vertices=corners,
          strict_side_pattern=side_pattern,proof='Four cyclic feasible finite corners and strict nonincident signs certify the whole quadrilateral. All remaining linear inequalities contain its convex hull.'),
       previous_compact_triangle_new_pairs=triangle,
       strict_escape=dict(name='E7',norm=str(L['q7']),physical_point_bounds=vec(E7),
                          all158_pairs=vec(ep),center_pair=saved(ipair(E7,H))),
       scope='Exact158 actual facets, at least184 finite vertices, five proper cusps, a whole compact quadrilateral and the old compact triangle. Strict E7 still proves infinite volume. The quadrilateral is not asserted to be a Coxeter polygon; ambient wall angles are Coxeter.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
