#!/usr/bin/env python3
"""Earlier changed-core plane experiments; exact, explicitly scoped failures.

Shares only rational input reconstruction with the standalone 813 checker.
The complete old cusp census is a hash-bound external coverage dependency.
"""
import argparse,json
from fractions import Fraction as F
from fresh_h22_core_replacement_20261013 import P,RAW,CORE,BANK,WALK,read,sha,setup,pair,add,scale,enc

PRE='fresh_h22_core_replacement_20261013_planes'
OUT='results/'+PRE+'.json'
WALK_AUDIT='results/fresh_h21_source_corner_walk_20261013_independent.json'

def rank(M):
    a=[list(map(F,row)) for row in M]; k=0
    for j in range(len(a[0])):
        p=next((i for i in range(k,len(a)) if a[i][j]),None)
        if p is None:continue
        a[k],a[p]=a[p],a[k]; q=a[k][j];a[k]=[x/q for x in a[k]]
        for i in range(k+1,len(a)):
            if a[i][j]:
                q=a[i][j];a[i]=[x-q*y for x,y in zip(a[i],a[k])]
        k+=1
    return k

def run():
    Q,W,C,h,M,R=setup(); norms=[pair(Q,w,w) for w in W]
    audit=read(WALK_AUDIT)
    assert audit['producer_sha256']==sha(WALK) and audit['queue_closed']
    short=[i for i,n in enumerate(norms) if n==2]
    assert len(short)==42 and len(W)==210
    w=W[0]; assert norms[0]==2
    g=[pair(Q,w,x) for x in W]
    removed=[i for i in short if g[i]==-1]
    assert len(removed)==5 and all(pair(Q,W[i],W[j])==0 for i in removed for j in removed if i!=j)
    assert sum(g[i]==0 for i in short)==36
    assert sorted(set(g[i] for i,n in enumerate(norms) if n==4))==[-2,0]
    assert all(pair(Q,h,W[i])==(-1 if norms[i]==2 else -8) for i in range(210))
    assert all(pair(Q,h,c)==-6 and pair(Q,w,c) in (0,-1) for c in C)
    ua=add(h,scale(w,F(1,2)));ub=scale(w,F(-1,2))
    assert pair(Q,ua,ua)+31*pair(Q,ub,ub)==1 and pair(Q,ua,ub)==0
    # u=h-t*w, t=(sqrt31-1)/2: t²+t=15/2 and u²=1.
    # r±=(u,±1/sqrt2) have norm2 and mutual pairing0.
    # The 36 unremoved short neighbors have pairing-1; the own wall has
    # pairing-sqrt31. Long pairings are-8 or-9+sqrt31.
    assert 31>4 and F(28,5)**2>31 and F(17,5)**2>8
    # Worst parent pairing(-13+sqrt31)/2+sqrt2<-2 via simple bounds.
    assert F(3,2)**2>2 and (-13+F(28,5))/2+F(3,2)<-2
    # Center=-29/2+sqrt31/2<0.
    assert F(28,5)<29
    keptshort=[i for i in short if i not in removed]
    shortsum=[sum(W[i][j] for i in keptshort) for j in range(22)]
    assert pair(Q,shortsum,shortsum)==-86
    cusps=[r for r in read(WALK)['corner_orbits'] if r['type']=='ideal']
    assert len(cusps)==18
    bounds=[]
    for row in cusps:
        d=row['primitive_ray']; hd=-pair(Q,h,d)
        assert pair(Q,d,d)==0 and hd>0
        ratios=[F(-pair(Q,W[i],d),hd) for i in short]
        assert all(x>=0 for x in ratios) and max(ratios)<=F(1,3)
        bounds.append({'orbit':row['id'],'max_short_ratio':str(max(ratios))})
    assert max(F(x['max_short_ratio']) for x in bounds)==F(1,3)
    # Any cap u=h-t*w acute/Coxeter to one retained orthogonal short has
    # q<=2. Then u²<=q gives 2t²+2t-16<=0, hence t<3.
    assert 2*3*3+2*3-16>0
    v,vp=read(CORE)['old_projections']; points,lines=read(CORE)['incidence_parts']
    assert pair(Q,v,v)==pair(Q,vp,vp)==-2 and pair(Q,v,vp)==-5 and add(v,vp)==h
    assert all(pair(Q,v,W[i])==-1 and pair(Q,vp,W[i])==0 for i in points)
    assert all(pair(Q,v,W[i])==0 and pair(Q,vp,W[i])==-1 for i in lines)
    assert all(pair(Q,v,W[i])==pair(Q,vp,W[i])==-4 for i,n in enumerate(norms) if n==4)
    assert set(points+lines)==set(short) and rank([W[i] for i in short])==22
    assert all(pair(Q,W[i],W[j])==0 for part in (points,lines) for i in part for j in part if i!=j)
    seen={short[0]}
    while True:
        newer=seen|{j for j in short if any(pair(Q,W[i],W[j])==-1 for i in seen)}
        if newer==seen:break
        seen=newer
    assert seen==set(short)
    # The connected short incidence graph has a unique bipartition up to
    # interchange. Since its walls span the old space, their pairing patterns
    # determine v,vprime uniquely. Every certified diagram symmetry therefore
    # fixes or interchanges these vectors and fixes h=v+vprime.
    ratio_rows=[]
    for row in cusps:
        d=row['primitive_ray']; a=-pair(Q,v,d); b=-pair(Q,vp,d)
        assert a>0 and b>0 and F(1,2)<=F(a,b)<=2
        ratio_rows.append({'orbit':row['id'],'minus_v_pair':str(a),'minus_vprime_pair':str(b)})
    # After deleting points, u=-s*v+vp. Retained lines force q<=2;
    # retained longs force s<=1. Cutting any old cusp requires s>1/2.
    # But u²=-2s²+10s-2 is increasing on[1/2,1] and its minimum5/2>2.
    assert -2*F(1,2)**2+10*F(1,2)-2==F(5,2)>2
    return {'status':'PASS_TWO_SCOPED_CHANGED_CORE_PLANE_OBSTRUCTIONS',
        'sources':{x:sha(x) for x in (RAW,CORE,BANK,WALK,WALK_AUDIT,
            'scripts/fresh_h22_core_replacement_20261013.py')},
        'five_star':{'deleted_short_walls':removed,'facets':809,
            'new_normals':'(h-((sqrt31-1)/2)*W0, +/-1/sqrt2)',
            'old_projection_square':'1','new_norms':['2','2'],
            'new_pair_counts':{'right':1,'pi3':72,'dotted':1542},
            'conditional_nonquasi_cycle':'56-9sqrt31','retained_short_ones_square':'-86',
            'all_old_cusps_strictly_survive':True},
        'complete_cusp_short_ratio_bounds':bounds,
        'one_short_plane':{'normal':'(s*h-t*W0,a)','retains':'all812 except five adjacent short walls',
            'result':'No nonzero acute-compatible Coxeter wall in this plane cuts any inherited old cusp, for any real coefficients and any finite label or dotted relation.',
            's_zero':'KeepingW0 forces t>=0; a retained long wall with W0 dot=-2 forces t<=0. Then both parent heights exclude a nonzero pure-vertical wall.',
            's_positive':'Normalize s=1. A retained short orthogonal toW0 forces q<=2, so t<3; all old cusp ratios are<=1/3.'},
        'incidence_center_cusp_ratios':ratio_rows,
        'one_part_plane':{'deleted':points,'retains':'all21 lines,168 longs,602 signed parent caps',
            'normal':'(-s*v+vprime,a)',
            'result':'No acute-compatible Coxeter wall in span(v,vprime,vertical) cuts any inherited old cusp.',
            'normalization':'A retained line forces the vprime coefficient>=0. If zero, old long acuteness makes the v coefficient>=0, which cannot cut old cusps. Otherwise scale vprime coefficient to1.',
            'proof':'Old long acuteness gives s<=1; an old cusp cut needs s>1/2. A retained line forces q<=2, but old projection square>=5/2 on that interval.'},
        'coverage_dependency':'The full395-orbit source census and241920wall symmetry group are previously independently certified. Every image of a short wall is one ofthe42shortwalls. The full symmetry group either fixes or swaps the two incidence-center projections.',
        'scope':'These two old-projection planes only. The new D20 one-wall replacement lies outside them. No obstruction to arbitrary old-facet deletion or arbitrary changed-core caps.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run()
    if args.verify:assert out==read(OUT)
    else:(P/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])
