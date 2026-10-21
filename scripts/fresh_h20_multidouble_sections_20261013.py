#!/usr/bin/env python3
"""Cusp-preserving sections of four/eight-tile ambient reflection unions.

Successive seams169,194,71 generate a nonspherical ideal-triangle subgroup.
Each individual double is certified convex by nonpositive seam pairings.
Only accepted actual ridge certificates, not numerical pair proposals,
are geometric conclusions.
"""
import json
import numpy as np
import sympy as sp
import fresh_h20_reflected_sections_double_20261013 as D

prior=D.prior; O=D.O; P=D.P; t=D.t; old=D.old; ROOT=D.ROOT
OUT=ROOT/'data/enumerations/fresh_h20_multidouble_sections_20261013.json'


def double(R,Q,seam):
    matches=np.flatnonzero(np.all(R==seam,axis=1));assert len(matches)==1
    k=int(matches[0]);q=int(seam@Q@seam);p=R@Q@seam
    assert all(int(x)<=0 for j,x in enumerate(p) if j!=k)
    normals=[];seen=set()
    for reflect in (False,True):
        for j,r in enumerate(R):
            if j==k:continue
            if reflect:
                assert (2*int(p[j]))%q==0
                v=r-(2*int(p[j])//q)*seam
            else:v=r
            key=tuple(map(int,v))
            if key not in seen:seen.add(key);normals.append(key)
    return np.array(normals,dtype=np.int64),{'seam_index_zero_based':k,
           'seam_square':q,'before':len(R),'after':len(normals)}


def build(Q,W,G):
    R=np.asarray(W,dtype=np.int64);stages=[];states={}
    for n,label in enumerate((169,194,71),1):
        R,record=double(R,Q,W[label-1]);record['original_seam_one_based']=label
        record['ambient_chambers']=2**n;stages.append(record)
        if n>=2:states[2**n]=R.copy()
    return stages,states


def center_and_data(Q,W,G,R,H,h,e,F):
    sd=prior.sd_from_e(Q,R,H,e,F)
    _,q,_,_=sd;he=old.pdot(Q,e,[P(int(x)) for x in h],F)
    C=[F.rem(q*int(h[d])-he*e[d]) for d in range(22)]
    if any(F.sign(old.pdot(Q,C,[P(int(x)) for x in w],F))>=0 for w in W):
        C=O.four.alternate_center(Q,W,h,prior.sd_from_e(Q,W,G,e,F),F)
    if C is None:return None,None
    assert F.sign(old.pdot(Q,C,C,F))<0
    assert all(F.sign(old.pdot(Q,C,[P(int(x)) for x in w],F))<0 for w in R)
    return C,sd


def test_row(Q,W,G,R,H,h,F,case,e,tau,co,name,target_index=None):
    C,sd=center_and_data(Q,W,G,R,H,h,e,F)
    if C is None:return None
    target=None
    if target_index is not None:
        for b in range(len(R)):
            if b==target_index or H[target_index,b]!=0:continue
            target=D.ridge(Q,R,h,sd,C,F,target_index,b)
            if target is not None:break
        if target is None:return None
    bad=None
    for a,b in D.candidates(Q,R,H,e,C):
        bad=D.ridge(Q,R,h,sd,C,F,a,b,True)
        if bad is not None:break
    link=prior.cusp_proof(Q,W,G,np.array(case['u'],dtype=object),
                         prior.sd_from_e(Q,W,G,e,F),F,case['dropped'])
    row={'name':name,'normal':[O.coeffs(p)for p in e],
         'translation':O.coeffs(tau),'component_coefficients':[O.coeffs(p)for p in co],
         'strict_center':[O.coeffs(p)for p in C],'normal_square':O.coeffs(sd[1]),
         'cusp_link_in_original_wall_labels':link,
         'orthogonal_new_boundary_index':None if target_index is None else target_index+1,
         'target_ridge':target,'actual_forbidden_ridge':bad,
         'status':'NON_COXETER_ACTUAL_RIDGE'if bad else'OPEN_COXETER_CANDIDATE'}
    print(name,len(R),row['status'],None if bad is None else bad['ridge_walls_one_based'],flush=True)
    return row


def run():
    Q,W,G,h=O.base.load();Q=np.asarray(Q,dtype=np.int64);W=np.asarray(W,dtype=np.int64);G=np.asarray(G,dtype=np.int64)
    stages,states=build(Q,W,G);source=json.loads(prior.OUT.read_text());case=source['cases'][0]
    F=O.Field(P(t*t-2),*map(sp.Rational,case['field_interval']))
    decode=lambda p:sp.Poly.from_list(list(map(sp.Rational,p)),t,domain=sp.QQ)
    e0=list(map(decode,case['base_normal']));u=np.array(case['u'],dtype=np.int64)
    assert all(int(W[j-1]@Q@u)<0 for j in (169,194,71))
    runs=[]
    for size,R in states.items():
        H=R@Q@R.T
        assert set(map(tuple,R[(R@Q@u)==0]))==set(map(tuple,W[(W@Q@u)==0]))
        original=case['candidates'][0]
        row=test_row(Q,W,G,R,H,h,F,case,list(map(decode,original['normal'])),
                     decode(original['translation']),list(map(decode,original['component_coefficients'])),
                     f'{size}_tiles_old_control')
        runs.append({'ambient_chambers':size,'boundary_normals':R.tolist(),'candidates':[row]})
    # Explicit fresh targets: a new translated boundary mirror perpendicular
    # to e(tau). The full rank19 cusp shape remains unchanged.
    R=states[8];H=R@Q@R.T;up=R@Q@u;base=prior.sd_from_e(Q,R,H,e0,F)
    old_translations={tuple(r['translation'])for c in source['cases']for r in c['candidates']}
    previous=json.loads(D.OUT.read_text());old_translations|={tuple(r['translation'])for r in previous['candidates']}
    prior_normals=set(map(tuple,states[4]));seen=set();proposals=[]
    for a in range(len(R)):
        if tuple(R[a]) in prior_normals or not up[a]:continue
        tau=F.rem(-base[2][a]/int(up[a]));key=tuple(O.coeffs(tau))
        if key in seen or key in old_translations:continue
        seen.add(key);e=[F.rem(p+tau*int(z))for p,z in zip(e0,u)]
        co=prior.component_coefficients(W,e,case['component'],F);di=case['component'].index(case['dropped']);sg=F.sign(co[di])
        if not sg or any(F.sign(x)*sg>=0 for j,x in enumerate(co)if j!=di):continue
        proposals.append((abs(float(tau.eval(2**.5))),a,tau,e,co))
    attempts=0
    for _,a,tau,e,co in sorted(proposals,key=lambda x:x[:2]):
        if attempts>=12 or len(runs[-1]['candidates'])>=7:break
        attempts+=1
        row=test_row(Q,W,G,R,H,h,F,case,e,tau,co,f'8_tiles_new_target_{a+1}',a)
        if row:runs[-1]['candidates'].append(row)
    data={'schema':'successive-infinite-subgroup-ambient-doubles-v1',
          'source':str(prior.OUT.relative_to(ROOT)),'source_case':0,
          'seams_one_based':[169,194,71],'seam_gram':G[np.ix_([168,193,70],[168,193,70])].tolist(),
          'stages':stages,'field_polynomial':['1','0','-2'],'field_interval':[str(F.lo),str(F.hi)],
          'runs':runs,'new_target_attempts':attempts,
          'scope':'Literal four/eight-chamber convex unions from successive doubles across original169,194,71. '
                  'These seams generate an infinite ideal-triangle subgroup. Two old-control checks and '
                  'the declared new cusp-preserving perpendicular-mirror targets; no general assembly exclusion.'}
    OUT.write_text(json.dumps(data,indent=2)+'\n');print('DONE',[(x['ambient_chambers'],len(x['boundary_normals']),len(x['candidates']))for x in runs],flush=True)


def continuation():
    data=json.loads(OUT.read_text());source=json.loads(prior.OUT.read_text());case=source['cases'][0]
    Q,W,G,h=O.base.load();Q=np.asarray(Q,dtype=np.int64);W=np.asarray(W,dtype=np.int64);G=np.asarray(G,dtype=np.int64)
    F=O.Field(P(t*t-2),*map(sp.Rational,data['field_interval']))
    decode=lambda p:sp.Poly.from_list(list(map(sp.Rational,p)),t,domain=sp.QQ)
    R=np.array(data['runs'][1]['boundary_normals'],dtype=np.int64)
    seam=R[395].copy();u=np.array(case['u'],dtype=np.int64)
    assert seam@Q@W[21]==0 and seam@Q@u<0 and W[21]@Q@u<0
    additions=[]
    for size,normal in ((16,seam),(32,W[21])):
        before=R.copy();R,stage=double(R,Q,normal);H=R@Q@R.T
        stage['seam_normal']=normal.tolist();stage['ambient_chambers']=size
        assert set(map(tuple,R[(R@Q@u)==0]))==set(map(tuple,W[(W@Q@u)==0]))
        selected=data['runs'][1]['candidates'][:1] if size==16 else data['runs'][1]['candidates']
        rows=[]
        for previous in selected:
            e=list(map(decode,previous['normal']));co=list(map(decode,previous['component_coefficients']));tau=decode(previous['translation'])
            row=test_row(Q,W,G,R,H,h,F,case,e,tau,co,f'{size}_tiles_{previous["name"]}')
            oldproof=previous['actual_forbidden_ridge'];pa,pb=[j-1 for j in oldproof['ridge_walls_one_based']]
            oldR=np.array(data['runs'][1]['boundary_normals'],dtype=np.int64)
            relevant=oldR[pb] if pa==0 else oldR[pa]
            if np.array_equal(relevant,seam) or np.array_equal(relevant,W[21]):
                v=W[0]-(2*int(W[0]@Q@relevant)//int(relevant@Q@relevant))*relevant
                aa=int(np.flatnonzero(np.all(R==W[0],axis=1))[0]);bb=int(np.flatnonzero(np.all(R==v,axis=1))[0])
                C=list(map(decode,row['strict_center']));sd=prior.sd_from_e(Q,R,H,e,F)
                obtuse=D.ridge(Q,R,h,sd,C,F,aa,bb,True)
                row['local_doubling_repair']={'old_ridge':oldproof['ridge_walls_one_based'],
                    'seam_normal':relevant.tolist(),'new_pair':[aa+1,bb+1],
                    'actual_obtuse_ridge':obtuse}
            rows.append(row)
        additions.append({'ambient_chambers':size,'stage':stage,'boundary_normals':R.tolist(),'candidates':rows})
    data['continuation']=additions
    data['scope']+=' Continuation doubles at eight-tile boundary396 and original22, orthogonal non-cusp seams, producing16/32 chambers; only the same displayed section planes are tested.'
    OUT.write_text(json.dumps(data,indent=2)+'\n')
    print('CONTINUATION DONE',[(x['ambient_chambers'],len(x['boundary_normals']))for x in additions],flush=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--continuation',action='store_true');args=parser.parse_args()
    continuation() if args.continuation else run()
