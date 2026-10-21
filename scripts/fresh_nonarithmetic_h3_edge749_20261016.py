#!/usr/bin/env python3
"""Recentered seven-profile H3 endpoint scout near old edge[7,49]."""
import argparse,json
from hashlib import sha256
import mpmath as mp
import numpy as np
import fresh_nonarithmetic_h3_a2_profiles_20261016 as A
C=A.C;ROOT=A.ROOT;TAG='fresh_nonarithmetic_h3_edge749_20261016';META={}

def setup():
    data=json.loads((ROOT/f'results/{C.TAG}.json').read_text());Q,N,H,E=C.P.geometry()
    old=C.P.P.B.read(C.H3.CAP);new=C.P.P.B.read(C.CAPS)
    N.append(mp.matrix([C.P.mid(x)for x in old['normal_interval']]))
    N.extend(mp.matrix([C.P.mid(x)for x in r['normal_interval']])for r in new['caps']);assert len(N)==55
    oldbasis=[mp.matrix([mp.mpf(x)for x in v])for v in data['basis']]
    ledger=json.loads((ROOT/'results/fresh_nonarithmetic_h3_edge_completion_20261016.json').read_text())
    ends=[r for r in ledger['initial_open_ends']if r['sides']==[7,49]];assert len(ends)==2
    point=mp.matrix([(mp.mpf(str(ends[0]['point'][i]))+mp.mpf(str(ends[1]['point'][i])))/2 for i in range(4)])
    point[0]+=mp.mpf('0.000001')
    T=sum((oldbasis[i]*point[i]for i in range(4)),mp.zeros(10,1));T/=mp.sqrt(-(T.T*Q*T)[0])
    W=mp.matrix([[N[j][i]for j in C.S]for i in range(10)]);J=(W.T*Q*W)**-1
    assert all((T.T*Q*n)[0]<0 for i,n in enumerate(N)if i not in C.S)
    basis=[T]
    for unused in range(3):
        choices=[]
        for n in N:
            v=n-W*J*W.T*Q*n;v+=((v.T*Q*T)[0])*T
            for z in basis[1:]:v-=((v.T*Q*z)[0])*z
            choices.append(((v.T*Q*v)[0],v))
        q,v=max(choices,key=lambda z:z[0]);assert q>mp.mpf('1e-30');basis.append(v/mp.sqrt(q))
    def coord(v):return np.array([float(-(v.T*Q*T)[0])]+[float((v.T*Q*z)[0])for z in basis[1:]])
    NC=np.array([coord(n)for n in N]);white=W*J*W.T*Q*H
    scale=max(mp.sqrt((white.T*Q*white)[0])+1,max((white.T*Q*n)[0]/(-(T.T*Q*n)[0])for i,n in enumerate(N)if i not in C.S))*2
    center=white+scale*T
    assert (center.T*Q*center)[0]<0 and all((center.T*Q*n)[0]<0 for n in N)
    a2=mp.matrix([[N[j][i]for j in(5,6)]for i in range(10)]);J2=(a2.T*Q*a2)**-1;profiles=[]
    for node,label in[(0,2)]+[(node,m)for node in(5,6)for m in(3,4,5)]:
        c=mp.matrix([-2*mp.cos(mp.pi/label)if k==node else 0 for k in(5,6)])
        p=a2*J2*c;offset=np.array([float((p.T*Q*n)[0])for n in N])
        profiles.append(dict(node=node,label=label,residual_square=float(2-(p.T*Q*p)[0]),native_offsets=offset,
            center_lower_bound=float((p.T*Q*center)[0]/scale),physical_projection=[mp.nstr(x,65)for x in p]))
    rays=[]
    for rec in data['actual_face_rays']:
        v=mp.matrix(C.P.P.null([list(Q*N[i])for i in rec['active']]));v/=-(v.T*Q*T)[0]
        cv=coord(v);assert abs(cv[0]-1)<1e-10
        rays.append(dict(mask=sum(1<<i for i in rec['active']),v=cv,kind=rec['kind']))
    META.update(focus_edge=[7,49],old_chart_timelike_seed=[mp.nstr(x,40)for x in point],basis=[[mp.nstr(x,65)for x in v]for v in basis],
        common_center=[mp.nstr(x,65)for x in center],center_square=mp.nstr((center.T*Q*center)[0],65),center_scale=mp.nstr(scale,65),
        old_center_replaced=True,old_projected_normals=NC,
        face_rays=[dict(active=[j for j in range(r['mask'].bit_length())if r['mask']>>j&1],point=r['v'],kind=r['kind'])for r in rays])
    return data,NC,profiles,rays

def run():
    # Reuse only the frozen bounded endpoint calculation. The geometry setup
    # and center are replaced within this process; no frozen file is edited.
    A.setup=setup;d=A.run(10)
    paths=[f'scripts/{TAG}.py',f'scripts/{A.TAG}.py',f'results/{A.TAG}.json',
        'results/fresh_nonarithmetic_h3_edge_completion_20261016.json']
    d['sources'].update({p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths})
    d['status']='NUMERICAL_RECENTERED_H3_SEVEN_PROFILE_ENDPOINT_SCOUT';d['recentered_geometry']=C.clean(META)
    d['scope']='Bounded seven-profile numerical endpoint scout in a new high-precision-constructed Lorentz chart near old[7,49]. Retains all55 old walls but replaces the common center by an explicitly checked strict future timelike vector. Continuous samples and isolated third-contact labels remain bounded2..10; no exact cap or nonexistence claim. The imported exact seven-profile lemma is unchanged.'
    return d

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved recentered profile scout')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
