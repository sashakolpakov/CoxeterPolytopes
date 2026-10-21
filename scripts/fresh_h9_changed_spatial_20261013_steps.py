#!/usr/bin/env python3
"""Complete finite profiles and unit intervals on the D9 corner staircase."""
from fresh_h9_changed_spatial_20261013 import *
import subprocess


def eye():return tuple(tuple(O if i==j else Z for j in range(10))for i in range(10))
def mm(A,B):return tuple(tuple(dot(row,col)for col in zip(*B))for row in A)
def boost(j,eps):
    z,x=eps;B=[list(row)for row in eye()]
    B[0][0]=B[j][j]=z;B[0][j]=-x;B[j][0]=-ALPHA*x
    return tuple(map(tuple,B))


def corner(k):
    out=[]
    r=[Z]*10;r[1]=2*O;out.append(tuple(r))
    for j in range(1,k):
        r=[Z]*10;r[j]=-O;r[j+1]=O;out.append(tuple(r))
    for j in range(k+1,9):
        r=[Z]*10;r[j]=-O;r[j+1]=O;out.append(tuple(r))
    if k<9:
        r=[Z]*10;r[9]=-2*O;out.append(tuple(r))
    assert len(out)==9
    return out


def profiles(k):
    free=k+1;roots=corner(k);edge=roots[:k]+roots[k+1:]
    assert all(r[free]==0 for r in edge)and len(edge)==8
    C=tuple(tuple(pair(r,G,s)for s in edge)for r in edge);ci=inverse(C)
    out=[];alphabets=[]
    for q in DYADIC:
        alph=[]
        for r in edge:
            choices={Z}
            for lam in ALLOWED:
                p=square_root(q*pair(r,G,r)*lam/4)
                if p is not None and p.is_integral and (2*p/q).is_integral:choices.add(p)
            alph.append(sorted(choices))
        alphabets.append({'q':q.to_json(),'magnitudes':[enc(a)for a in alph]})
        def visit(pre):
            n=len(pre)
            if n and pair(pre,tuple(row[:n]for row in ci[:n]),pre)>=q:return
            if n==8:
                cs=mv(ci,tuple(-v for v in pre))
                r0=tuple(sum((c*r[j]for c,r in zip(cs,edge)),Z)for j in range(10))
                assert r0[0]==r0[free]==0
                if not all(v.is_integral for v in r0):return
                N=q-pair(r0,G,r0)
                if any(N.sign(e)<=0 for e in range(4)):return
                assert N.is_integral and N.norm().denominator==1
                out.append({'q':q.to_json(),'magnitudes':enc(pre),'fixed':enc(r0),
                            'target':N.to_json(),'absolute_norm':int(N.norm())});return
            for a in alph[n]:visit((*pre,a))
        visit(())
    return {'status':'COMPLETE_SPHERICAL_EDGE_PROFILES','stage':k+1,'free_coordinate':free,
            'edge':menc(edge),'gram':menc(C),'alphabets':alphabets,'profiles':out,
            'absolute_norms':sorted({r['absolute_norm']for r in out})}


def run(k,verify=False):
    assert 1<=k<=8
    tag=f'{PREFIX}_step{k+1}'
    pp=ROOT/'results'/f'{tag}_profiles.json';rp=ROOT/'results'/f'{tag}_raw.txt'
    gp=ROOT/'scripts'/f'{tag}.gp';op=ROOT/'results'/f'{tag}.json'
    pd=profiles(k)
    if verify:assert pd==json.loads(pp.read_text())
    else:
        pp.write_text(json.dumps(pd,indent=2)+'\n')
        script=GP.read_text().replace('ns=vector(13,j,2^(j-1));',f'ns={pd["absolute_norms"]};')
        gp.write_text(script)
        with rp.open('w') as f:
            completed=subprocess.run(['gp','-q',str(gp)],stdout=f,stderr=subprocess.PIPE,text=True,timeout=60)
            assert completed.returncode==0,completed.stderr
    raw=rp.read_text()
    assert 'CERTIFY 1\n'in raw and'BASE_CERTIFY 1\n'in raw and'CLASS 1\n'in raw
    absr={};counts={}
    for line in raw.splitlines():
        if line.startswith('NORM '):
            _,n,c=line.split();absr[int(n)]=[];counts[int(n)]=int(c)
        if line.startswith('SOLUTION '):
            _,n,s=line.split(' ',2);absr[int(n)].append(tuple(map(K.from_json,json.loads(s))))
    assert set(absr)==set(pd['absolute_norms'])and all(len(v)==counts[n]for n,v in absr.items())
    _,_,_,_,_,eps=read_interface();back=(eps[0],-eps[1]);A=eye()
    for j in range(1,k+1):A=mm(boost(j,eps),A)
    AI=inverse(A);assert mm(mm(tuple(zip(*A)),G),A)==G
    assert all(z.is_integral for row in A+AI for z in row)
    free=k+1;upper=ALPHA*eps[1]/eps[0]
    targets={};visited=[];valid=[]
    for i,row in enumerate(pd['profiles']):
        q=K.from_json(row['q']);N=K.from_json(row['target']);r0=tuple(map(K.from_json,row['fixed']))
        if N not in targets:
            orbits=[];excluded=[]
            for z,x in absr[int(N.norm())]:
                assert (z*z-ALPHA*x*x).norm()==N.norm()
                ratio=(z*z-ALPHA*x*x)/N
                if not ratio.is_integral or not ratio.inverse().is_integral:
                    excluded.append({'representative':enc((z,x)),'ratio':ratio.to_json()});continue
                assert all(ratio.sign(e)>0 for e in range(4))
                unit=square_root(ratio);assert unit is not None and unit.is_integral and unit.inverse().is_integral
                z/=unit;x/=unit
                if z<0:z=-z;x=-x
                assert z*z-ALPHA*x*x==N
                orbits.append({'representative':enc((z,x)),'unit_correction':unit.to_json()})
            targets[N]={'target':N.to_json(),'orbits':orbits,'nonunit_rejections':excluded}
        for j,orb in enumerate(targets[N]['orbits']):
            u=tuple(map(K.from_json,orb['representative']));power=0
            while u[1]<=0:u=emul(u,eps);power+=1;assert abs(power)<30
            while emul(u,back)[1]>0:u=emul(u,back);power-=1;assert abs(power)<30
            previous=emul(u,back);first=power;pts=[]
            while ALPHA*u[1]<=upper*u[0]:
                z,x=u;rr=list(r0);rr[0]=x;rr[free]=z;r=mv(AI,rr)
                why=root_reason(r,q)
                row={'profile':i,'orbit':j,'power':power,'root':enc(r),'transformed_root':enc(rr),
                     'q':q.to_json(),'parameter':(ALPHA*x/z).to_json(),'reason':why}
                visited.append(row);pts.append(power)
                if why is None:valid.append(row)
                u=emul(u,eps);power+=1;assert power-first<30
            orb.setdefault('intervals',[]).append({'profile':i,'first_power':first,'last_power':power-1,
                'previous':enc(previous),'next':enc(u),'visited_powers':pts})
    valid.sort(key=lambda r:K.from_json(r['parameter']))
    assert valid,'The proposed expected staircase cap was not a root.'
    nearest=[r for r in valid if r['parameter']==valid[0]['parameter']]
    # The ray of the new facet must imply the other simultaneous crossing
    # mirrors together with the eight retained edge facets. Test exact cone
    # coefficients in each nine-normal basis; no numerical selection.
    edge=[mv(AI,tuple(map(K.from_json,r)))for r in pd['edge']]
    simple=[]
    for row in nearest:
        r=tuple(map(K.from_json,row['root']));frame=edge+[r]
        gram=tuple(tuple(pair(a,G,b)for b in frame)for a in frame);inv=inverse(gram)
        coefficients=[]
        for other in nearest:
            s=tuple(map(K.from_json,other['root']));cs=mv(inv,tuple(pair(s,G,e)for e in frame))
            if any(c<0 for c in cs):break
            assert tuple(sum((c*e[t]for c,e in zip(cs,frame)),Z)for t in range(10))==s
            coefficients.append(enc(cs))
        else:simple.append((row,coefficients))
    assert len(simple)==1
    chosen,cone_coeff=simple[0];root=tuple(map(K.from_json,chosen['root']))
    prior=ROOT/'results'/f'{PREFIX}.json'if k==1 else ROOT/'results'/f'{PREFIX}_step{k}.json'
    old=json.loads(prior.read_text());walls=[tuple(map(K.from_json,r))for r in old['all_walls']]
    assert all(compatible(root,w)for w in walls)
    walls.append(root);pprime=[Z]*10;pprime[0]=O;pprime[free]=K.from_json(chosen['parameter'])
    point=mv(AI,pprime);point=tuple(z/point[0]for z in point)
    assert pair(point,G,point)<0 and all(pair(point,G,w)<=0 for w in walls)
    assert all(pair(point,G,r)==0 for r in edge+[root])
    expected=[Z]*10;expected[0]=eps[1];expected[k]=-O;expected[free]=eps[0]
    expected=mv(AI,expected)
    continues=(root==expected and K.from_json(chosen['parameter'])==upper)
    result={'status':'PASS_COMPLETE_D9_NEXT_EDGE_WITH_CERTIFIED_PARI_INTERFACE',
        'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(pp,rp,gp,prior,OUT)},
        'stage':k+1,'profiles':len(pd['profiles']),'absolute_norm_counts':{str(n):c for n,c in counts.items()},
        'frame_to_axis':menc(A),'frame_from_axis':menc(AI),'relative_unit':enc(eps),
        'upper_parameter':upper.to_json(),'targets':list(targets.values()),'visited_points':visited,
        'primitive_reflective_survivors':valid,'nearest_candidates':nearest,
        'chosen_root':enc(root),'chosen_norm':pair(root,G,root).to_json(),
        'simultaneous_root_cone_coefficients':cone_coeff,'endpoint':enc(point),
        'endpoint_norm':pair(point,G,point).to_json(),'endpoint_incident':menc(edge+[root]),
        'all_walls':menc(walls),'wall_count':len(walls),'expected_staircase_continues':continues,
        'scope':'Complete first-contact census on this finite segment, using all integral spherical profiles and certified full relative norm orbits. Exact positive-cone comparisons select the full local chamber when several mirrors cross together. No compatibility rejection is used before selecting the first contact; all inherited walls are checked afterward. Norm/unit interface remains an explicit PARI dependency.'}
    if verify:assert result==json.loads(op.read_text())
    else:op.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'stage',k+1,'profiles',len(pd['profiles']),'points',len(visited),
          'first mirrors',len(nearest),'continues',continues,flush=True)
    return continues


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',type=int,default=2);p.add_argument('--verify',action='store_true')
    args=p.parse_args();run(args.stage-1,args.verify)
