#!/usr/bin/env python3
"""High-precision re-solving of selected bounded numerical contact recipes."""
import json,argparse
from pathlib import Path
import mpmath as mp
import fresh_h20_nonright_s61_20261014_scout as S
def refine(d,row):
    c=mp.matrix([mp.mpf(x)for x in row['contacts']])
    # Reconstruct exact named contact cosines rather than retain float values.
    bank=json.loads((S.ROOT/f'results/{S.TAG}.json').read_text());pr=next(x for x in bank['profiles']if x['id']==row['profile'])
    for j,m in zip(pr['support'],pr['labels']):c[j]=2*mp.cos(mp.pi/m)
    mass=2-(c.T*d['J']*c)[0];ii,jj=row['anchors'];L=mp.matrix([[d['P'][i,j]for j in range(3)]for i in[ii,jj]])
    wp=d['W'].T*c;rhs=mp.matrix([wp[i]-2*mp.cos(mp.pi/m)for i,m in zip([ii,jj],row['anchor_labels'])])
    cross=[L[0,1]*L[1,2]-L[0,2]*L[1,1],L[0,2]*L[1,0]-L[0,0]*L[1,2],L[0,0]*L[1,1]-L[0,1]*L[1,0]]
    free=max(range(3),key=lambda i:abs(cross[i]));other=[i for i in range(3)if i!=free]
    M=mp.matrix([[L[i,j]for j in other]for i in range(2)]);p=mp.zeros(3,1);v=mp.zeros(3,1);v[free]=1
    pp=mp.lu_solve(M,rhs);vv=mp.lu_solve(M,-L[:,free])
    for j,x,y in zip(other,pp,vv):p[j]=x;v[j]=y
    dot=lambda a,b:-a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
    qa=dot(v,v);qb=2*dot(p,v);qc=dot(p,p)-mass;disc=qb*qb-4*qa*qc
    assert disc>0
    candidates=[p+(-qb+sg*mp.sqrt(disc))/(2*qa)*v for sg in[-1,1]]
    approx=mp.matrix(row['residual_frame']);branch=min(range(2),key=lambda i:mp.norm(candidates[i]-approx));y=candidates[branch]
    gp=d['P']*y-wp;cut=(d['ep'].T*y+d['AE'].T*c)[0];assert cut>0
    failures=[]
    for i,(g,rel)in enumerate(zip(gp,row['relations'])):
        if rel=='dotted':good=g<-2
        elif rel=='parallel':good=abs(g+2)<mp.mpf('1e-65')
        elif rel=='right':good=abs(g)<mp.mpf('1e-65')
        else:good=abs(g+2*mp.cos(mp.pi/int(rel.split('/')[1])))<mp.mpf('1e-65')
        if not good:failures.append((i,rel,mp.nstr(g,40)))
    return dict(profile=pr,branch=[-1,1][branch],free=free,p=p,v=v,y=y,contacts=c,pairs=gp,cut=cut,failures=failures)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('indices',nargs='*',type=int,default=[45,63]);a=ap.parse_args();d=S.load();bank=json.loads((S.ROOT/f'results/{S.TAG}.json').read_text())
    for i in a.indices:
        r=refine(d,bank['proposals'][i]);print(i,'cut',mp.nstr(r['cut'],40),'branch',r['branch'],'free',r['free'],'y',[mp.nstr(x,50)for x in r['y']],'failures',r['failures'],flush=True)
