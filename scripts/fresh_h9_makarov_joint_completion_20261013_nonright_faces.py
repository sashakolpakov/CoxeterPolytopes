#!/usr/bin/env python3
"""Small numerical nonright-PD7 feasibility pass at the remaining remote gap."""
from pathlib import Path
import json,time
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize
from fresh_h9_makarov_joint_completion_20261013_profiles import setup,ROOT,P
from fresh_h9_joint_frontier_20261013 import spherical_profiles,K,Z,PHI,I
OUT=ROOT/f'results/{P}_nonright_faces.json'
def run():
    start=time.time();src,Q,W,c,E=setup();data=json.loads((ROOT/'results/fresh_h9_joint_frontier_20261013_pd7.json').read_text())
    target=np.array(json.loads((ROOT/f'results/{P}_all_relaxed_volume.json').read_text())['best_point']);target=.999*target+.001*c/(-c@Q@c)
    # Fixed small scope: nearest sixteen exact PD7 faces, one nonright contact.
    supports=[r['support']for r in data['faces']];supports.sort(key=lambda ids:sum(abs(W[ids]@Q@target)))
    G=[[None]*27 for _ in W]
    for i in range(27):G[i][i]=K(2)
    for row in src['all_pair_gates']:
        if row['type']=='dotted':continue
        i,j=row['pair'];G[i][j]=G[j][i]={'m2':Z,'m3':K(-1),'m5':I.neg(PHI)}[row['type']]
    old={tuple(r['support']):r['normal']for r in json.loads((ROOT/f'results/{P}_faces_latest.json').read_text())['caps']}
    rows=[];counts={'supports':16,'profiles':0,'feasible_numerical':0,'target_cut':0}
    for supp in supports[:16]:
        packet=spherical_profiles([[G[i][j]for j in supp]for i in supp],supp)
        profiles=[r['contacts']for r in packet['finite_profiles']if len(r['contacts'])==1]
        profiles += [[[r['wall'],m]]for r in packet['dihedral_tail_profiles']for m in[3,4,5]]
        V=null_space(W[supp]@Q);H=V.T@Q@V;outside=sorted(set(range(27))-set(supp));A=W[outside]@Q@V;ee=target@Q@V;cc=c@Q@V
        Gi=np.linalg.inv(W[supp]@Q@W[supp].T);p=c-(W[supp]@Q@c)@Gi@W[supp];p/=np.sqrt(-(p@Q@p))
        spatial=V[:,0]+(V[:,0]@Q@p)*p
        if spatial@Q@spatial<1e-8:spatial=V[:,1]+(V[:,1]@Q@p)*p
        spatial/=np.sqrt(spatial@Q@spatial)
        for contacts in profiles:
            counts['profiles']+=1;rhs=np.zeros(7)
            for wid,m in contacts:rhs[supp.index(wid)]=-2*np.cos(np.pi/m)
            u=rhs@Gi@W[supp];R=2-rhs@Gi@rhs;assert R>1e-8
            starts=[]
            if tuple(supp)in old:starts.append(V.T@np.array(old[tuple(supp)])*np.sqrt(R/2))
            for boost in[1.,3.]:
                for sign in[-1,1]:starts.append(V.T@(np.sqrt(R)*(np.sinh(boost)*p+sign*np.cosh(boost)*spatial)))
            best=None;d=W[outside]@Q@u;te=target@Q@u;ct=c@Q@u
            for x in starts:
                sol=minimize(lambda z:-ee@z,x,jac=lambda z:-ee,method='SLSQP',constraints=[{'type':'eq','fun':lambda z:z@H@z-R,'jac':lambda z:2*H@z},{'type':'ineq','fun':lambda z:-2-d-A@z,'jac':lambda z:-A},{'type':'ineq','fun':lambda z:-1e-8-ct-cc@z,'jac':lambda z:-cc}],options={'maxiter':120,'ftol':1e-10})
                z=sol.x
                if abs(z@H@z-R)>1e-5 or max(d+A@z)>-2+1e-6 or ct+cc@z>=0:continue
                if best is None or te+ee@z>best[0]:best=(float(te+ee@z),u+V@z)
            if best:
                counts['feasible_numerical']+=1;cut,r=best
                if cut>1e-8:counts['target_cut']+=1;print('CUT',supp,contacts,cut,flush=True)
                rows.append({'support':supp,'contacts':contacts,'normal':r.tolist(),'target_pairing':cut,'all27_pairings':(W@Q@r).tolist()})
        print('support',supp,'counts',counts,flush=True)
    out={'status':'BOUNDED_NUMERICAL_NONRIGHT_PD7_SCOUT','counts':counts,'selected_supports':supports[:16],'target':target.tolist(),'caps':rows,'seconds':time.time()-start,'scope':'Sixteen nearest exact PD7 supports; one nonright contact label3/4/5 admitted by exact spherical-profile Schur tests; four or five numerical optimizer starts per profile. All other old contacts constrained dotted. No exhaustive continuous failure claim.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print('DONE',counts,time.time()-start,flush=True)
if __name__=='__main__':run()
