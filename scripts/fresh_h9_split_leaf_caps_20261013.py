#!/usr/bin/env python3
"""Numerical scout: break S8 leaf symmetry into S4 x S4 before adding caps.

This is discovery code, not an exhaustive algebraic or volume certificate.
The exact selected construction is checked in a separate script.
"""
import argparse, itertools, json
from pathlib import Path
import numpy as np
from scipy.spatial import HalfspaceIntersection, QhullError

ROOT=Path(__file__).resolve().parents[1]
PHI=(1+5**.5)/2
EPS=np.array([1.]*4+[-1.]*4)
LABELS=np.array([0.]+[2*np.cos(np.pi/m) for m in range(3,121)])

def types():
    out=[]
    for l,h in itertools.combinations_with_replacement(range(3),2):
        for k in ([4] if l==h else range(5)):
            out.append((l,h,k))
    return out

def values(typ):
    l,h,k=typ;z=[0.,1.,PHI]
    x=np.array([z[l]]*k+[z[h]]*(4-k)); y=np.array([z[h]]*k+[z[l]]*(4-k))
    return x,y

def equation(a,typ):
    x1,y1=values(typ[0]); x2,y2=values(typ[1])
    x=np.r_[x1,x2]; y=np.r_[y1,y2]; delta=x-y
    t=((x1+y1)[0]+(x2+y2)[0])/2
    h=((x1+y1)[0]-(x2+y2)[0])/2
    d=float(sum(delta)); e=float(EPS@delta)
    L=float(delta@delta)/8-t*t/(2*a)-a*d*d/(32*(1+2*a))-2
    return [64*L-e*e,32*L-32*h*h,-16*h*h]

def coeff(a,b,x,y):
    delta=x-y; t=((x+y)[0]+(x+y)[4])/2;h=((x+y)[0]-(x+y)[4])/2
    w=-delta/8+a*sum(delta)/(32*(1+2*a))+b*(EPS@delta)/(32*(1+2*b))*EPS
    c=t/(16*a)+h/(16*b)*EPS
    return np.r_[c+w,c-w]

def seed(a,b):
    C=-2*np.eye(8)-a*np.ones((8,8))-b*np.outer(EPS,EPS)
    return np.block([[2*np.eye(8),C],[C,2*np.eye(8)]])

def variants(typ):
    l,h,k=typ; z=[0.,1.,PHI]
    for inds in itertools.combinations(range(4),k):
        x=np.full(4,z[h]); y=np.full(4,z[l]);x[list(inds)]=z[l];y[list(inds)]=z[h]
        yield x,y

def bank(a,b,T):
    G=seed(a,b);rows=[];refs=[]
    for typ in itertools.product(T,repeat=2):
        q=equation(a,typ)
        if abs(np.polyval(q,b))>1e-7: continue
        for (x1,y1),(x2,y2) in itertools.product(variants(typ[0]),variants(typ[1])):
            x=np.r_[x1,x2];y=np.r_[y1,y2]
            c=coeff(a,b,x,y)
            if abs(c@G@c-2)>1e-7: raise ValueError('norm')
            rows.append(c);refs.append([list(typ[0]),list(typ[1]),x.tolist(),y.tolist()])
    return G,np.array(rows),refs

def compatible(H):
    dist=np.min(abs((-H)[:,:,None]-LABELS[None,None,:]),axis=2)
    out=(H<=-2+1e-8)|(dist<1e-8)
    np.fill_diagonal(out,False)
    return out

def geometry(G,R,selection):
    vals,Q=np.linalg.eigh(G);idx=list(np.where(vals< -1e-7)[0])+list(np.where(vals>1e-7)[0])
    N=Q[:,idx]*np.sqrt(abs(vals[idx]))
    if sum(N[:,0])<0:N[:,0]*=-1
    walls=np.vstack([N,R[selection]@N]) if len(selection) else N
    assert min(walls[:,0])>0
    halfspaces=np.c_[walls[:,1:],-walls[:,0]]
    try:
        verts=HalfspaceIntersection(halfspaces,np.zeros(9)).intersections
    except QhullError:
        verts=HalfspaceIntersection(halfspaces,np.zeros(9),qhull_options='Qs').intersections
    q=np.einsum('ij,ij->i',verts,verts)
    return {'vertices':len(verts),'outside_vertices':int(sum(q>1+1e-7)),
            'ideal_vertices':int(sum(abs(q-1)<1e-7)),'max_klein_square':float(max(q)),
            'worst_klein_vertex':verts[np.argmax(q)].tolist()}

def run(out):
    T=types();records=json.loads(out.read_text())['records'] if out.exists() else [];counter=0
    for a in [PHI+.5,2.25,2.5,3.,4.]:
        proposals=[]
        for typ in itertools.product(T,repeat=2):
            q=equation(a,typ)
            if max(abs(x) for x in q)<1e-8: continue
            roots=np.roots(np.trim_zeros(np.array(q),'f'))
            for r in roots:
                if abs(r.imag)>1e-8:continue
                b=float(r.real)
                if not (-min(.5,a-2)+1e-7 < b < -1e-7):continue
                if not any(abs(b-t[0])<1e-8 for t in proposals): proposals.append((b,typ,q))
        print({'a':a,'new_parameters':len(proposals)},flush=True)
        for b,typ,eq in proposals:
            if any(abs(r['a']-a)<1e-8 and abs(r['b']-b)<1e-8 for r in records):continue
            G,R,refs=bank(a,b,T);counter+=1
            H=R@G@R.T;compat=compatible(H)
            selections=set()
            for first in range(len(R)):
                selected=[first];available=np.flatnonzero(compat[first]).tolist()
                while available:
                    nxt=max(available,key=lambda j:sum(compat[j,k] for k in available))
                    selected.append(nxt);available=[j for j in available if compat[nxt,j]]
                selections.add(tuple(sorted(selected)))
            # All maximal cliques from greedy roots are evaluated, including small ones.
            best=None
            for selection in selections:
                geom=geometry(G,R,list(selection))
                score=(geom['max_klein_square'],geom['outside_vertices'])
                if best is None or score<best[0]: best=(score,selection,geom)
            rec={'a':a,'b':b,'parameter_type':[list(x) for x in typ],'parameter_equation':eq,
                 'cap_count':len(R),'compatibility_edges':int(np.sum(compat)//2),
                 'greedy_selections':len(selections),'max_selected':max(map(len,selections)),
                 'selected_indices':list(best[1]),'selected_caps':[refs[i] for i in best[1]],
                 'geometry':best[2]}
            records.append(rec)
            print({k:rec[k] for k in ['a','b','cap_count','compatibility_edges','max_selected','geometry']},flush=True)
            out.write_text(json.dumps({'scope':'Numerical S4xS4 scout; finite labels in leaf pairings 2,3,5, five fixed a values, all quadratic b roots from324 orbit equations. Greedy compatible cap sets, not exhaustive clique completion.','records':records},indent=2)+'\n')
    print({'parameter_points':counter,'best':min(r['geometry']['max_klein_square'] for r in records)},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='data/enumerations/fresh_h9_split_leaf_caps_20261013_scout.json');args=p.parse_args()
    run(ROOT/args.output)
