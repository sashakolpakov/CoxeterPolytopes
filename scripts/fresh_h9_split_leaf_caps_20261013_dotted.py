#!/usr/bin/env python3
"""Constructive numerical search beyond the exact finite-leaf dictionary."""
import json,itertools
from pathlib import Path
import numpy as np
import networkx as nx
from scripts.fresh_h9_split_leaf_caps_20261013 import PHI,EPS,bank,types,coeff,compatible,geometry
ROOT=Path(__file__).resolve().parents[1]

def v(pair,k):
    l,h=pair
    return np.array([l]*k+[h]*(4-k)),np.array([h]*k+[l]*(4-k))
def allv(pair,k):
    l,h=pair
    for inds in itertools.combinations(range(4),k):
        x=np.full(4,h);y=np.full(4,l);x[list(inds)]=l;y[list(inds)]=h
        yield x,y

def run():
    a=PHI+.5;b=-a/(12*a+1);G,R,refs=bank(a,b,types());rows=list(R);profiles=[r[2:4] for r in refs]
    seen={tuple(np.round(np.r_[*r],9)) for r in profiles};new_refs=[]
    vals=[0.,1.,PHI];equations=[]
    for l,h,k in types():
        for low,j in itertools.product(vals,range(5)):
            for swap in [False,True]:
                fixed=(vals[l],vals[h])
                template=lambda d,fixed=fixed,k=k,low=low,j=j,swap=swap: ((fixed,k),((low,d),j))[::(-1 if swap else 1)]
                equations.append(template)
    for l,h,k,j in itertools.product(vals,vals,range(5),range(5)):
        equations.append(lambda d,l=l,h=h,k=k,j=j:(((l,d),k),((h,d),j)))
    for template in equations:
        def q(d):
            (p,k),(r,j)=template(d);x,y=v(p,k);xx,yy=v(r,j);c=coeff(a,b,np.r_[x,xx],np.r_[y,yy]);return c@G@c-2
        f0,fp,fm=q(0.),q(1.),q(-1.);poly=[(fp+fm)/2-f0,(fp-fm)/2,f0]
        for root in np.roots(np.trim_zeros(np.array(poly),'f')):
            if abs(root.imag)>1e-8 or root.real<2-1e-8:continue
            d=float(root.real);(p,k),(r,j)=template(d)
            for (x,y),(xx,yy) in itertools.product(allv(p,k),allv(r,j)):
                x=np.r_[x,xx];y=np.r_[y,yy];key=tuple(np.round(np.r_[x,y],9))
                if key in seen:continue
                seen.add(key);c=coeff(a,b,x,y);assert abs(c@G@c-2)<1e-6
                new_refs.append({'x':x.tolist(),'y':y.tolist(),'dotted':d,'equation':poly})
                profiles.append([x.tolist(),y.tolist()]);rows.append(c)
    R=np.array(rows);print({'all_caps':len(R),'new_dotted_caps':len(new_refs),'equations':len(equations)},flush=True)
    # Build angle graph in blocks to avoid a 3D dense temporary.
    H=R@G@R.T;C=np.zeros(H.shape,dtype=bool)
    for offset in range(0,len(R),128):
        X=H[offset:offset+128];D=np.full(X.shape,np.inf)
        for m in [2]+list(range(3,121)):
            D=np.minimum(D,abs(X+(0. if m==2 else 2*np.cos(np.pi/m))))
        C[offset:offset+128]=(X<=-2+1e-8)|(D<1e-8)
    np.fill_diagonal(C,False)
    graph=nx.from_numpy_array(C)
    print({'edges':graph.number_of_edges()},flush=True)
    W=[]
    for ii,jj in itertools.product(itertools.combinations(range(4),2),itertools.combinations(range(4,8),2)):
        eta=-np.ones(8);eta[list(ii)+list(jj)]=1
        W.append(np.ones(16)*(1/(2*a)+1/1024)+np.r_[eta,-eta])
    W=np.array(W);cuts=R@G@W.T>1e-8
    print({'caps_cutting_balanced_escapes':int(np.any(cuts,axis=1).sum()),
           'new_caps_cutting_balanced_escapes':int(np.any(cuts[len(refs):],axis=1).sum())},flush=True)
    rng=np.random.default_rng(20261013);selections=set()
    for trial in range(120):
        available=np.arange(len(R));chosen=[]
        covered=np.zeros(len(W),dtype=bool)
        while len(available):
            scores=C[np.ix_(available,available)].sum(axis=1)+10000*np.sum(cuts[available][:,~covered],axis=1)
            order=np.argsort(scores)[-min(5,len(available)):]
            nxt=int(available[int(rng.choice(order))]);chosen.append(nxt)
            covered|=cuts[nxt]
            available=available[C[nxt,available]]
        selections.add(tuple(sorted(chosen)))
        if trial%20==0:print({'greedy_trial':trial,'chosen':len(chosen),'unique_selections':len(selections)},flush=True)
    best=None;records=[]
    # Largest compatible sets first; evaluate all distinct greedy sets.
    for selected in sorted(selections,key=lambda c:(-len(c),c)):
        feasible=np.flatnonzero(np.all(R[list(selected)]@G@W.T<0,axis=0))
        if len(feasible):
            wi=int(feasible[0]);q=float(W[wi]@G@W[wi]);
            geom={'numerical_strict_balanced_escape':wi,'lorentz_square':q};score=(q,1)
        else:
            geom=geometry(G,R,list(selected));score=(geom['max_klein_square'],geom['outside_vertices'])
        records.append({'selected':list(selected),'geometry':geom})
        if best is None or score<best:
            best=score;print({'walls':16+len(selected),'geometry':geom},flush=True)
    result={'scope':'Numerical only. Exact seed a=phi+1/2,b=-a/(12a+1). Quadratic dotted leaf caps; finite leaf values0,1,phi, one or both block high values share a dotted parameter.120 seeded randomized greedy clique runs, no exhaustive clique or finite-volume conclusion.',
            'a':a,'b':b,'finite_bank_size':len(refs),'new_dotted_caps':new_refs,
            'all_profiles':profiles,'compatible_edges':graph.number_of_edges(),'greedy_selection_count':len(selections),'records':records}
    (ROOT/'data/enumerations/fresh_h9_split_leaf_caps_20261013_dotted_scout.json').write_text(json.dumps(result,indent=2)+'\n')
    print({'complete':True,'selections':len(selections),'best':best},flush=True)
if __name__=='__main__':run()
