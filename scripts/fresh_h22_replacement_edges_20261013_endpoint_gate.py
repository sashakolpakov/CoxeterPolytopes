#!/usr/bin/env python3
"""All-label spherical/affine endpoint obstruction for the complete22star.

Independently replays the two rational cone identities, then uses exact
Schur bounds. No LP solver, finite-label scan, or floating arithmetic.
"""
import argparse,json
from fractions import Fraction as F
import fresh_h22_replacement_edges_20261013 as P
R=P.R;TAG=P.TAG+'_endpoint_gate';OUT='results/'+TAG+'.json'
DUAL='results/fresh_h22_replacement_escape_dual_20261013.json'

def run(verify=False):
    src,M,W,cov,h,v,inc,U=P.setup();d=R.read(DUAL)
    for p,s in d['sources'].items():assert R.sha(p)==s,p
    escape=R.ratvec(list(map(F,d['escape'])))
    assert escape==R.ratvec(src['strict_escape']['coordinates'])
    assert R.pair(escape,R.mv(M,escape))==R.C(1586)
    assert R.sign(R.pair(escape,R.mv(M,h)))<0
    mandatory=[7,168]
    assert d['mandatory_strictly_negative_pairs_for_any_escape_cutter']==mandatory
    for target in mandatory:
        terms={int(i):F(a)for i,a in d['single_negative_identities'][str(target)].items()}
        assert terms[target]<0 and all(a>0 for i,a in terms.items()if i!=target)
        assert all(i in W for i in terms)
        value=[R.psum(R.scale(W[i][j],a)for i,a in terms.items())for j in range(23)]
        assert value==escape
    A=[]
    for i in inc:
        rr=[]
        for j in inc:
            a=R.pair(W[i],cov[j]);assert set(a)<={1};rr.append(a.get(1,F(0)))
        A.append(rr)
    R.pd(A);AI=P.inverse(A);rows=[]
    for k,released in enumerate(inc):
        ii=[i for i in range(22)if i!=k];edge=[inc[i]for i in ii]
        G=[[A[i][j]for j in ii]for i in ii]
        B=[[AI[i][j]-AI[i][k]*AI[k][j]/AI[k][k]for j in ii]for i in ii]
        assert all(x>=0 for rr in B for x in rr)
        assert all(sum(B[i][k]*G[k][j]for k in range(21))==int(i==j)for i in range(21)for j in range(21))
        present=[i for i in mandatory if i in edge];assert present
        magnitude={i:R.sqrtq(A[inc.index(i)][inc.index(i)]/2)for i in present}
        lower=R.psum(R.scale(R.mul(magnitude[i],magnitude[j]),B[edge.index(i)][edge.index(j)])for i in present for j in present)
        margin=R.add(lower,R.C(-2));assert R.sign(margin)>0
        assert R.sign(R.add(lower,R.C(-3)))>=0
        rows.append({'released':released,'edge':edge,'mandatory_contacts_in_edge':present,
                     'minimum_nonright_pairing_magnitudes':{str(i):R.enc(a)for i,a in magnitude.items()},
                     'projected_norm_lower_bound':R.enc(lower),'excess_over_norm_two':R.enc(margin),
                     'inverse_nonnegative':True,'conclusion':'NO_FINITE_OR_IDEAL_COXETER_ENDPOINT_THAT_CUTS_FIXED_ESCAPE'})
    assert len(rows)==22
    out={'status':'PASS_EXACT_ALL_LABEL22_EDGE_ENDPOINT_CUTTER_OBSTRUCTION',
         'sources':{p:R.sha(p)for p in(P.SOURCE,DUAL,'results/'+P.TAG+'.json','scripts/'+TAG+'.py','scripts/'+P.TAG+'.py','scripts/fresh_h22_replacement_escape_dual_20261013.py','scripts/fresh_h22_multiquadratic_20261013_review.py')},
         'mandatory_nonright_walls':mandatory,'independently_replayed_vector_identities':2,
         'edges_checked':22,'minimum_projected_norm_lower_bound':'3','new_normal_squared_norm':'2','rows':rows,
         'scope':'For any real positive-norm Coxeter normal acute to the813retained walls and cutting the savednorm1586escape, both wall7 and wall168 contacts are nonright. On every outgoing21wall edge support from the specified actualfinite vertex, minimalCoxeter nonright pairings already give projectednorm at least3 after normalizing the newroot to norm2. Therefore its22Gram cannot bePSD: the cap cannot meet any of these22edges at a finite or ideal endpoint. No arithmetic or finite-label bound. This does not rule out intermediate caps that leave this escape feasible, caps elsewhere, subsequent vertices, retained-wall replacements, or all finite-volume completions.'}
    if verify:assert out==R.read(OUT)
    else:(R.P/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'22supports;min3>2',flush=True);return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');run(ap.parse_args().verify)
