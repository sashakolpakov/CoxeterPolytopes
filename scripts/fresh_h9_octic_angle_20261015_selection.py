#!/usr/bin/env python3
"""Bounded 152-cutter compatibility and common-escape diagnostic, numerical only."""
import argparse,gzip,json,math
from hashlib import sha256
from pathlib import Path
import mpmath as mp
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_octic_angle_20261015_selection'
BASE='results/fresh_h9_fourth_escape_extension_20261015.json'
RESCORE='results/fresh_h9_after_fourth_rescore_20261015.json'
BANK='data/enumerations/fresh_h9_after_third_vertex_free_scout_20261015_bank.json.gz'
ESCAPE=[[['57/4096', '25/4096'], ['17/16', '7/16']], [['23/512', '5/256'], ['13/4', '3/2']], [['17/256', '1/32'], ['81/16', '37/16']], [['91/1024', '5/128'], ['13/2', '3']], [['135/1024', '15/256'], ['39/4', '9/2']], [['55/512', '3/64'], ['127/16', '59/16']], [['21/256', '9/256'], ['6', '11/4']], [['57/1024', '3/128'], ['65/16', '29/16']], [['29/1024', '3/256'], ['17/8', '7/8']], [['-1025/1024', '0'], ['0', '0']]]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def sc(x):
    s=mp.sqrt(5);d=6+16*s/5
    return mp.mpf(x[0][0])+mp.mpf(x[0][1])*s+(mp.mpf(x[1][0])+mp.mpf(x[1][1])*s)*mp.sqrt(d)
def out(x):return float(format(float(x),'.13g'))
def run():
    mp.mp.dps=100;d,r=read(BASE),read(RESCORE)
    for doc in[d,r]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    bank=json.loads(gzip.decompress((ROOT/BANK).read_bytes()))
    Q=mp.matrix([[sc(x)for x in row]for row in d['form']])
    N=[mp.matrix([sc(x)for x in row])*mp.sqrt(2/sc(d['norms'][i]))for i,row in enumerate(d['normals'])]
    M=mp.matrix([[x for x in Q*N[i]]for i in r['independent_pairing_basis']]);J=M**-1
    ids=[x['source_bank_index']for x in r['cutters']];assert len(ids)==152
    X=[J*mp.matrix([mp.mpf(str(bank[i]['normalized_wall_pairings'][k]))for k in r['independent_pairing_basis']])for i in ids]
    P=mp.matrix([[x for x in v]for v in X]);G=P*Q*P.T;gf=np.array(G.tolist(),float)
    ii,jj=np.triu_indices(len(ids),1);vv=gf[ii,jj];assert len(vv)==11476
    finite=-2*np.cos(np.pi/np.arange(2,257));err=np.min(abs(vv[:,None]-finite),axis=1)
    ok=(vv<-2-1e-7)|(abs(vv+2)<1e-7)|(err<1e-7)
    assert not any(ok)
    # The numerical lower bound excludes every m>=7 by monotonicity, too.
    tail=-2*mp.cos(mp.pi/7);assert min(vv)>float(tail)+1e-7
    E=mp.matrix([sc(x)for x in ESCAPE]);H=mp.matrix([sc(x)for x in d['center']])
    ep=[(v.T*Q*E)[0]for v in X];old=[(n.T*Q*E)[0]for n in N]
    eq=(E.T*Q*E)[0];eh=(E.T*Q*H)[0]
    assert eq>0 and eh<0 and max(ep)<-0.2 and max(old)<0
    k=int(np.argmin(vv));near=int(np.argmin(err))
    paths=[BASE,RESCORE,BANK,f'scripts/{TAG}.py']
    return dict(status='PASS_BOUNDED_NUMERICAL_152_SELECTION_NO_COMPATIBLE_PAIR_COMMON_ESCAPE',
        sources={p:sha(p)for p in paths},arithmetic='100-decimal mpmath frame, coordinates and pair matrix from saved binary64 proposal data; binary64 fit comparisons. No exact cap certificate.',
        candidate_count=152,pair_count=11476,compatible_pair_count=0,maximum_compatible_batch_size=1,
        tolerance=1e-7,finite_fit_labels=[2,256],
        numerical_all_label_tail=dict(first_label=7,threshold=out(tail),smallest_pairing=out(min(vv)),
            margin=out(mp.mpf(float(min(vv)))-tail),reason='-2cos(pi/m) decreases for m>=2; all pairs are greater than the m7 value.'),
        mutual_sign_counts=dict(negative=int(sum(vv<-1e-7)),near_zero=int(sum(abs(vv)<=1e-7)),positive=int(sum(vv>1e-7))),
        largest_pairing=out(max(vv)),
        smallest_pair=dict(indices=[ids[ii[k]],ids[jj[k]]],pairing=out(vv[k])),
        closest_finite_fit=dict(indices=[ids[ii[near]],ids[jj[near]]],pairing=out(vv[near]),error=out(err[near]),
            label=int(np.arange(2,257)[np.argmin(abs(vv[near]-finite))])),
        normalized_proposal_norm_max_error=out(max(abs(gf[i,i]-2)for i in range(len(ids)))),
        escape=dict(point=ESCAPE,norm=out(eq),center_pairing=out(eh),old29_max_pairing=out(max(old)),
            proposal_max_pairing=out(max(ep)),proposal_pairings=[dict(source_bank_index=i,pairing=out(g))for i,g in zip(ids,ep)]),
        scope='Numerical only, on exactly the152 saved cutters of the earlier29-wall escape. The common escape remains strict even for their entire union with29 walls while ignoring incompatibility. No all-profile or all-support nonexistence theorem; no new walls added.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();outpath=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(outpath.read_text());print('PASS saved152 numerical compatibility/common escape')
    else:outpath.write_text(json.dumps(d,indent=2)+'\n');print(d['status'])
