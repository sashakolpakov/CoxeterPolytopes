#!/usr/bin/env python3
"""Exact interval conditional escape for the previous literal cap recipes.

Each recipe defines the unique norm-two root in its rational proposal box.
Every root either misses one common strict escape or violates a retained-cap
Coxeter pair. All finite labels are covered by explicit2..120 intervals and
a conservative tail. No numerical graph or finite-volume heuristic is used.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,time
import fresh_h9_split_multiple_edge_escape_20261013 as I
from fresh_h9_split_leaf_3034_safe_graph_20261013 import angle_enclosures,AngleGate,PAIR_SCALE
from fresh_h9_makarov_face_extension_20261013 import K,Z,O,inv,mv,pair,inertia

ROOT=Path(__file__).resolve().parents[1];P='fresh_h9_makarov_joint_completion_20261013'
CORE=ROOT/'results/fresh_h9_makarov_joint_20261013_exact.json'
BANKS=[ROOT/f'results/fresh_h9_makarov_joint_20261013{s}.json'for s in ['', '_profiles']]
OUT=ROOT/f'results/{P}_old_bank_gate.json'
ESCAPE=[153240957355,189417557874,130887763098,72360464506,33886289902,10763857490,-240266086393,-296985201824,-205211687438,-113438150973]
def dec(x):return K(*x)
def ki(x):return I.add(I.rat(x.a),I.mul(I.rat(x.b),I.sqrt(I.rat(5))))
def enc(x):return [str(a)for a in x]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run():
    start=time.time();src=json.loads(CORE.read_text());Q=[[dec(x)for x in row]for row in src['form']];W=[[dec(x)for x in row]for row in src['seed_normals']];c=[dec(x)for x in src['strict_center']]
    assert inertia(Q)==[9,1]
    qi=[[ki(x)for x in row]for row in Q];wi=[[ki(x)for x in row]for row in W];ci=[ki(x)for x in c]
    core=list(wi)
    for row in src['caps']:
        d=I.sqrt(ki(dec(row['radicand'])))
        core.append([I.add(ki(dec(x[0])),I.mul(ki(dec(x[1])),d))for x in row['normal']])
    cov=[I.mv(qi,r)for r in core];e=[K(x)for x in ESCAPE];ei=[ki(x)for x in e];qe=pair(e,e,Q);ce=pair(e,c,Q)
    assert qe.sign()>0 and ce.sign()<0
    oldpairs=[I.dot(r,ei)for r in cov];assert all(p[1]<0 for p in oldpairs)
    assert all(I.dot(r,ci)[1]<0 for r in cov)
    proof=angle_enclosures();gate=AngleGate(proof)
    scale=I.S//PAIR_SCALE;assert I.S==scale*PAIR_SCALE
    angles={row['label']:(int(row['lower'])*scale,int(row['upper'])*scale)for row in proof['finite_angles']}
    angles[2]=I.Z;angles[3]=I.neg(I.O)
    qnorm=[pair(w,w,Q)for w in W];factor=[I.sqrt(ki(x/2))for x in qnorm]
    factor=[I.O if x==K(2)else y for x,y in zip(qnorm,factor)]
    cache={};records=[];counts={};branch_tol=F(1,100000)
    for bi,path in enumerate(BANKS):
        bank=json.loads(path.read_text())
        for index,row in enumerate(bank['caps']):
            if 'old_support'in row:
                ids=row['old_support']+list(range(7,11));labs=row['old_labels']+[2]*4
                for j,m in zip(row['transverse_support'],row['transverse_labels']):labs[5+j]=m
            else:
                ids=row['support']+[row['anchor']];ls=dict(row['contacts']);ls[row['anchor']]=row['anchor_label'];labs=[ls.get(i,2)for i in ids]
            pairs=sorted(zip(ids,labs));key=tuple(i for i,m in pairs);assert len(set(key))==9
            if key not in cache:
                M=[mv(Q,W[i])for i in key];B=None
                for free in range(10):
                    try:B=inv(M+[[K(int(j==free))for j in range(10)]]);break
                    except (AssertionError,StopIteration):pass
                assert B is not None
                z=[r[-1]for r in B];T=[[ki(x)for x in r[:-1]]for r in B]
                cache[key]=(T,[ki(x)for x in z],ki(pair(z,z,Q)),[ki(x)for x in mv(Q,z)])
            T,z,A,Qz=cache[key]
            rhs=[I.mul(angles[m],factor[i])for i,m in pairs];u=I.mv(T,rhs)
            B=I.mul(I.rat(2),I.dot(u,Qz));C=I.sub(I.pair(u,qi,u),I.rat(2))
            assert A[0]>0 or A[1]<0
            D=I.sub(I.mul(B,B),I.mul(I.rat(4),I.mul(A,C)));assert D[0]>0,(bi,index,'disc')
            target=[F(str(x))for x in row['normal']];matches=[]
            for branch in [-1,1]:
                t=I.div(I.add(I.neg(B),I.mul(I.rat(branch),I.sqrt(D))),I.mul(I.rat(2),A))
                r=[I.add(x,I.mul(t,y))for x,y in zip(u,z)]
                inside=all(F(lo,I.S)>a-branch_tol and F(hi,I.S)<a+branch_tol for(lo,hi),a in zip(r,target))
                if inside:matches.append((branch,r))
                else:assert any(F(hi,I.S)<a-branch_tol or F(lo,I.S)>a+branch_tol for(lo,hi),a in zip(r,target)),(bi,index,'branch box')
            assert len(matches)==1,(bi,index,'nonunique branch')
            branch,r=matches[0];ep=I.dot(r,I.mv(qi,ei));norm=I.pair(r,qi,r);assert norm[0]<=2*I.S<=norm[1]
            tag='misses_escape';blocker=None;value=ep
            if ep[1]>=0:
                for j,s in enumerate(cov[11:],11):
                    p=I.dot(r,s);reason=gate.rejection(p[0]//scale,I.ceildiv(p[1],scale))
                    if reason:
                        tag=reason;blocker=j;value=p;break
                assert blocker is not None,(bi,index,'no exact gate',ep)
            counts[tag]=counts.get(tag,0)+1
            records.append({'bank':bi,'entry':index,'branch':branch,'gate':tag,'retained_blocker':blocker,'gate_interval':enc(value)})
            if len(records)%2000==0:print('checked',len(records),'seconds',round(time.time()-start,1),flush=True)
    assert len(records)==25101
    dependencies=[CORE,*BANKS,ROOT/'scripts/fresh_h9_split_multiple_edge_escape_20261013.py',ROOT/'scripts/fresh_h9_split_leaf_3034_safe_graph_20261013.py',ROOT/'scripts/fresh_h9_makarov_face_extension_20261013.py']
    return {'status':'EXACT_CONDITIONAL_OLD_BANK_ESCAPE_WITH_RETAINED_27_FACETS','source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in dependencies},
            'interval_scale':str(I.S),'proposal_box_radius':str(branch_tol),'angle_intervals':proof,'count':len(records),
            'distinct_linear_systems':len(cache),'gate_counts':counts,'escape':{'point':ESCAPE,'norm':qe.enc(),'center_pairing':ce.enc(),'retained_pair_intervals':list(map(enc,oldpairs))},
            'entries':records,'scope':'For all25101 literal archived root recipes, either the exact escape pairing is strict negative or some distinct retained cap has a forbidden Coxeter pairing. Hence no finite acute-compatible subset of this bank can give finite volume while all27 facets are retained. Other roots/labels/recipes and replacements are not excluded.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();result=run()
    if a.verify:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',result['gate_counts'],flush=True)
