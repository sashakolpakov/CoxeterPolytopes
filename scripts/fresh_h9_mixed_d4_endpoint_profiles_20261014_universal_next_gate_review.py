#!/usr/bin/env python3
"""Independent frame/coverage arithmetic; no imports from project/producer.

The inherited interval box engine has a separate proof/cosine audit. Its full
saved replay is run separately; this is not a second box implementation.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_next_gate_review'
PRIMARY='fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_next_gate'
BASE='fresh_h9_mixed_d4_endpoint_pair_20261014'
OLD='fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate_review'
Z=(F(0),F(0)); O=(F(1),F(0)); R=(F(16,11),F(4,11))
def qa(x,y): return (x[0]+y[0],x[1]+y[1])
def qn(x): return (-x[0],-x[1])
def qs(x,y): return qa(x,qn(y))
def qm(x,y): return (x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def qi(x):
    d=x[0]*x[0]-5*x[1]*x[1]; assert d
    return (x[0]/d,-x[1]/d)
def qsign(x):
    a,b=x; sa=(a>0)-(a<0); sb=(b>0)-(b<0)
    if not sa:return sb
    if not sb or sa==sb:return sa
    d=a*a-5*b*b;return sa*((d>0)-(d<0))
def q(v):return tuple(F(x) for x in v)
def e(v):return (q(v[0]),q(v[1]))
def scalar(x):return ((F(x),F(0)),Z)
ZERO=(Z,Z)
def ea(x,y):return (qa(x[0],y[0]),qa(x[1],y[1]))
def en(x):return (qn(x[0]),qn(x[1]))
def es(x,y):return ea(x,en(y))
def em(x,y):return (qa(qm(x[0],y[0]),qm(R,qm(x[1],y[1]))),qa(qm(x[0],y[1]),qm(x[1],y[0])))
def sign(x):
    a,b=x; sa=qsign(a);sb=qsign(b)
    if not sa:return sb
    if not sb or sa==sb:return sa
    return sa*qsign(qs(qm(a,a),qm(R,qm(b,b))))
def dot(x,y):
    out=ZERO
    for a,b in zip(x,y):
        if a!=ZERO and b!=ZERO:out=ea(out,em(a,b))
    return out
def matvec(M,x):return [dot(r,x) for r in M]
def inverse_pd(B):
    n=len(B);a=[list(r)+[O if i==j else Z for j in range(n)] for i,r in enumerate(B)]
    # Positive elimination pivots independently certify positive definiteness.
    for k in range(n):
        assert qsign(a[k][k])>0
        pivot=qi(a[k][k]);a[k]=[qm(v,pivot) for v in a[k]]
        for i in range(n):
            if i==k:continue
            c=a[i][k]
            if c!=Z:a[i]=[qs(v,qm(c,w)) for v,w in zip(a[i],a[k])]
    assert all(a[i][j]==(O if i==j else Z) for i in range(n) for j in range(n))
    return [r[n:] for r in a]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):return json.loads((ROOT/p).read_text())

def run():
    paths=[f'results/{PRIMARY}.json',f'results/{BASE}.json',f'results/{BASE}_next_profiles.json',f'results/{OLD}.json']
    gate,data,profiles,old=[load(p) for p in paths]
    seen=set()
    def manifest(p):
        if p in seen:return
        seen.add(p)
        for s,h in load(p).get('sources',{}).items():
            assert sha(s)==h,s
            if s.endswith('.json'):manifest(s)
    for p in paths:manifest(p)
    assert old['status']=='PASS_UNIVERSAL_ENDPOINT_GATE_PROOF_COVERAGE_COSINE_AUDIT'
    assert gate['status']=='PASS_ALL_71_ENDPOINT_FAMILIES_MISS_SAME_RATIONAL_GAP'
    oldgate=load('results/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.json')
    assert gate['escape']['coordinates']==oldgate['escape']['coordinates']
    for k in ['parameterization','allowed_pairing_relaxation','cosine_enclosures','cosine_error_proof']:
        assert gate[k]==oldgate[k]
    M=[[e(v) for v in row] for row in data['form']]
    N=[[e(v) for v in row] for row in data['normals']]
    H=[e(v) for v in data['center']]
    E=[scalar(F(v)) for v in gate['escape']['coordinates']]
    assert len(M)==10 and len(N)==14 and all(len(v)==10 for v in N)
    assert all(M[i][j]==M[j][i] for i in range(10) for j in range(10))
    MN=[matvec(M,n) for n in N];MH=matvec(M,H);ME=matvec(M,E)
    norm=[dot(n,mn) for n,mn in zip(N,MN)]
    assert all(sign(v)>0 for v in norm)
    assert sign(dot(H,MH))<0 and sign(dot(E,ME))>0 and sign(dot(E,MH))<0
    assert dot(E,ME)==e(gate['escape']['norm'])
    assert dot(E,MH)==e(gate['escape']['future_center_pairing'])
    assert [dot(E,m) for m in MN]==[e(v) for v in gate['escape']['old_pairings']]
    assert all(sign(dot(E,m))<0 and sign(dot(H,m))<0 for m in MN)
    # Rebuild normalized finite pair entries exactly; no square-root approximation.
    candidates=[Z,(F(-1),F(0)),(F(-1,2),F(-1,2))]
    B14={}
    for i in range(14):
        B14[i,i]=(F(2),F(0))
        for j in range(i):
            raw=dot(N[i],MN[j]);assert sign(raw)<=0
            found=[]
            for v in candidates:
                if sign(raw)==qsign(v) and em(scalar(4),em(raw,raw))==em((qm(v,v),Z),em(norm[i],norm[j])):found.append(v)
            if found:B14[i,j]=B14[j,i]=found[0]
    rays=data['complete_ray_census']['rays'];dirs=data['finite_corner_directions']
    open_ids={i for i,r in enumerate(dirs) if r['kind']=='open_to_light_cone'}
    assert len(dirs)==135 and len(open_ids)==71
    assert {int(k) for k in profiles['profiles']}==open_ids=={x['edge_id'] for x in gate['edges']}
    finite_ids=[i for i,r in enumerate(rays) if r['kind']=='finite']
    corners={}
    for k in finite_ids:
        incident=rays[k]['incident'];assert len(incident)==9
        B=[[B14[i,j] for j in incident] for i in incident]
        J=inverse_pd(B);assert all(qsign(v)>=0 for row in J for v in row)
        corners[k]=incident
    assert len(corners)==15
    new_incidents={tuple(v['incident']) for v in data['new_finite_vertices']}
    assert len(new_incidents)==3 and new_incidents <= {tuple(v) for v in corners.values()}
    counts=Counter();newcounts=Counter();frame_records=[]
    for entry in gate['edges']:
        i=entry['edge_id'];r=profiles['profiles'][str(i)];d=dirs[i]
        for key in ['support','start_ray','other_cone_ray','released_wall']:assert r[key]==d[key]
        assert r['direction_index']==r['edge_id']==i
        assert entry['from_new_finite_vertex']==r['from_new_finite_vertex']
        assert r['from_new_finite_vertex']==(tuple(corners[r['start_ray']]) in new_incidents)
        supp=r['support'];outside=set(range(14))-set(supp)
        assert len(supp)==8 and sorted(supp+[r['released_wall']])==corners[r['start_ray']]
        B=[[q(v) for v in row] for row in r['normalized_Gram']]
        J=inverse_pd(B)
        assert J==[[q(v) for v in row] for row in r['normalized_inverse']]
        assert B==[[B14[ii,jj] for jj in supp] for ii in supp]
        X=[e(v) for v in rays[r['start_ray']]['vector']]
        Q=[e(v) for v in rays[r['other_cone_ray']]['vector']]
        MX=matvec(M,X);MQ=matvec(M,Q);xx=en(dot(X,MX));xq=dot(X,MQ)
        Y=[ea(em(xx,v),em(xq,u)) for u,v in zip(X,Q)]
        MY=matvec(M,Y);yy=dot(Y,MY)
        assert sign(xx)>0 and sign(yy)>0 and dot(X,MY)==ZERO and sign(dot(Q,MQ))>0
        assert sign(dot(X,MH))<0 and all(sign(dot(X,m))<=0 for m in MN)
        assert all(dot(X,MN[j])==dot(Y,MN[j])==ZERO for j in supp)
        assert dot(X,MN[r['released_wall']])==ZERO and sign(dot(Y,MN[r['released_wall']]))<0
        # Exact A/sqrt(xx)+C/sqrt(yy)<0 test, avoiding any radical interval.
        for mm in [MN[j] for j in sorted(outside)]+[ME]:
            aa=dot(X,mm);cc=dot(Y,mm);assert sign(aa)<=0
            if sign(cc)>0:
                assert sign(aa)<0
                assert sign(es(em(em(aa,aa),yy),em(em(cc,cc),xx)))>0
            else:assert sign(aa)<0 or sign(cc)<0
        assert sign(dot(X,ME))<0 and sign(dot(Y,ME))>0
        expected=[]
        for p in r['finite_profiles']:
            active=[supp.index(w) for w,m in p['contacts']]
            assert all(J[a][b]==Z for a in active for b in active if a!=b)
            cost=Z
            for w,m in p['contacts']:
                c2={3:O,4:(F(2),F(0)),5:(F(3,2),F(1,2))}[m]
                cost=qa(cost,qm(J[supp.index(w)][supp.index(w)],c2))
            rho=qs((F(2),F(0)),cost);assert qsign(rho)>0 and rho==q(p['residual_norm'])
            expected.append(('finite',tuple(map(tuple,p['contacts']))))
        for p in r['dihedral_tail_profiles']:
            j=supp.index(p['wall']);assert J[j][j]==(F(1,2),F(0))
            assert all(J[j][k]==Z for k in range(8) if k!=j)
            expected.append(('tail',p['wall']))
        actual=[(v['kind'],tuple(map(tuple,v['contacts'])) if v['kind']=='finite' else v['wall']) for v in entry['profiles']]
        assert actual==expected
        for p in entry['profiles']:
            co=p['counts'];leaves=sum(v for k,v in co.items() if k!='boxes')
            assert co['boxes']==2*leaves-1 and p['h_bound']>=1 and p['maximum_depth']<180
            assert all(k in ('boxes','escape') or (k.startswith('wall:') and int(k[5:]) in outside) for k in co)
            counts.update(co);counts[p['kind']+'_profiles']+=1
            if r['from_new_finite_vertex']:newcounts[p['kind']+'_profiles']+=1
        frame_records.append(dict(edge_id=i,start_ray=r['start_ray'],other_cone_ray=r['other_cone_ray'],support=supp))
    assert dict(counts)==gate['counts'] and counts['finite_profiles']==2221 and counts['tail_profiles']==92
    assert dict(newcounts)==gate['new_vertex_counts']==dict(finite_profiles=459,tail_profiles=21)
    assert counts['boxes']==11093 and sum(v for k,v in counts.items() if k=='escape' or k.startswith('wall:'))==6703
    paths+=[f'scripts/{PRIMARY}.py',f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_71_FRAME_PROFILE_AND_QUEUE_REVIEW',sources={p:sha(p) for p in paths},
        recursive_manifests_checked=len(seen),finite_corners_checked=len(corners),edge_frames=frame_records,
        finite_profiles=2221,continuous_tails=92,new_vertex_profile_counts=dict(newcounts),closed_boxes=11093,
        terminal_boxes=6703,arithmetic='Fresh stdlib Fraction Q(sqrt5,sqrt((16+4sqrt5)/11)); no project imports',
        scope='Independent exact frames, full-corner inverse positivity, normalized support Grams, profile costs, escape signs, queue/profile coverage and source/proof review. Full interval saved replay is separate; box decisions and inherited complete profile classification are not independently reimplemented.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();r=run()
    path=ROOT/f'results/{TAG}.json'
    if a.verify:assert r==json.loads(path.read_text())
    else:path.write_text(json.dumps(r,indent=2)+'\n')
    print(r['status'],len(r['edge_frames']),r['finite_profiles'],r['continuous_tails'])
