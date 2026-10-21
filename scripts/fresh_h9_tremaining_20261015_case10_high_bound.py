#!/usr/bin/env python3
"""Exact all-label metric bound for the saved T10 family: its high block has m<=11."""
import argparse
import json
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_h20_joint_route_20261015_exact_review as R

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_tremaining_20261015_case10_high_bound'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
STATE='results/fresh_h9_tremaining_20261015_search_case_10_batch_001.json'
DOMAINS='results/fresh_h9_tremaining_20261015_case10_domains.json'
HELPER='scripts/fresh_h20_joint_route_20261015_exact_review.py'
Z=R.ZERO;O=R.ONE
COS={2:Z,3:R.dec('1/2'),4:R.dec([0,'1/2',0,0]),6:R.dec([0,0,'1/2',0])}


def run():
    def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
    seen=set()
    def read(p):
        d=json.loads((ROOT/p).read_text())
        if p not in seen:
            seen.add(p)
            for q,h in d.get('sources',{}).items():
                assert sha(q)==h,q
                if q.endswith('.json'):read(q)
        return d
    geo,state,domains=map(read,(GEOM,STATE,DOMAINS))
    assert domains['cross_block_right']
    g=geo['representatives'][10]['geometry'];V=[set(v['incident'])for v in g['vertices']]
    labels={(a,b):m for a,b,m in state['solution']['pair_labels']}
    blocks=[[3,4],[6,7],[9,10],[12,13]];S=sum(blocks,[])
    assert [b for b in g['blocks']if len(b)==2]==blocks
    assert [labels[tuple(b)]for b in blocks]==[4,4,7,4]
    assert [labels[tuple(sorted((x,y)))]for x in(11,16)for y in(12,13)]==[3]*4
    assert labels[9,16]==2 and labels[10,16]==3
    required=[list(p)for p in combinations(S,2)if list(p)not in blocks]
    required += [[9,11],[10,11],[3,16],[4,16],[6,16],[7,16],[11,16]]
    assert len(required)==31 and len(set(map(tuple,required)))==31
    rows={tuple(r['pair']):r for r in domains['all_domains']}
    failures=[]
    for p in required:
        row=rows[tuple(p)];assert row['domain']==[2]
        witness=next(r for r in row['failures']if r['label']==3)
        support=witness['support'];assert set(support)==V[witness['actual_vertex']]
        M=[]
        for a in support:
            rr=[]
            for b in support:
                if a==b:rr.append(O);continue
                pair=tuple(sorted((a,b)))
                m=3 if pair==tuple(p)else labels.get(pair,2)
                if m==7:m=6
                rr.append(R.neg(COS[m]))
            M.append(rr)
        A=[r[:]for r in M];pivots=[]
        while A:
            pivot=A[0][0];pivots.append(pivot)
            if R.sign(pivot)<=0:break
            q=R.inv(pivot)
            A=[[R.sub(A[i][j],R.mul(R.mul(A[i][0],q),A[0][j]))
                for j in range(1,len(A))]for i in range(1,len(A))]
        assert R.sign(pivots[-1])<=0
        failures.append(dict(pair=p,actual_vertex=witness['actual_vertex'],support=support,
                             relaxed_high_label=6,pivots=[[str(x)for x in t]for t in pivots]))
    # Low D-block inverse cost of the common covector(-1/2,-1/2).
    s=R.dec([0,'1/2',0,0]);u=R.add(O,s)
    invD=[[R.dec(2),R.dec([0,1,0,0])],[R.dec([0,1,0,0]),R.dec(2)]]
    G=[[O,R.neg(s)],[R.neg(s),O]]
    assert [[R.dot(row,col)for col in zip(*invD)]for row in G]==[[O,Z],[Z,O]]
    p=[R.dec('-1/2')]*2
    assert R.dot(p,R.mv(invD,p))==u
    # d=1/(4 sin²(pi/m)); m>=12 gives d>=2+sqrt3.
    d12=R.dec([2,0,1,0])
    assert R.mul(d12,R.dec([2,0,-1,0]))==O
    margin=R.sub(R.mul(s,R.sub(d12,R.dec(2))),O)
    assert margin==R.dec([-1,0,0,'1/2']) and R.sign(margin)>0
    assert R.sub(R.mul(s,R.add(s,d12)),R.mul(u,u))==margin
    paths=[GEOM,STATE,DOMAINS,HELPER,f'scripts/{TAG}.py']
    return dict(status='PASS_T10_SAVED_FAMILY_HIGH_LABEL_AT_MOST_ELEVEN',
                sources={p:sha(p)for p in paths},fixed_mandatory_family=state['solution']['pair_labels'],
                positive_support=S,forced_right_guards=failures,common_own_block=[12,13],
                residual_normals=[11,16],
                projection_data=dict(s='sqrt2/2',u='1+s',d='1/(4sin²(pi/m))',
                                     norm11='-s-F, F>=0',norm16='-s-d',inner='-(1+s)'),
                residual_determinant='(s+F)(s+d)-(1+s)^2 >= s(d-2)-1',
                m_at_least12_margin='sqrt6/2-1>0',remaining_integer_high_labels=[7,8,9,10,11],
                scope='This exact signature obstruction applies to the displayed mandatory-label family '
                      'for every high-block integer m>=12, with all omitted meeting labels and dotted '
                      'distances free. The31 actual-vertex guards are independently checked after relaxing '
                      'm to6; increasing actual off-diagonal magnitudes preserves non-PD by Perron monotonicity. '
                      'The surviving labels7..11 and the whole incidence remain open; no continuation is run.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['forced_right_guards']))
