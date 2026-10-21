#!/usr/bin/env python3
"""Four literal triangles and two marked forests force T >= 285/8."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import fresh_sqrt5_h3_all_actual_count_split_20261016 as E
M=E.M;ROOT=M.ROOT;TAG='fresh_sqrt5_h4_triangle_compensation_barrier_20261016'
TRI='results/fresh_sqrt5_triangle_extension_complete_20261016_rank22.json'

def run():
    E.install()
    data=json.loads((ROOT/TRI).read_text());faces=[];bonuses=[]
    specs=[(152,(0,0,0,2,0,0,0),5,'A7'),
           (276,(1,0,0,2,0,0,0),-3,'A4'),
           (385,(1,0,0,2,0,0,0),-3,'A2'),
           (386,(1,0,0,2,0,0,0),-3,'A2')]
    prior=json.loads((ROOT/f'results/{E.TAG}.json').read_text())
    for index,ext,mark,filler in specs:
        rec=next(r for r in data['records']if r['index']==index);row=rec['row']
        es=[(a,b,m)for(a,b),m in zip(combinations(range(row['nodes']),2),row['labels'])]
        corners=[M.compile_corner(row['white']+[a,b],es,[a,b],[mark])for a,b in combinations(row['black'],2)]
        profile,tables=M.add_bonus(corners,ext)
        for table in tables:profile=M.plus(profile,table[mark])
        js=[j for j in profile if j>=2];assert len(js)==1;bonus=js[0]
        assert profile=={-1:F(2),1:F(-1),bonus:F(1)},(index,profile)
        assert prior['variable_order'][bonus]['type']=='H4'
        assert row['nodes']+8+ext[0]+len(M.MODELS[filler])==22
        bonuses.append(bonus)
        faces.append(dict(index=index,row=row,external_counts=list(ext),facet_root=mark,
            explicit_spectators=['H4','H4']+['A1']*ext[0]+[filler],
            profile=[[j,str(z)]for j,z in sorted(profile.items())],bonus_index=bonus,
            bonus_definition=prior['variable_order'][bonus],corners=rec['corners']))
    assert len(set(bonuses))==4
    fronts=[dict(context=[1,0,0,0,0,0,0,0],marked_hub='A1',marked_component='A1',facet_root=0,residual_components=['D20']),
            dict(context=[1,0,0,3,0,0,0,0],marked_hub='A1',marked_component='A1',facet_root=0,residual_components=['E8'])]
    v,w=[M.literal_profile(f)for f in fronts]
    assert v[-1]==38 and v[0]==v[1]==-19
    assert w[-1]==35 and w[0]==w[1]==-7
    assert all(w[j]==3 for j in bonuses)
    assert all(z>=0 for p in(v,w)for j,z in p.items()if j>=2)
    columns={M.N:F(-1)}
    for p,scale in[(v,F(5,24)),(w,F(19,24))]:
        for j,z in p.items():
            if j>=0:columns[j]=columns.get(j,F(0))+scale*z
    # Each triangle contributes delta - beta <= 0, with multiplier 19/8.
    for j in bonuses:
        columns[1]+=F(19,8);columns[j]-=F(19,8)
    # The model separately fixes distance-two discount to zero.
    columns[0]+=F(19,2)
    assert columns[0]==columns[1]==0 and all(columns[j]==0 for j in bonuses)
    assert all(z>=0 for j,z in columns.items()if j!=M.N)
    lower=F(5,24)*38+F(19,24)*35;assert lower==F(285,8)>35
    paths=[f'scripts/{TAG}.py',f'scripts/{E.TAG}.py',f'results/{E.TAG}.json',TRI,
           'scripts/fresh_sqrt5_marked_h3_center_cross_20261016.py']
    return dict(status='PASS_EXACT_H4_TRIANGLE_COMPENSATION_MODEL_BARRIER',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        dimension=21,facet_dimension=20,target='35',lower_bound=str(lower),faces=faces,
        vertices=[dict(metadata=f,profile=[[j,str(z)]for j,z in sorted(p.items())])for f,p in zip(fronts,[v,w])],
        reduced_vertex_inequalities=['T >= 38 - 19 delta','T >= 35 + 5 delta'],
        exact_combination=dict(vertex_multipliers=['5/24','19/24'],triangle_multipliers=['19/8']*4,
            distance2_zero_multiplier='19/2',remaining_columns=[[j,str(z)]for j,z in sorted(columns.items())if z],rhs=str(-lower)),
        interpretation='The same four H4 coefficients are required by triangles whose H4 corners have A4/A6/isolated-A1 complementary components, and charged three times at A1+3H4+E8. Actual H3 count does not distinguish these contexts. Recording E8 presence in the H4 compensation context would separate this particular witness, without proving feasibility of such a successor.',
        scope='Exact obstruction for the current payment family, with all ordinary H3 cross coefficients already split by actual count and the E6 two-edge exemption. It is not a geometric nonexistence theorem. Literal triangle rows and spherical corners are inherited from the full arithmetic-filtered catalogue; this report does not independently prove an ambient compact realization or complete paid-quadrilateral coverage.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['lower_bound'],flush=True)
