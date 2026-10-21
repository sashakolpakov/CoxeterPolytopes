#!/usr/bin/env python3
"""Order16 ideal reductions and a combined finite cross-component budget."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_ideal7_finite_cross_payments_20261017 as P

ROOT=P.ROOT
TAG='fresh_order16_finite_cross_payments_20261017'
EMBED='fresh_order16_preferred_embeddings_20261017'
CYCLE='fresh_cycle16_half_coupon_transfer_20261017'

def coupon(native,ids,a,b):
    old=P.coupon(native,ids,a,b)
    if old:return old
    ca,cb=P.component(native,ids,a),P.component(native,ids,b)
    if ca==cb:return F()
    for small,large,x,y in((ca,cb,a,b),(cb,ca,b,a)):
        if(len(small),len(large))not in((1,15),(2,13)):continue
        if any(native[i][j]not in(0,-1)for i,j in combinations(small|large,2)):continue
        adj={i:{j for j in large if j!=i and native[i][j]}for i in large}
        forks=[i for i in large if len(adj[i])==3]
        if len(forks)!=1 or any(len(adj[i])>3 for i in large):continue
        f=forks[0];leaves=[j for j in adj[f]if len(adj[j])==1]
        if len(leaves)==2 and y in[f,*leaves]:return F(1,2)
    return F()

def exceptional():
    embed=json.loads((ROOT/f'results/{EMBED}.json').read_text())
    old=json.loads((ROOT/f'results/{P.OLD}.json').read_text())
    diagrams={r['index']:r for r in json.loads((ROOT/f'results/{P.CENSUS}.json').read_text())['diagrams']}
    g=embed['targets']['D']['target'];adj,dist,_=P.distance_data(g)
    dropped={p for p in combinations(range(16),2)if dist[p[0]][p[1]]>=12 or len(adj[p[0]])==len(adj[p[1]])==3}
    assert len(dropped)==9
    affected=[r['case']for r in embed['targets']['D']['case_embeddings']if set(map(tuple,r['marked_images']))&dropped]
    assert affected==[384,394,401,648,663,676,814,829,842]
    native,gram,sc,mp=P.R.combined(diagrams,25,25,((8,9),(9,8)))
    assert sc==[1]*18 and P.R.inertia(gram)==[17,1,0]
    inverse=P.R.inverse(gram)
    assert inverse[3][3]==inverse[13][13]==4 and inverse[3][13]==-4
    affine=[i for i in range(18)if i not in(3,13)]
    block=[[gram[i][j]for j in affine]for i in affine]
    assert P.R.inertia(block)==[15,0,1]
    kernel=P.R.nullspace(block);assert len(kernel)==1
    z=kernel[0]
    if all(v<0 for v in z):z=[-v for v in z]
    assert all(v>0 for v in z)
    records=[];needed=[]
    for index in affected:
        r=old['isolation']['cases'][index]
        assert r['left'][0]==r['right'][0]==25 and r['match']==[[8,9],[9,8]]
        assert r['method']=='lorentz_residual'and r['possible_contacts']==[3,13]and not r['single_contact_candidates']
        left=r['left'][1];right=[mp[i]for i in r['right'][1]];black=left+right
        corners=P.case_payments(native,sc,black,left,right);total=F();base=F()
        white=[i for i in range(18)if i not in black]
        for c in corners:
            a,b=c['marks'];extra=coupon(native,white+[a,b],a,b)if c['finite']else F()
            c['cross_coupon']=str(extra);c['weight']=str(F(c['base'])+extra)
            base+=F(c['base']);total+=F(c['weight'])
        assert total>=1
        if base<1:needed.append(index)
        records.append(dict(case=index,black=black,corners=corners,base_total=str(base),total=str(total)))
    assert needed==[384,401,814,842]
    return dict(gram=[[str(v)for v in row]for row in gram],signature=[17,1,0],
        inverse_contact_block=[[str(inverse[i][j])for j in(3,13)]for i in(3,13)],
        affine_nodes=affine,affine_positive_kernel=list(map(str,z)),dropped_pairs=[list(p)for p in sorted(dropped)],
        affected_cases=affected,coupon_needed_cases=needed,cases=records,
        remaining_D_pairs=[p for p in embed['targets']['D']['paid_pairs']if tuple(p)not in dropped],
        sealing='Every extra supporting white is orthogonal to the affine factor. Contact with both3/13 pays a path2. Single contact has projection norm>=4>2 in a Lorentzian core with positive orthogonal complement, impossible. Thus an otherwise unpaid star has every extra white orthogonal to the whole core.')

def finite_budget():
    rows=[]
    for n in range(87,109):
        count=0;maximum=F()
        # t counts D15,D16,A15, v counts D13,D14. Rank/slack lower bounds
        # overestimate both the feasible count domain and the coupon cost.
        for t in range(n//15+1):
            for v in range((n-15*t)//13+1):
                rem=n-15*t-13*v
                for b in range(rem//2+1):
                    for a in{0,rem-2*b}:
                        diff=F(3*a*t,2)+3*b*v-9*a-17*b-24*t-35*v
                        assert diff<=0;maximum=max(maximum,diff);count+=1
        rows.append(dict(dimension=n,endpoint_checks=count,maximum_excess=str(maximum)))
    return dict(rows=rows,formula='cost<=9n+a*(3t/2-9)+b*(3v-17)-24t-35v',
        rank_constraint='a+2b+15t+13v<=n',
        component_meaning='a=A1 count,b=A2 count,t=D15+D16+A15 counts,v=D13+D14 counts',
        analytic_cases=dict(t7_remaining_rank=3,v6_bmax=15,v7_bmax=8,v8_bmax=2),
        scope='Includes all retained91 finite coupons, new D15/A1,D13/A2 fork coupons and A15-center/A1 half-coupons. Other whole finite components have cost<=9rank. Endpoint enumeration corroborates the written proof for every n<=108.')

def run():
    exc=exceptional();finite=finite_budget()
    embed=json.loads((ROOT/f'results/{EMBED}.json').read_text())
    assert{t:r['cost']for t,r in embed['targets'].items()}==dict(A=112,B=107,C=103,D=120)
    assert len(exc['remaining_D_pairs'])==111
    # Every inherited affine linear branch has cost/m<=217/29. Its finite
    # extrema and decreasing infinite branch give the all-rank certificate.
    old=json.loads((ROOT/f'results/{P.TAG}.json').read_text());lam=F(217,29);branches=[]
    for r in old['affine_cost']['all_rank_linear_branches']:
        slope=lam-r['slope'];lo,hi=r['lower_order'],r['upper_order']
        assert hi is not None or slope>0
        at=hi if slope<0 else lo;slack=slope*at-r['intercept'];assert slack>=0
        branches.append(dict(kind=r['kind'],lower_order=lo,upper_order=hi,slope=r['slope'],intercept=r['intercept'],
            minimum_at=at,minimum_slack=str(slack)))
    assert lam>7 and F(127,17)<lam
    stems=(EMBED,CYCLE,P.TAG,P.OLD,P.CENSUS)
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    for p in paths:assert p.is_file(),p
    return dict(status='PASS_ORDER16_FINITE_CROSS_PAYMENTS_AND_PIECEWISE_IDEAL_BUDGET',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        exceptional_D15=exc,finite_cost=finite,order16_costs=dict(A=112,B=107,C=103,D=111),
        large_order_linear_branches=branches,large_normalized_bound=str(lam),
        ideal_rank_envelope={'1_to_14':'r/2','15':'7','16':'127/17','17_and_larger':'217/29'},
        scope='Complete local transfer using inherited91 rule, order16 preferred-component embeddings, separate universal antipodal-cycle replacement and nine sealed D15 exceptions. Every triangle receives2 and quadrilateral1, finite cost<=9n for n<=108, and the displayed all-rank ideal envelope holds. Dimension transfer is separate.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=json.loads(json.dumps(run()))
    p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
