#!/usr/bin/env python3
"""Independent interval proof for four entire orthogonal H20 edge families.

Only the independent145 physical certificate supplies geometric intervals.
The producer gate is compared after all inequalities have been proved.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_h4_edge_escape_gate_20261016_review'
GEOMETRY='results/fresh_h20_residual_h4_corner_20261016_review.json'
PRIMARY='results/fresh_h20_h4_edge_escape_gate_20261016.json'


class Interval:
    def __init__(self,lo=0,hi=None):self.lo=F(lo);self.hi=F(lo if hi is None else hi);assert self.lo<=self.hi
    def __neg__(self):return Interval(-self.hi,-self.lo)
    def __add__(self,other):
        other=lift(other);return Interval(self.lo+other.lo,self.hi+other.hi)
    __radd__=__add__
    def __sub__(self,other):return self+-lift(other)
    def __rsub__(self,other):return lift(other)+-self
    def __mul__(self,other):
        other=lift(other);vs=[a*b for a in (self.lo,self.hi) for b in (other.lo,other.hi)]
        return Interval(min(vs),max(vs))
    __rmul__=__mul__
    def __truediv__(self,other):
        other=lift(other);assert other.lo>0 or other.hi<0
        return self*Interval(1/other.hi,1/other.lo)
    def saved(self):return list(map(str,(self.lo,self.hi)))


def lift(x):return x if isinstance(x,Interval) else Interval(x)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()


def run():
    g=json.loads((ROOT/GEOMETRY).read_text())
    assert g['status']=='PASS_INDEPENDENT_145_WALL_RESIDUAL_FINITE_CORNER_AND_ESCAPE_RECTANGLE'
    for p,h in g['sources'].items():assert digest(p)==h,p
    q=list(map(F,g['four_norms']));c=list(map(F,g['corner']['chart']))
    metric=[F(-1),F(1),F(1),F(1,2),F(1,2)]
    dot=lambda v,w:sum(a*b*d for a,b,d in zip(v,w,metric))
    cq=dot(c,c);assert str(cq)==g['corner']['norm'] and cq<0
    assert all(v>0 for v in q)
    assert g['four_gram']==[[str(q[i] if i==j else F(0)) for j in range(4)] for i in range(4)]
    assert len(g['corner']['support16_pivots'])==16 and all(F(x)>0 for x in g['corner']['support16_pivots'])
    e=[F(1),F(-1,100),F(-49,50),F(-3,10),F(-1,20)]
    ce=dot(c,e);assert ce<0
    pairs={tuple(row['walls']):Interval(*row['bounds']) for row in g['new_pairs']}
    cp=lambda j:Interval(*g['corner']['all145_pairs'][j])
    upper=[Interval(*g['strict_escape']['all145_pairs'][i])/ce for i in range(141,145)]
    assert all(v.lo>0 for v in upper)
    direct=[]
    for k in range(3):
        a,b=pairs[(133,141+k)],cp(133)
        assert a.hi<0 and b.hi<0
        lower=a/b
        assert lower.lo>upper[k].hi
        direct.append(dict(axis=141+k,blocking_wall=133,acute_lower=lower.saved(),
                           nonnegative_escape_upper=upper[k].saved(),gap=(lower-upper[k]).saved()))
    a,b=pairs[(132,144)],cp(132)
    assert a.hi<0 and b.hi<0
    lower=a/b;U=upper[3]
    assert 0<lower.lo and lower.hi<U.lo
    qU=q[3]+cq*U*U;gU=-a+U*b
    assert qU.lo>0 and gU.hi<0
    margin=2*qU-4*gU*gU;assert margin.lo>0
    # For lower < lambda <= U, -g132 grows from zero and q decreases.
    # Thus 0 < 4*g132^2/(2*q) <= 4*g132(U)^2/(2*q(U)) < 1.
    # At the sole right-contact endpoint lambda=lower, root31 is forbidden.
    endpoint_pair=-pairs[(31,144)]+lower*cp(31)
    endpoint_norm=q[3]+cq*lower*lower
    assert endpoint_pair.hi<0 and endpoint_norm.lo>0
    endpoint_margin=2*endpoint_norm-4*endpoint_pair*endpoint_pair
    assert endpoint_margin.lo>0
    primary=json.loads((ROOT/PRIMARY).read_text())
    assert primary['status']=='PASS_ALL_PARAMETER_FOUR_ORTHOGONAL_EDGE_E4_GATES'
    assert F(primary['corner_E4_pair'])==ce
    for v,raw in zip(upper,primary['E4_cut_lambda_uppers']):assert v.lo<=F(raw)<=v.hi
    paths=[f'scripts/{TAG}.py',GEOMETRY,PRIMARY]
    return dict(status='PASS_INDEPENDENT_ALL_PARAMETER_FOUR_EDGE_ESCAPE_GATES',
                sources={p:digest(p) for p in paths},corner_norm=str(cq),corner_escape_pair=str(ce),
                direct_acute_obstructions=direct,
                fourth_axis=dict(acute_lower=lower.saved(),nonnegative_escape_upper=U.saved(),
                                 norm_at_upper=qU.saved(),pair132_at_upper=gU.saved(),
                                 forbidden_cycle_margin=margin.saved(),endpoint_pair31=endpoint_pair.saved(),
                                 endpoint_norm=endpoint_norm.saved(),endpoint_forbidden_cycle_margin=endpoint_margin.saved()),
                proof='A positive rank19 orthogonal support leaves span(n_i,C). Acuteness and positive norm force r=alpha*(-n_i+lambda*C), alpha>0. Nonnegative E4 pairing bounds lambda above. Three axes violate old133 acuteness. On the fourth, old132 acuteness bounds lambda below; the entire positive-norm interval has forbidden cycle strictlybetween0and1, with old31 covering its unique right endpoint.',
                scope='Every acute Coxeter additional wall right to S16 and three of the four residual caps has strictlynegative E4 pairing. All real parameters, scales and labels; no retained-center condition. Bent contacts and walls leaving S16 remain open. Independent geometry intervals and separate standard-library arithmetic; no producer arithmetic import.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)
