#!/usr/bin/env python3
"""Independent cofactor-determinant check of the balanced E8 obstruction."""
import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nikulin_balanced_height_20261015_review'
PRIMARY='fresh_nikulin_balanced_height_20261015'

def matrix(n,edges):
    G=2*s.eye(n)
    for a,b in edges:G[a,b]=G[b,a]=-1
    return G

def costs(G):
    out=[]
    for i in range(G.rows):
        ids=[j for j in range(G.rows) if j!=i];B=G.extract(ids,ids)
        assert all(B[:k,:k].det()>0 for k in range(1,B.rows+1))
        det=B.det()
        out.append({i:s.S(0),**{j:s.cancel(B.minor_submatrix(k,k).det()/det) for k,j in enumerate(ids)}})
    return out

def run():
    primary=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    assert primary['status']=='PASS_BALANCED_HEIGHT_UNIT_CHARGE_OBSTRUCTION'
    G=matrix(9,[(0,1),(0,2),(2,3),(0,4),(4,5),(5,6),(6,7),(7,8)])
    delta=s.Matrix([6,3,4,2,5,4,3,2,1]);assert G*delta==s.zeros(9,1)
    assert G==s.Matrix(primary['E8']['gram']) and list(map(str,delta))==primary['E8']['kernel']
    C=costs(G);rows=[];legal_count=0
    for i,j in combinations(range(9),2):
        old=next(x for x in primary['E8']['distinct_omission_pairs'] if x['omissions']==[i,j])
        permitted=[];assert len(old['unit_common_rows'])==9
        for k in range(9):
            a,b=C[i][k],C[j][k];total=a+b;legal=max(a,b)<=2
            other=next(x for x in old['unit_common_rows'] if x['node']==k)
            unit=[int(l==k) for l in range(9)]
            assert other['row_p']==other['row_q']==unit
            assert other['heights']==[str(delta[k])]*2
            assert other['M']=='1' and other['weighted_energy']=='0'
            assert list(map(s.Rational,other['endpoint_costs']))==[a,b]
            assert s.Rational(other['total_cost'])==total
            assert other['endpoint_budget_legal']==bool(legal)
            assert other['violates_unit_charge']==bool(legal and total<s.Rational(16,9))
            legal_count+=int(bool(legal))
            if legal:permitted.append((total,k,a,b))
        total,k,a,b=min(permitted);assert 0<total<s.Rational(16,9)
        alpha=s.Rational(16,9)/total
        assert s.Rational(old['necessary_charge_from_common_rows'])==alpha
        node=old['selected_witness']['node']
        assert max(C[i][node],C[j][node])<=2 and C[i][node]+C[j][node]<s.Rational(16,9)
        rows.append(dict(omissions=[i,j],node=k,costs=list(map(str,(a,b))),necessary_coefficient=str(alpha)))
    assert len(rows)==36 and legal_count==179
    best=min(s.Rational(x['necessary_coefficient']) for x in rows)
    minimizers=[x['omissions'] for x in rows if s.Rational(x['necessary_coefficient'])==best]
    assert best==s.Rational(8,7) and minimizers==[[1,6]]
    assert all(C[i][i]==0 for i in range(9))
    E=matrix(8,[(0,1),(0,2),(2,3),(3,4),(0,5),(5,6),(6,7)])
    EC=costs(E);assert E==s.Matrix(primary['current_E7']['gram'])
    assert EC[3][4]==s.Rational(1,2) and EC[6][4]==1
    assert EC[3][4]+EC[6][4]<s.Rational(14,9)
    digest=lambda p:sha256((ROOT/p).read_bytes()).hexdigest()
    files=[f'scripts/{TAG}.py',f'scripts/{PRIMARY}.py',f'results/{PRIMARY}.json',f'results/{PRIMARY}.md']
    for p,h in primary['sources'].items():assert digest(p)==h
    return dict(status='PASS_INDEPENDENT_COFACTOR_BALANCED_HEIGHT_OBSTRUCTION',
                sources={p:digest(p) for p in files},distinct_pairs=rows,common_rows_checked=324,
                endpoint_budget_legal_rows=legal_count,identical_omissions_trivially_excluded=9,
                best_necessary_charge=str(best),best_pairs=minimizers,
                current_E7_endpoint_costs=['1/2','1'],upper_lemma_proved=False,
                scope='Separate SymPy determinant-cofactor arithmetic, no project imports. Refutes the preset-omission unit-charge local lemma even at equal heights. The8/7 value is a necessary coefficient from common-row witnesses, not a proved upper lemma or new dimension bound; no global geometric realizability of these local rows is asserted.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
