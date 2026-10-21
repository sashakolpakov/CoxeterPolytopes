#!/usr/bin/env python3
"""Independent exact reflection word for the transported McLeod matrix.

Standalone Fraction pairs; no producer or project arithmetic imports.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = "fresh_h9_beta11_mcleod_word_20261014"
Z, O = (F(0), F(0)), (F(1), F(0))
WORD = [2,0,3,2,4,3,1,5,4,3,2,0,6,5,4,3,1,2,3,4,5,7,6,5,4,3,1,2,0,3,2,4,3,1,5,4,3,2,0,8,9,8,9,8,9,8,9,11,8,9,8,11,8,9,0,2,3,1,4,3,2,0,5,4,3,1,2,3,4,5,6,5,4,3,1,2,0,3,2,4,3,1,5,4,3,2,0,7,6,5,4,3,2,8]


def add(x,y): return x[0]+y[0],x[1]+y[1]
def neg(x): return -x[0],-x[1]
def sub(x,y): return add(x,neg(y))
def scale(x,q): return q*x[0],q*x[1]
def mul(x,y): return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def conj(x): return x[0]+x[1],-x[1]
def norm(x): return x[0]*x[0]+x[0]*x[1]-x[1]*x[1]
def inv(x):
    n=norm(x);assert n
    return scale(conj(x),1/n)
def sign(x):
    a,b=2*x[0]+x[1],x[1]
    sg=lambda q:(q>0)-(q<0)
    if not a:return sg(b)
    if not b or sg(a)==sg(b):return sg(a)
    return sg(a)*sg(a*a-5*b*b)
def integral(x):return x[0].denominator==x[1].denominator==1
def total(xs):
    out=Z
    for x in xs:out=add(out,x)
    return out
def dot(x,y):return total(mul(a,b)for a,b in zip(x,y))
def transpose(A):return list(map(list,zip(*A)))
def mm(A,B):return [[dot(row,col)for col in zip(*B)]for row in A]
def mv(A,v):return[dot(row,v)for row in A]
def eye(n):return [[O if i==j else Z for j in range(n)]for i in range(n)]
def dec(x):return(F(x['a']),F(x['b']))if isinstance(x,dict)else tuple(map(F,x))
def enc(x):return list(map(str,x))
def em(A):return [[enc(x)for x in row]for row in A]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()


def run():
    paths=['results/fresh_h9_beta11_literature_20261014.json',
           'results/fresh_h9_beta11_vinberg_20261014.json']
    seen={}
    def read(p):
        if p in seen:return seen[p]
        d=json.loads((ROOT/p).read_text());seen[p]=d
        for s,h in d.get('sources',{}).items():
            assert sha(s)==h,s
            if s.endswith('.json'):read(s)
        return d
    literature,prefix=map(read,paths)
    assert prefix['height_bound']=='100' or prefix['height_bound']==100
    assert prefix['status']=='PASS_COMPLETE_BOTH_NORM_VINBERG_PREFIX'
    G=[[dec(x)for x in row]for row in prefix['form']]
    assert G==[[dec(x)for x in row]for row in literature['gram']]
    g=[[dec(x)for x in row]for row in literature['mcleod_symmetry_proposal']['transported_matrix10']]
    assert mm(transpose(g),mm(G,g))==G
    assert all(integral(x)for row in g for x in row)
    assert sign(sub(total(g[i][i]for i in range(10)),(F(10),F(0))))>0
    roots=[[dec(x)for x in r]for r in prefix['roots']]
    H=list(map(dec,prefix['basepoint']))
    assert sign(dot(H,mv(G,H)))<0
    assert len(WORD)==94 and set(WORD)==set(range(10))|{11}
    reflections={};records=[]
    for index in sorted(set(WORD)):
        r=roots[index];gr=mv(G,r);q=dot(r,gr)
        assert q==dec(prefix['gram'][index][index])
        assert q in [(F(2),F(0)),(F(6),F(2))]
        assert sign(q)>0 and sign(conj(q))>0
        assert all(integral(x)for x in r)
        coefficients=[mul(scale(x,2),inv(q))for x in gr]
        assert all(integral(x)for x in coefficients)
        M=[[sub(O if i==j else Z,mul(r[i],coefficients[j]))for j in range(10)]for i in range(10)]
        assert all(integral(x)for row in M for x in row)
        assert mm(M,M)==eye(10)
        assert mm(transpose(M),mm(G,M))==G
        assert mv(M,r)==list(map(neg,r))
        assert total(M[i][i]for i in range(10))==(F(8),F(0))
        reflections[index]=M
        records.append(dict(index=index,root=list(map(enc,r)),norm=enc(q),
                            integral_reflection_matrix=em(M)))
    current=[row[:]for row in g]
    at54=None
    for step,index in enumerate(WORD,1):
        current=mm(reflections[index],current)
        if step==54:
            assert mv(current,H)==H
            assert current!=eye(10)
            at54=current
    assert current==eye(10)
    # Independent multiplication direction: g equals S_first ... S_last.
    rebuilt=eye(10)
    for index in WORD:rebuilt=mm(rebuilt,reflections[index])
    assert rebuilt==g
    paths += [f'scripts/{TAG}.py']
    return dict(status='PASS_MCLEOD_MATRIX_IS_LITERAL_REFLECTION_WORD',
        sources={p:sha(p)for p in paths},recursive_manifests=len(seen),
        arithmetic='Standalone Fraction pairs a+b*phi; no producer, project arithmetic or symbolic imports.',
        coordinate_order='x,E8_0,...,E8_7,t; same frozen height100 and literature basis',
        height_bound=100,word=WORD,word_length=len(WORD),
        word_convention='S[word[-1]] ... S[word[0]] g = I; equivalently g = S[word[0]] ... S[word[-1]].',
        used_indices=sorted(reflections),occurrences=dict(sorted(Counter(WORD).items())),
        used_reflections=records,transported_matrix=em(g),
        after54_residual_matrix=em(at54),after54_fixes_basepoint=True,
        after54_identity=False,after94_matrix=em(current),
        subgroup_statement='g is in the subgroup generated by the eleven listed known integral reflections, including the norm2beta wall9.',
        conclusion='This infinite-order candidate cannot be a nontrivial stabilizer of a chamber of the full integral reflection arrangement. This rejects the candidate, not reflectivity or nonreflectivity of L.',
        scope='Exact word certificate only. No new root enumeration, Vinberg height-completeness replay, or full-chamber assertion. The equality already uses known roots, so later walls do not change subgroup membership.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    # JSON object keys normalize the integer occurrence indices to strings.
    d=json.loads(json.dumps(d))
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['word_length'])
