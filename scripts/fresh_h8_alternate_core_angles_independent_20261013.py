#!/usr/bin/env python3
"""Complete finite-angle audit on the 54-node compact-cap obstruction core.

Pair/sign arithmetic uses Fraction only. SymPy is used solely to construct
real-cyclotomic polynomials and rational isolating intervals by exact algebra.
No numeric proposal filter or cap-producer code is imported.
"""
from fractions import Fraction as F
from itertools import combinations
import argparse
import json
from pathlib import Path
import sympy as sp
from scripts.fresh_h4_blocks_independent_verify import (
    Z,O,P,parse,enc,add,sub,neg,scale,mul,div,sign,total,dot,mv,matrix)

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'results/fresh_h8_alternate_core_angles_independent_20261013.json'

def qa(x,y):return add(x[0],y[0]),add(x[1],y[1])
def qn(x):return neg(x[0]),neg(x[1])
def qs(x,y):return qa(x,qn(y))
def qm(x,y,d):return add(mul(x[0],y[0]),mul(d,mul(x[1],y[1]))),add(mul(x[0],y[1]),mul(x[1],y[0]))
def qsign(x,d):
    if d==Z:return sign(x[0])
    a,b=x;sa,sb=sign(a),sign(b)
    if not sb:return sa
    if not sa or sa==sb:return sb
    return sa*sign(sub(mul(a,a),mul(d,mul(b,b))))
def bsign(x,d,e):
    a=(x[0],x[1]);b=(x[2],x[3]);sa=qsign(a,d)
    if e==Z:return sa
    sb=qsign(b,d)
    if not sb:return sa
    if not sa or sa==sb:return sb
    return sa*qsign(qs(qm(a,a,d),tuple(mul(e,t) for t in qm(b,b,d))),d)
def bmul(x,y,d,e):
    out=[Z]*4
    for i,a in enumerate(x):
        if a==Z:continue
        for j,b in enumerate(y):
            if b==Z:continue
            t=mul(a,b);overlap=i&j
            if overlap&1:t=mul(t,d)
            if overlap&2:t=mul(t,e)
            out[i^j]=add(out[i^j],t)
    return tuple(out)
def bsubscalar(x,s):return sub(x[0],s),*x[1:]
def totient(m):
    n=m;result=m;p=2
    while p*p<=n:
        if n%p==0:
            result-=result//p
            while n%p==0:n//=p
        p+=1
    if n>1:result-=result//n
    return result

def angle_polynomials():
    # phi(m)^2/m >=1/2, proved prime-by-prime, gives m<=512 for
    # phi(m)<=16. The real cyclotomic degree must divide a power of two <=8.
    labels=[m for m in range(2,513) if m==2 or totient(m) in (2,4,8,16)]
    assert labels==[2,3,4,5,6,8,10,12,15,16,17,20,24,30,32,34,40,48,60]
    z,v=sp.symbols('z v');out=[]
    for m in labels:
        poly=sp.Poly(v,v) if m==2 else sp.Poly(sp.resultant(sp.cyclotomic_poly(m,z),z*z-(v-2)*z+1,z),v).sqf_part().monic()
        intervals=poly.intervals(eps=sp.Rational(1,10**24));lo,hi=intervals[-1][0]
        assert poly.degree()==(1 if m==2 else totient(m)//2)
        assert sum(multiplicity for _,multiplicity in intervals)==poly.degree()
        out.append((m,tuple(F(str(c)) for c in poly.all_coeffs()),F(str(lo)),F(str(hi))))
    return out

def finite_label(value,d,e,polys):
    for label,coefficients,lo,hi in polys:
        if bsign(bsubscalar(value,scale(O,lo)),d,e)<0:continue
        if bsign(bsubscalar(value,scale(O,hi)),d,e)>0:continue
        out=(Z,Z,Z,Z)
        for c in coefficients:
            out=bmul(out,value,d,e);out=(add(out[0],scale(O,c)),*out[1:])
        if bsign(out,d,e)==0:return label
    return None

def load():
    graph=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_selection_20261013.json').read_text())
    core=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_graph_core_20261013.json').read_text())
    base=json.loads((ROOT/graph['base_source']).read_text());polar=json.loads((ROOT/graph['polar_source']).read_text());caps=json.loads((ROOT/graph['coupled_source']).read_text())
    G=matrix(base['seed_gram']);nodes=[]
    for src in graph['nodes']:
        kind=src['source'];i=src['index']
        if kind=='base':
            row=base['caps'][i];a=tuple(map(parse,row['seed_coefficients']));b=(Z,)*10;d=Z;q=parse(row['norm'])
        elif kind=='polar':
            row=polar['polar_caps'][i];a=tuple(map(parse,row['seed_coefficients']));b=(Z,)*10;d=Z;q=parse(row['norm'])
        else:
            assert kind=='coupled';row=caps['caps'][i];a=tuple(map(parse,row['coefficients_rational']));b=tuple(map(parse,row['coefficients_radical']));d=parse(row['discriminant']);q=scale(O,2)
        if d!=Z:assert sign(d)>0
        nodes.append({'source':src,'a':a,'b':b,'pa':mv(G,a),'pb':mv(G,b),'d':d,'q':q})
    assert len(nodes)==493 and len(core['union_nodes'])==54
    return graph,core,nodes

def pairing(a,b):return dot(a['pa'],b['a']),dot(a['pb'],b['a']),dot(a['pa'],b['b']),dot(a['pb'],b['b'])

def run(verify=False):
    graph,core,nodes=load();polys=angle_polynomials();checks=[];missing=[];counts={}
    for i,j in combinations(core['union_nodes'],2):
        a,b=nodes[i],nodes[j];d,e=a['d'],b['d'];p=pairing(a,b)
        ps=bsign(p,d,e)
        if ps>0:kind='positive';label=None
        else:
            value=tuple(mul(div(scale(O,4),mul(a['q'],b['q'])),x) for x in bmul(p,p,d,e))
            cmp=bsign(bsubscalar(value,scale(O,4)),d,e)
            if cmp>0:kind='dotted';label=None
            elif cmp==0:kind='parallel';label=None
            else:
                label=finite_label(value,d,e,polys);kind='finite' if label is not None else 'noncoxeter'
        old=graph['edges'].get(f'{i},{j}');compatible=kind in ('finite','dotted')
        if old is not None:
            assert compatible and old==(str(label) if kind=='finite' else 'dotted'),(i,j,old,kind,label)
        elif compatible:missing.append({'i':i,'j':j,'kind':kind,'label':label})
        counts[kind]=counts.get(kind,0)+1
        checks.append({'i':i,'j':j,'kind':kind,'label':label,'present_in_saved_graph':old is not None})
    result={'status':'COMPLETE_CORE_COMPACT_COMPATIBILITY_GRAPH' if not missing else 'MISSING_COMPACT_COMPATIBLE_EDGES',
        'source_graph':'data/enumerations/fresh_h8_alternate_selection_20261013.json',
        'source_core':'data/enumerations/fresh_h8_alternate_graph_core_20261013.json',
        'core_clause_indices':core['clause_indices'],'core_nodes':core['union_nodes'],
        'checked_pairs':len(checks),'class_counts':counts,'missing_compatible_edges':missing,
        'complete_safe_label_pool':[row[0] for row in polys],
        'real_cyclotomic_polynomials':[{'label':m,'coefficients_descending':list(map(str,p)),'largest_root_interval':[str(lo),str(hi)]} for m,p,lo,hi in polys],
        'pair_checks':checks,
        'scope':'Every pair among the54 nodes covering the five exact witness clauses has been classified using all possible finite Coxeter angles in a field of degree2,4 or8. Parallel pairs are excluded only for compact completion. No height, proposal, or full-group completeness claim.'}
    if verify:assert json.loads(OUTPUT.read_text())==result
    else:OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k in ('status','checked_pairs','class_counts','missing_compatible_edges')}))
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
