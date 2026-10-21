#!/usr/bin/env python3
"""Exact compatibility supergraph with all possible degree-dividing-eight labels.

Each archived norm-two cap lies in K(sqrt(D)), K=Q(phi). A cap pairing
therefore lies in K(sqrt(D),sqrt(E)), whose rational degree divides eight.
For a finite Coxeter label m, [Q(2cos(pi/m)):Q]=totient(2m)/2. The
elementary exponent-budget enumeration below leaves fourteen possible
labels, including 17 conservatively. No large-label tail is needed.

The earlier all-label graph is preserved. Every rejected pair is replayed
with general four-product interval multiplication, and all coefficient
enclosures are checked using an independent exact quadratic-field sign.
This graph alone makes no SAT or finite-volume assertion.
"""
import argparse
from bisect import bisect_left
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
import time

import fresh_h9_split_leaf_3034_safe_graph_20261013 as S

ROOT=S.ROOT
OUT=ROOT/'data/enumerations/fresh_h9_split_leaf_3034_degree_graph_20261013.json'
AUDIT=ROOT/'results/fresh_h9_split_leaf_3034_degree_graph_20261013_verify.json'
LABELS=(2,3,4,5,6,8,10,12,15,16,17,20,24,30)


def degree_labels():
    """Exhaustive arithmetic enumeration, not a search up to an m cutoff.

    totient(2m) must be a power of two <=16. Any odd prime p dividing
    2m has p-1 in {2,4,8,16}, hence p in {3,5,17}; its exponent is one.
    Write 2m=2^a*3^b*5^c*17^d, with a>=1 and b,c,d in {0,1}.
    totient(2m)=2^(a-1+b+2c+4d); this exponent must lie in [1,4].
    """
    rows=[]
    for a in range(1,6):
        for b in range(2):
            for c in range(2):
                for d in range(2):
                    exponent=a-1+b+2*c+4*d
                    if not 1<=exponent<=4:continue
                    m=(2**a*3**b*5**c*17**d)//2
                    assert m>=2
                    rows.append({'label':m,'exponents_2_3_5_17':[a,b,c,d],
                                 'totient_2m':2**exponent,'real_cyclotomic_degree':2**(exponent-1)})
    rows.sort(key=lambda x:x['label'])
    assert tuple(x['label']for x in rows)==LABELS
    return rows


class Gate:
    def __init__(self,proof):
        rows=[x for x in proof['finite_angles']if x['label']in LABELS]
        assert tuple(x['label']for x in rows)==LABELS
        rows.reverse()
        self.lower=[int(x['lower'])for x in rows]
        self.upper=[int(x['upper'])for x in rows]
    def rejection(self,lo,hi):
        assert lo<=hi
        if lo>0:return'positive'
        # Preserve any interval which could be parallel or ultraparallel.
        if lo<=-2*S.PAIR_SCALE:return None
        i=bisect_left(self.upper,lo)
        if i<len(self.upper)and self.lower[i]<=hi:return None
        return'finite_angle_gap'


def setup():
    source,enc,normal,covector,angles=S.source_and_bounds()
    old=json.loads(S.OUT.read_text())
    assert old['source_sha256']==hashlib.sha256(S.SOURCE.read_bytes()).hexdigest()
    assert old['angle_proof']==angles
    assert old['coefficient_scale']==str(S.SCALE)
    assert S.digest(old['adjacency_hex'])==old['adjacency_sha256']
    oldadj=[int(x,16)for x in old['adjacency_hex']]
    assert len(oldadj)==len(normal)==3040
    return source,enc,normal,covector,angles,old,oldadj


def build():
    start=time.time()
    source,enc,normal,covector,angles,old,oldadj=setup()
    label_proof=degree_labels();gate=Gate(angles)
    n=len(normal);adj=[0]*n;bits=[1<<i for i in range(n)]
    counts=Counter();additional=[]
    print('Exact degree alphabet and coefficient enclosures ready',LABELS,flush=True)
    for i in range(n):
        ai=adj[i];a=normal[i];bit=bits[i]
        for j in range(i+1,n):
            lo,hi=S.interval_dot_fast(a,covector[j]);reason=gate.rejection(lo,hi)
            if reason is None:
                assert oldadj[i]&bits[j],(i,j,'unexpected added edge')
                ai|=bits[j];adj[j]|=bit;counts['possible_edges']+=1
            else:
                counts[reason]+=1
                if oldadj[i]&bits[j]:additional.append([i,j])
        adj[i]=ai
        if i%500==0:print('rows',i,'seconds',round(time.time()-start,1),dict(counts),flush=True)
    hexadecimal=[format(x,'x')for x in adj]
    assert sum(counts.values())==n*(n-1)//2
    assert old['pair_counts']['possible_edges']-counts['possible_edges']==len(additional)
    out={'status':'COMPLETE_EXACT_DEGREE_EIGHT_COMPATIBILITY_SUPERGRAPH',
         'entry_count':n,'exact_distinct_mirrors':3034,
         'source':str(S.SOURCE.relative_to(ROOT)),'source_sha256':old['source_sha256'],
         'earlier_graph':str(S.OUT.relative_to(ROOT)),
         'earlier_adjacency_sha256':old['adjacency_sha256'],
         'coefficient_scale':str(S.SCALE),'fine_scale':str(S.FINE),
         'phi_fine_interval':list(map(str,enc.phi)),
         'degree_enumeration':label_proof,'allowed_finite_labels':list(LABELS),
         'angle_proof':angles,'pair_counts':dict(counts),
         'additional_nonedge_count':len(additional),'additional_nonedges':additional,
         'adjacency_hex':hexadecimal,'adjacency_sha256':S.digest(hexadecimal),
         'scope':'Every compatible pair of distinct archived norm-two cap mirrors is retained. Pairings lie in a field of degree dividing eight, so the fourteen listed finite labels exhaust all possibilities; m=17 is retained conservatively. Dotted/parallel or unresolved intervals remain edges. The six duplicate profile IDs remain, but duplicate inequalities can be omitted. No clique, SAT, volume, or global polytope-existence conclusion follows from this graph alone.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n')
    print('COMPLETE',dict(counts),'additional_nonedges',len(additional),
          'seconds',round(time.time()-start,1),flush=True)


def verify():
    import fresh_h9_split_leaf_3040_independent_20261013 as Exact
    start=time.time()
    source,enc,normal,covector,angles,old,oldadj=setup()
    saved=json.loads(OUT.read_text());n=len(normal)
    assert saved['status']=='COMPLETE_EXACT_DEGREE_EIGHT_COMPATIBILITY_SUPERGRAPH'
    assert saved['entry_count']==n and saved['exact_distinct_mirrors']==3034
    assert saved['source']==str(S.SOURCE.relative_to(ROOT))
    assert saved['source_sha256']==hashlib.sha256(S.SOURCE.read_bytes()).hexdigest()
    assert saved['earlier_graph']==str(S.OUT.relative_to(ROOT))
    assert saved['earlier_adjacency_sha256']==old['adjacency_sha256']
    assert saved['coefficient_scale']==str(S.SCALE)and saved['fine_scale']==str(S.FINE)
    assert saved['phi_fine_interval']==list(map(str,enc.phi))
    assert saved['degree_enumeration']==degree_labels()
    assert saved['allowed_finite_labels']==list(LABELS)
    assert saved['angle_proof']==angles
    seen=set()
    for raw,nr,cg in zip(source['caps'],normal,covector):
        d=Exact.parse(raw['radicand']);field=Exact.Quadratic(d)
        for vals,bounds in((raw['normal'],nr),(raw['covector'],cg)):
            for value,(lo,hi)in zip(vals,bounds):
                x=tuple(map(Exact.parse,value));key=(x,d,lo,hi)
                if key in seen:continue
                seen.add(key)
                assert field.sign(field.sub(x,((F(lo,S.SCALE),F(0)),Exact.Z)))>=0
                assert field.sign(field.sub(((F(hi,S.SCALE),F(0)),Exact.Z),x))>=0
    print('Independent coefficient-bound checks PASS',len(seen),flush=True)
    hexadecimal=saved['adjacency_hex']
    assert len(hexadecimal)==n and S.digest(hexadecimal)==saved['adjacency_sha256']
    adj=[int(x,16)for x in hexadecimal];mask=(1<<n)-1
    assert all(a>=0 and not a&~mask and not(a>>i&1)for i,a in enumerate(adj))
    assert all(not(a&~b)for a,b in zip(adj,oldadj))
    counts=Counter();additional=[];gate=Gate(angles)
    for i in range(n):
        for j in range(i+1,n):
            edge=bool(adj[i]>>j&1);assert edge==bool(adj[j]>>i&1)
            if edge:counts['possible_edges']+=1;continue
            # General interval multiplication independently replays ALL nonedges,
            # including every rejection newly justified by the field degree.
            lo,hi=S.interval_dot_general(normal[i],covector[j])
            reason=gate.rejection(lo,hi);assert reason is not None,(i,j)
            counts[reason]+=1
            if oldadj[i]>>j&1:additional.append([i,j])
        if i%500==0:print('verified rows',i,'seconds',round(time.time()-start,1),flush=True)
    assert dict(counts)==saved['pair_counts']
    assert additional==saved['additional_nonedges']
    assert len(additional)==saved['additional_nonedge_count']
    assert old['pair_counts']['possible_edges']-counts['possible_edges']==len(additional)
    # Negative controls on the gate: no allowed exact angle or dotted pair may
    # be rejected, whereas a tightly enclosed degree-three label 7 is rejected.
    assert gate.rejection(0,0)is None
    assert gate.rejection(-2*S.PAIR_SCALE,-2*S.PAIR_SCALE)is None
    for row in angles['finite_angles']:
        lo,hi=int(row['lower']),int(row['upper'])
        if row['label']in LABELS:assert gate.rejection(lo,hi)is None
        else:assert gate.rejection(lo,hi)=='finite_angle_gap'
    result={'status':'PASS_DEGREE_EIGHT_SAFE_SUPERGRAPH_REPLAY','standard_library_only':True,
            'entry_count':n,'exact_distinct_mirrors':3034,
            'distinct_coefficient_enclosures_checked_by_independent_exact_field_signs':len(seen),
            'allowed_finite_labels':list(LABELS),'pair_counts':dict(counts),
            'additional_nonedge_count':len(additional),
            'nonedges_checked_with_general_interval_products':counts['positive']+counts['finite_angle_gap'],
            'strict_subgraph_of_preserved_all_label_graph':len(additional)>0,
            'adjacency_sha256':saved['adjacency_sha256'],
            'all_allowed_angle_enclosures_retained':True,
            'all_excluded_labels_through_120_rejected':True,
            'scope':saved['scope']}
    AUDIT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    verify()if args.verify else build()
