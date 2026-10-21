#!/usr/bin/env python3
"""Exhaust the eight diagram symmetries of the new A4+A4+A1 corner."""
import argparse,hashlib,json
from itertools import product
from fresh_h9_new_host_next_20261013_a9_independent import a,G,P
ROOT=a.ROOT;SOURCE=ROOT/'results/fresh_h9_new_host_next_20261013_a9_endpoint.json';OUT=ROOT/'results/fresh_h9_new_host_next_20261013_a9_symmetry.json'
mm=lambda A,B:tuple(tuple(a.dot(r,c)for c in zip(*B))for r in A)
I=tuple(tuple(a.O if i==j else a.Z for j in range(10))for i in range(10))
enc=lambda M:[[a.enc(v)for v in row]for row in M]

def run(verify=False):
    d=json.loads(SOURCE.read_text());walls=[a.vec(r)for r in d['all_walls']];incident=[walls[i]for i in d['actual_finite_vertex']['indices']]
    gram=tuple(tuple(P(r,s)for s in incident)for r in incident);gi=a.minverse(gram);cov=[a.mv(G,r)for r in incident]
    point=a.vec(d['endpoint']);origin=(a.O,*(a.Z,)*9)
    J=tuple(tuple(a.O if j==([0,*range(9,0,-1)][i])else a.Z for j in range(10))for i in range(10));rows=[]
    for swap,r1,r2 in product((False,True),repeat=3):
        b1=list(range(4));b2=list(range(4,8))
        if r1:b1.reverse()
        if r2:b2.reverse()
        if swap:b1,b2=b2,b1
        perm=b1+b2+[8]
        assert tuple(tuple(gram[perm[i]][perm[j]]for j in range(9))for i in range(9))==gram
        H=tuple(tuple(a.add(I[i][j],a.dot(tuple(a.sub(incident[perm[k]][i],incident[k][i])for k in range(9)),
          a.mv(gi,tuple(cov[l][j]for l in range(9)))))for j in range(10))for i in range(10))
        assert mm(mm(tuple(zip(*H)),G),H)==G and a.mv(H,point)==point and a.mv(H,origin)==origin
        assert [a.mv(H,r)for r in incident]==[incident[i]for i in perm]
        integral=all(a.integral(z)for row in H for z in row)
        if integral:assert H in(I,J)
        else:assert all(a.integral(a.scale(z,5))for row in H for z in row)
        rows.append({'permutation':perm,'matrix':enc(H),'integral':integral,'fixes_initial_time_axis':True})
    assert sum(r['integral']for r in rows)==2 and mm(J,J)==I
    result={'status':'COMPLETE_SCOPED_A9_CORNER_DIAGRAM_SYMMETRY_TEST','source':str(SOURCE.relative_to(ROOT)),
      'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'endpoint_diagram_automorphisms':rows,
      'integral_count':2,'integral_group_order':2,
      'scope':'All eight positive oriented diagram automorphisms of A4+A4+A1 fix the unique A1 cap and the timelike endpoint, hence the initial time axis. Exactly identity and the initial A9 reversal are integral; the other six have denominator5. Thus this two-corner diagram-symmetry mechanism gives only a group of order2, not an infinite stabilizer. Other corners or symmetries remain untested.'}
    if verify:assert json.loads(OUT.read_text())==result
    else:OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],result['integral_count'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
