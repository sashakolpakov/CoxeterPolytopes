#!/usr/bin/env python3
"""Exact test of the two natural corner symmetries; no global exclusion."""
import argparse,hashlib,json
import fresh_h9_new_host_20261013_independent as a
from fresh_h9_glued_host_20261013_structure import norm

SOURCE=a.ROOT/'results/fresh_h9_new_host_20261013_first.json'
OUT=a.ROOT/'results/fresh_h9_new_host_20261013_symmetry.json'

def mm(A,B):return tuple(tuple(a.dot(row,col)for col in zip(*B))for row in A)
def identity(n):return tuple(tuple(a.O if i==j else a.Z for j in range(n))for i in range(n))
def enc(M):return [[a.enc(v)for v in row]for row in M]

def run(verify=False):
    data=json.loads(SOURCE.read_text());host=a.Host(1);I=identity(6)
    J=tuple(tuple(a.O if j==([0,4,3,2,1,5][i])else a.Z for j in range(6))for i in range(6))
    edge=[a.vec(r)for r in data['edge_roots']];root=a.vec(data['selected_choice']['vector']);incident=edge+[root]
    delta=tuple(a.sub(x,y)for x,y in zip(edge[1],edge[2]));cov=a.mv(host.g,delta);assert host.pair(delta,delta)==a.scale(a.O,6)
    H=tuple(tuple(a.sub(I[i][j],a.scale(a.mul(delta[i],cov[j]),a.Q(1,3)))for j in range(6))for i in range(6))
    assert mm(H,H)==I and mm(J,J)==I
    assert mm(mm(tuple(zip(*J)),host.g),J)==host.g and mm(mm(tuple(zip(*H)),host.g),H)==host.g
    start=a.vec(data['start']);end=a.vec(data['end']);assert a.mv(J,start)==start and a.mv(H,start)==start and a.mv(H,end)==end
    assert [a.mv(H,r)for r in incident]==[incident[i]for i in(0,2,1,3,4)]
    assert all(a.integral(z)for row in J for z in row)
    bad=[{'row':i,'column':j,'entry':a.enc(z)}for i,row in enumerate(H)for j,z in enumerate(row)if not a.integral(z)]
    assert bad and all(a.Q(3)*z==int(a.Q(3)*z)for row in H for v in row for z in v)
    product=mm(J,H);assert a.mv(product,start)==start
    trace=a.Z
    for i in range(6):trace=a.add(trace,product[i][i])
    assert a.sign(a.sub(a.scale(a.O,6),trace))>=0
    assert norm(a.PI)==2 # pi cannot be a K-square, since its rational norm is2.
    result={'status':'SCOPED_CORNER_SYMMETRY_TEST_COMPLETE','source':str(SOURCE.relative_to(a.ROOT)),
        'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'initial_F4_reversal':enc(J),'initial_integral_with_integral_inverse':True,
        'endpoint_A2_swap':enc(H),'endpoint_A2_swap_inverse':enc(H),
        'endpoint_A2_oriented_permutation':[0,2,1,3,4],'nonintegral_entries':bad,
        'product':enc(product),'product_trace':a.enc(trace),'both_fix_initial_time_axis':True,
        'I2_swap_obstruction':'The I2(16) simple norms are2 and2/pi. A K-isometry swapping their oriented rays would require a scalar whose square is pi or1/pi. Norm(pi)=2 proves no such K-scalar exists.',
        'scope':'The natural new-endpoint A2 diagram swap is a rational form isometry but not an integral lattice isometry; its inverse has the same denominator3. Together with initial F4 reversal it fixes the initial timelike axis, so this tested pair does not yield a loxodromic stabilizer. This is not a theorem excluding all possible chamber symmetries or proving reflectivity/nonreflectivity.'}
    if verify:assert json.loads(OUT.read_text())==result
    else:OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],'nonintegral entries',len(bad),'trace',a.enc(trace),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
