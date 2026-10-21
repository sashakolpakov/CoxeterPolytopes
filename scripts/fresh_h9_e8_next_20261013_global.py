#!/usr/bin/env python3
"""Exact infinite full-chamber stabilizer for <-3(T-1)>+E8+<1>."""
from fresh_h9_e8_first_20261013 import *
TAG='fresh_h9_e8_next_20261013_global';OUT=ROOT/'results'/f'{TAG}.json'

def mm(A,B):return tuple(tuple(sum((A[i][k]*B[k][j]for k in range(len(B))),Z)for j in range(len(B[0])))for i in range(len(A)))
def transport(src,target):return mm(tuple(zip(*target)),inverse(tuple(zip(*src))))
def closure(frame):
    fg=tuple(tuple(pair(r,G,s)for s in frame)for r in frame);n=len(frame)
    seen={tuple(O if i==j else Z for i in range(n))for j in range(n)};todo=list(seen)
    while todo:
        v=todo.pop()
        for j in range(n):
            c=2*sum((v[k]*fg[k][j]for k in range(n)),Z)/fg[j][j];w=list(v);w[j]-=c;w=tuple(w)
            if w not in seen:seen.add(w);todo.append(w)
        assert len(seen)<=242
    assert all(tuple(-x for x in v)in seen for v in seen)
    return len(seen)//2

def run(verify=False):
    names=[PREFIX,PREFIX+'_independent',PREFIX+'_e7',PREFIX+'_e7_independent','fresh_h9_e8_next_20261013','fresh_h9_e8_next_20261013_independent','fresh_h9_e8_next_20261013_d6','fresh_h9_e8_next_20261013_d6_independent','fresh_h9_e8_next_20261013_d7','fresh_h9_e8_next_20261013_d7_independent']
    paths=[ROOT/'results'/f'{n}.json'for n in names];a,ai,e,ei,s,si,d,di,j,ji=[json.loads(p.read_text())for p in paths]
    for obj in(a,ai,e,ei,s,si,d,di,j,ji):
        assert obj['status'].startswith('PASS')
        for p,h in obj['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h
    assert(ei['profiles'],ei['complete_interval_points'],di['profiles'],di['complete_interval_points'],ji['profiles'],ji['complete_interval_points'])==(12,7,24,17,14,18)
    basis=tuple(tuple(map(K.from_json,r))for r in a['lattice_basis']);B=tuple(zip(*basis));BI=inverse(B)
    walls=[tuple(map(K.from_json,r))for r in j['combined_walls']]
    origin=(O,*(Z for _ in range(9)));old=walls[:9]
    dp=tuple(map(K.from_json,d['endpoint']));jp=tuple(map(K.from_json,j['endpoint']));ep=tuple(map(K.from_json,e['endpoint']))
    dframe=[walls[i]for i in[2,3,4,5,6,7,8,10,12]];jframe=[walls[i]for i in[1,2,3,4,5,6,7,8,13]]
    assert dframe==[tuple(map(K.from_json,r))for r in d['endpoint_incident']]
    assert jframe==[tuple(map(K.from_json,r))for r in j['endpoint_incident']]
    hperm=[0,1,2,3,5,4,6,8,7];jperm=[8,0,1,2,3,4,6,5,7]
    scale=square_root(-BETA/pair(jp,G,jp));assert scale is not None and scale>0
    H=transport(dframe+[dp],[dframe[i]for i in hperm]+[dp])
    J=transport(old+[origin],[jframe[i]for i in jperm]+[tuple(scale*x for x in jp)])
    I=tuple(tuple(O if i==j else Z for j in range(10))for i in range(10));checks=[]
    for name,M,src,target,point,topoint in [('H',H,dframe,[dframe[i]for i in hperm],dp,dp),('J',J,old,[jframe[i]for i in jperm],origin,tuple(scale*x for x in jp))]:
        MI=inverse(M);module=mm(BI,mm(M,B));modulei=mm(BI,mm(MI,B))
        assert all(x.is_integral for row in module for x in row)and all(x.is_integral for row in modulei for x in row)
        assert mm(M,M)==I and mm(tuple(zip(*M)),mm(G,M))==G
        assert all(mv(M,r)==s for r,s in zip(src,target))and mv(M,point)==topoint
        assert point[0]>0 and topoint[0]>0 and pair(point,G,point)<0
        checks.append({'name':name,'matrix':menc(M),'inverse':menc(MI),'lattice_matrix':menc(module),'lattice_inverse':menc(modulei),'permutation':hperm if name=='H'else jperm,'positive_root_scales':enc((O,)*9),'point_scale':O.to_json()if name=='H'else scale.to_json()})
    F=mm(H,J);FI=mm(J,H);tr=sum((F[i][i]for i in range(10)),Z);assert tr>10 and mm(F,FI)==I
    assert all(x.is_integral for row in mm(BI,mm(F,B))for x in row)
    counts=[closure(old),closure(dframe),closure(jframe)];assert counts==[121,33,121]
    # Check actual primitive reflective simple frames and oriented chamber paths.
    def integral(r):return all(x.is_integral for x in mv(BI,r))
    for r in set(old+dframe+jframe):
        q=pair(r,G,r);assert q>0 and integral(r)
        assert not any(integral(tuple(x/p for x in r))for p in(T,K(3)))
        assert all((2*pair(r,G,b)/q).is_integral for b in basis)
    paths_checked=[]
    for label,start,end,retained,released,new in [('INITIAL_E7',origin,ep,[0,2,3,4,5,6,7,8],1,10),('E7_TO_D6',ep,dp,[2,3,4,5,6,7,8,10],0,12),('INITIAL_D7',origin,jp,list(range(1,9)),0,13)]:
        assert all(pair(p,G,p)<0 and p[0]>0 for p in(start,end))
        assert all(pair(p,G,walls[i])==0 for p in(start,end)for i in retained)
        assert pair(start,G,walls[released])==pair(end,G,walls[new])==0
        assert pair(start,G,walls[new])<0 and pair(end,G,walls[released])<0
        assert all(pair(p,G,r)<=0 for p in(start,end)for r in walls)
        paths_checked.append({'name':label,'start':enc(start),'end':enc(end),'retained_indices':retained,'released_index':released,'new_index':new,'all_fourteen_walls_nonpositive':True})
    # The two first-edge paths identify the complete incident sectors of one
    # full W(L) chamber. Each positively oriented transport therefore preserves
    # its germ and hence that entire chamber. Trace >10 gives infinite order.
    result={'status':'PASS_FULL_E8_PLUS_ONE_NONREFLECTIVITY','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in(*paths,Path(__file__).resolve(),ROOT/'scripts/fresh_h9_glued_host_20261013_arithmetic.py',Path(square_root.__code__.co_filename).resolve())},
        'coefficient':BETA.to_json(),'lattice_basis':menc(basis),'generators':checks,'product':menc(F),'product_inverse':menc(FI),'trace':tr.to_json(),'trace_minus_dimension':(tr-10).to_json(),
        'complete_corner_mirror_counts':counts,'complete_corner_types':['E8+A1','D6+3A1','E8+A1'],'oriented_same_chamber_paths':paths_checked,
        'chosen_d7_simultaneous_roots':len(j['primitive_reflective_survivors']),'chosen_d7_simple_roots':1,
        'new_edges':3,'new_profiles':si['profiles']+di['profiles']+ji['profiles'],'new_interval_points':si['complete_interval_points']+di['complete_interval_points']+ji['complete_interval_points'],
        'all_certified_edges':5,'all_profiles':ai['safe_profiles']+ei['profiles']+si['profiles']+di['profiles']+ji['profiles'],'all_interval_points':ai['complete_interval_points']+ei['complete_interval_points']+si['complete_interval_points']+di['complete_interval_points']+ji['complete_interval_points'],
        'finite_prefix_facets':14,'finite_prefix_strict_escape':j['strict_escape'],'finite_prefix_escape_norm':j['escape_norm'],
        'theorem':'The full root reflection chamber of L=<-3(T-1)>+E8+<1> over Z[T], T^4-4T^2+2=0 at its largest real root, has infinite volume. H and J preserve L in both directions and map complete oriented finite-corner sectors along certified first-edge paths within one full chamber. Their product has trace greater than dimension ten and hence infinite order. An infinite stabilizer in the discrete arithmetic lattice forces infinite chamber volume.',
        'scope':'Full-host nonreflectivity, not merely prefix nonclosure. The shared certified PARI norm/unit interface remains an explicit dependency of the complete first-edge certificates. No identification with previously closed lattices or arbitrary real H9 Coxeter configurations is asserted.'}
    if verify:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'trace>10',tr>10,'corners',counts,'edges',result['all_certified_edges'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');aa=p.parse_args();run(aa.verify)
