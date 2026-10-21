#!/usr/bin/env python3
"""Exact oriented B3+B6 corner transport in the tripled D9 lattice."""
from fresh_h9_tripled_d9_20261013 import *

OUT=ROOT/'results'/f'{PREFIX}_symmetry.json'


def run(verify=False):
    source=ROOT/'results'/f'{PREFIX}_step3.json';target=ROOT/'results'/f'{PREFIX}_edge3.json'
    sp=json.loads(source.read_text());tp=json.loads(target.read_text())
    old=[tuple(map(K.from_json,r))for r in sp['endpoint_incident']]
    new=[tuple(map(K.from_json,r))for r in tp['endpoint_incident']]
    p=tuple(map(K.from_json,sp['endpoint']));q=tuple(map(K.from_json,tp['endpoint']))
    perm=[8,1,2,3,4,5,6,7,0]
    assert all(pair(r,G,s)==pair(new[perm[i]],G,new[perm[j]])for i,r in enumerate(old)for j,s in enumerate(old))
    scale=square_root(pair(p,G,p)/pair(q,G,q));assert scale is not None and scale>0
    columns=tuple(zip(*(old+[p])));images=tuple(zip(*([new[i]for i in perm]+[tuple(scale*x for x in q)])))
    A=mm(images,inverse(columns));AI=inverse(A)
    assert mm(mm(tuple(zip(*A)),G),A)==G and mm(A,AI)==mm(AI,A)==eye()
    assert all(mv(A,r)==new[perm[i]]for i,r in enumerate(old))
    assert mv(A,p)==tuple(scale*x for x in q)and mv(A,(O,*(Z for _ in range(9))))[0]>0
    basis=tuple(zip(*DBASIS));bi=inverse(basis);module=mm(mm(bi,A),basis);mi=mm(mm(bi,AI),basis)
    failures=[]
    for direction,M in(('forward',A),('inverse',AI)):
        for i,b in enumerate(DBASIS):
            if not lattice(mv(M,b)):failures.append({'direction':direction,'basis':i,'image':enc(mv(M,b))})
    tr=sum((A[i][i]for i in range(10)),Z)
    print('D9 failures',len(failures),'trace>10',tr>10,flush=True)
    if not failures:
        assert all(x.is_integral for M in(module,mi)for row in M for x in row)
    # Reconnect the exact same oriented sectors. Firstness itself is the
    # separately reproduced all-profile/all-unit interval theorem.
    firstp=ROOT/'results'/f'{PREFIX}_edge1.json';secondp=ROOT/'results'/f'{PREFIX}_step2.json'
    first=json.loads(firstp.read_text());second=json.loads(secondp.read_text())
    f1=INITIAL[1:]+[tuple(map(K.from_json,first['chosen_root']))]
    f2=[tuple(map(K.from_json,r))for r in second['endpoint_incident']]
    origin=(O,*(Z for _ in range(9)))
    vertices=[origin,tuple(map(K.from_json,first['endpoint'])),tuple(map(K.from_json,second['endpoint'])),p]
    frames=[INITIAL,f1,f2,old];paths=[]
    for start,sv,end,ev in[*[(frames[i],vertices[i],frames[i+1],vertices[i+1])for i in range(3)],(INITIAL,origin,new,q)]:
        edge=set(start)&set(end);arrive=set(end)-set(start);release=set(start)-set(end)
        assert len(edge)==8 and len(arrive)==len(release)==1
        assert all(pair(sv,G,r)==0 for r in start)and all(pair(ev,G,r)==0 for r in end)
        assert pair(sv,G,next(iter(arrive)))<0 and pair(ev,G,next(iter(release)))<0 and pair(sv,G,ev)<0
        paths.append({'start':enc(sv),'end':enc(ev),'edge_normals':menc(sorted(edge,key=lambda r:tuple(x.coeffs for x in r))),
                      'arriving':enc(next(iter(arrive))),'released':enc(next(iter(release)))})
    audits=[ROOT/'results'/f'{PREFIX}_independent.json',ROOT/'results'/f'{PREFIX}_chain_independent.json']
    ia,ca=[json.loads(x.read_text())for x in audits]
    assert ia['status']=='PASS_INDEPENDENT_TRIPLED_D9_COMPLETE_INITIAL_EDGES'
    assert ca['status']=='PASS_INDEPENDENT_TRIPLED_D9_CHAIN_EDGES'
    for audit in(ia,ca):
        for path,h in audit['sources'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h
    mirror_counts=[]
    for roots in(old,new):
        seen=set(roots);todo=list(roots)
        while todo:
            s=todo.pop()
            for r in roots:
                c=2*pair(s,G,r)/pair(r,G,r);v=tuple(x-c*y for x,y in zip(s,r))
                if v not in seen:seen.add(v);todo.append(v)
            assert len(seen)<=90
        assert len(seen)==90 and all(tuple(-x for x in r)in seen for r in seen)
        mirror_counts.append(len(seen)//2)
    result={'status':'PASS_EXACT_TRIPLED_D9_INFINITE_CHAMBER_TRANSPORT'if not failures and tr>10 else'EXACT_TRIPLED_D9_TRANSPORT_TEST_ONLY',
        'sources':{str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest()for x in(source,target,firstp,secondp,STRUCTURE,*audits)},
        'coefficient':BETA.to_json(),'permutation':perm,'positive_point_scale':scale.to_json(),
        'matrix':menc(A),'inverse':menc(AI),'module_matrix':menc(module),'module_inverse':menc(mi),
        'D9_failures':failures,'trace':tr.to_json(),'trace_exceeds10':tr>10,'oriented_paths':paths,
        'full_corner_types':['B3+B6','B3+B6'],'finite_mirror_counts':mirror_counts,
        'scope':'An integral inverse and forward map plus oriented complete simple-root transport between corners joined by four independently certified first edges preserve the full reflection chamber. Trace>10 then proves infinite order and infinite chamber volume. This implication is used only when D9_failures is empty and trace_exceeds10 is true; no rational-space transfer is used.'}
    if verify:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();run(args.verify)
