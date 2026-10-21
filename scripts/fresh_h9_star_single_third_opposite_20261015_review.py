#!/usr/bin/env python3
"""Independent literal-domain and Fraction polynomial review; no producer imports."""
import argparse, ast, json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import comb
from pathlib import Path
import fresh_h9_star_lastblock_20261014_review as R
import fresh_h9_star_single_third_20261014_review as P
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_star_single_third_opposite_20261015'
S=[3,4,6,7,9,10,12,13]; CORE=sorted(S+[11,14,17])
Z,O=R.Z,R.O
k=P.k;pa=P.pa;ps=P.ps;pm=P.pm;pp=P.pp

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(v):return list(map(str,v))
def dec(v):return tuple(map(F,v))
def parse(text,variables):
    one={(0,0,0):O}
    def go(t):
        if isinstance(t,ast.Name):return variables[t.id]
        if isinstance(t,ast.Constant):return ps(one,k(t.value))
        if isinstance(t,ast.UnaryOp):assert isinstance(t.op,ast.USub);return ps(go(t.operand),k(-1))
        if isinstance(t,ast.Call):
            assert t.func.id=='sqrt' and t.args[0].value in(2,5,10)
            v=[F(0)]*4;v[{2:1,5:2,10:3}[t.args[0].value]]=F(1);return{(0,0,0):tuple(v)}
        assert isinstance(t,ast.BinOp)
        if isinstance(t.op,ast.Add):return pa(go(t.left),go(t.right))
        if isinstance(t.op,ast.Sub):return pa(go(t.left),ps(go(t.right),k(-1)))
        if isinstance(t.op,ast.Mult):return pm(go(t.left),go(t.right))
        if isinstance(t.op,ast.Pow):return pp(go(t.left),t.right.value)
        assert isinstance(t.op,ast.Div);den=go(t.right);assert list(den)==[(0,0,0)]
        return ps(go(t.left),R.inverse(den[(0,0,0)]))
    return go(ast.parse(text,mode='eval').body)

def divide_square(poly,axis):
    # Divide by (1-variable^2)^2, using an independent leading-term loop.
    n=(0,0,0);e2=list(n);e2[axis]=2;e4=list(n);e4[axis]=4
    divisor={n:O,tuple(e2):k(-2),tuple(e4):O};out={};rem=dict(poly)
    while rem:
        e=max(rem,key=lambda z:(z[axis],z));assert e[axis]>=4
        m=list(e);m[axis]-=4;m=tuple(m);v=rem[e]
        out=pa(out,{m:v});rem=pa(rem,ps(pm({m:O},divisor),R.neg(v)))
    return out

def compose(poly,subs):
    out={}
    for ex,v in poly.items():
        term={(0,0,0):v}
        for i,n in enumerate(ex):term=pm(term,pp(subs[i],n))
        out=pa(out,term)
    return out

def run():
    seen={}
    def read(p):
        if p in seen:return seen[p]
        d=json.loads((ROOT/p).read_text());seen[p]=d
        for f,h in d.get('sources',{}).items():
            assert sha(f)==h,f
            if f.endswith('.json'):read(f)
        return d
    d=read(f'results/{TAG}.json');assert d['status']=='PASS_ALL_LABEL_SINGLE_THIRD_OPPOSITE_ANCHOR_CORE'
    g=read('results/fresh_h9_tree_cube_wedges_20261014.json')['representatives'][0]['geometry']
    V=[set(v['incident'])for v in g['vertices']];MN=set(map(frozenset,g['minimal_nonfaces']))
    assert len(V)==378 and len(MN)==23 and d['rank_core']==CORE
    assert any(set(S)<=v for v in V) and all(any(set(p)<=v for v in V)for p in combinations(CORE,2))
    faces={frozenset(v&set(CORE))for v in V};maximal={v for v in faces if not any(v<w for w in faces)}
    assert len(maximal)==9 and maximal==set(map(frozenset,d['complete_maximal_core_faces']))
    base={(3,4):4,(6,7):7,(9,10):7,(12,13):7,(3,14):2,(4,14):5,
          (11,12):2,(11,13):3,(12,17):3,(13,17):2}
    expected=[]
    def expect(kind,nodes,edges,conditions):
        expected.append((kind,tuple(nodes),tuple(tuple(x)for x in edges),json.dumps(conditions,sort_keys=True)))
    def path(kind,nodes,labels,conditions):expect(kind,nodes,[(i,i+1,m)for i,m in enumerate(labels)],conditions)
    def star(kind,nodes,conditions):expect(kind,nodes,[(0,1,4),(0,2,3),(0,3,3)],conditions)
    def fork(kind,nodes,conditions):expect(kind,nodes,[(0,1,3),(0,2,3),(0,3,3),(3,4,4)],conditions)
    for i,j in combinations(range(4),2):
        for a,b in product(S[2*i:2*i+2],S[2*j:2*j+2]):
            aa=next(x for x in S[2*i:2*i+2]if x!=a);bb=next(x for x in S[2*j:2*j+2]if x!=b)
            path('cross_block_right',[aa,a,b,bb],[4,3,4],dict(trial_nonright=[a,b]))
    for a in(11,17):
        for pair in((6,7),(9,10)):
            for r in pair:path('high_middle_anchor_right',[next(x for x in pair if x!=r),r,a],[6,3],dict(trial_nonright=[r,a]))
        for r in(3,4):path('owner_foreign_requires_internal4or5',[7-r,r,a],[6,3],dict(trial_nonright=[r,a],trial_owner_internal_at_least=6))
        path('third_anchor_label_at_most_three',[4,14,a],[4,4],dict(trial_label_at_least_four=[14,a]))
        expect('nonright_owner_foreign_forces_own_third_right',[4,14,a],[(0,1,3),(1,2,3),(0,2,3)],dict(owner_foreign=[4,a],trial_nonright=[14,a]))
    path('mutual_anchor_label_at_most_five',[13,11,17],[3,6],dict(trial_mutual_label_at_least=6))
    for r in(3,4):star('same_owner_two_foreign_forbidden',[r,7-r,11,17],dict(foreign_contacts=[[r,11],[r,17]]))
    for r,t in((3,4),(4,3)):path('different_owner_two_foreign_forbidden',[13,11,r,t,17],[3,3,4,3],dict(foreign_contacts=[[r,11],[t,17]]))
    star('both_third_contacts_forbidden',[14,4,11,17],dict(third_contacts=[[14,11],[14,17]]))
    own={11:13,17:12}
    for a,b in((11,17),(17,11)):
        fork('mutual_and_third_contact_forbidden',[a,own[a],b,14,4],dict(third_contact=[14,a],trial_mutual_nonright=[11,17]))
        for r in(3,4):fork('mutual_and_owner_foreign_forbidden',[a,own[a],b,r,7-r],dict(foreign_contact=[r,a],trial_mutual_nonright=[11,17]))
        path('nonright_owner_foreign_forces_other_third_right',[a,4,14,b,own[b]],[3,4,3,3],dict(foreign_contact=[4,a],trial_third_nonright=[14,b]))
    rows=d['positive_kernel_guards'];saved=[(r['kind'],tuple(r['ordered_support']),tuple(tuple(x)for x in r['edges']),json.dumps(r['conditional'],sort_keys=True))for r in rows]
    assert Counter(expected)==Counter(saved)and len(rows)==54
    guardcounts=Counter()
    for row in rows:
        nodes=row['ordered_support'];assert len(nodes)==len(set(nodes))
        assert set(nodes)==set(row['support'])<=V[row['actual_vertex']]and set(nodes)<=set(CORE)
        facts=dict(base)
        for name,value in row['conditional'].items():
            if name=='trial_owner_internal_at_least':facts[3,4]=value
            elif name=='trial_mutual_label_at_least':facts[11,17]=value
            elif name in('foreign_contacts','third_contacts'):
                for pair in value:facts[tuple(sorted(pair))]=3
            else:
                assert name in('trial_nonright','trial_label_at_least_four','owner_foreign','third_contact','trial_mutual_nonright','foreign_contact','trial_third_nonright')
                facts[tuple(sorted(value))]=4 if name=='trial_label_at_least_four' else 3
        edges={tuple(sorted((nodes[i],nodes[j]))):m for i,j,m in row['edges']}
        assert all(facts.get(p,0)>=m for p,m in edges.items())
        rad=3 if any(m==6 for m in edges.values()) else 2
        def qa(a,b):return a[0]+b[0],a[1]+b[1]
        def qm(a,b):return a[0]*b[0]+rad*a[1]*b[1],a[0]*b[1]+a[1]*b[0]
        def qs(a,v):return a[0]*v,a[1]*v
        def qp(text):
            def f(t):
                if isinstance(t,ast.Constant):return F(t.value),F(0)
                if isinstance(t,ast.UnaryOp):assert isinstance(t.op,ast.USub);return qs(f(t.operand),-1)
                if isinstance(t,ast.Call):assert t.func.id=='sqrt'and t.args[0].value==rad;return F(0),F(1)
                if isinstance(t.op,ast.Add):return qa(f(t.left),f(t.right))
                if isinstance(t.op,ast.Sub):return qa(f(t.left),qs(f(t.right),-1))
                if isinstance(t.op,ast.Mult):return qm(f(t.left),f(t.right))
                assert isinstance(t.op,ast.Div)and isinstance(t.right,ast.Constant);return qs(f(t.left),F(1,t.right.value))
            return f(ast.parse(text,mode='eval').body)
        co={3:(F(1,2),F(0)),4:(F(0),F(1,2)),6:(F(0),F(1,2))}
        G=[[(F(1),F(0))if i==j else qs(co[edges[tuple(sorted((a,b)))]],-1)if tuple(sorted((a,b)))in edges else(F(0),F(0))for j,b in enumerate(nodes)]for i,a in enumerate(nodes)]
        assert G==[[qp(x)for x in r]for r in row['gram']]
        w=list(map(qp,row['positive_kernel']));assert all(R.sign_quad(*x,rad)>0 for x in w)
        for r in G:
            total=(F(0),F(0))
            for a,b in zip(r,w):total=qa(total,qm(a,b))
            assert total==(0,0)
        guardcounts[row['kind']]+=1
    budgets={11:(12,13),17:(13,12),14:(3,4)}
    for row in d['projection_budgets']:
        leaf=row['leaf'];drop,keep=budgets.pop(leaf)
        assert(row['removed_own_root'],row['retained_own_root'])==(drop,keep)
        assert set(row['support'])==(set(S)-{drop})|{leaf}and len(row['support'])==8
        assert set(row['support'])<=V[row['actual_vertex']]
    assert not budgets and R.sign(R.sub(R.mul(R.COS[5],R.COS[5]),k(F(1,2))))>0
    # Budgets: anchors have own cost1/4, so foreign contacts are at most one
    # pi3 in a block4/5. Third has own cost>1/2, so any foreign pi3 costs>=1/2.
    assert F(1,4)+2*F(1,2)>1 and F(1,4)+4*F(1,4)>1
    patterns=[]
    for a,b in product((None,3,4),repeat=2):
        for dx,dy,mu in product((0,1), (0,1), (2,3,4,5)):
            if a is not None and b is not None:continue
            if dx and dy:continue
            if mu!=2 and (a is not None or b is not None or dx or dy):continue
            if (a==4 or b==4)and(dx or dy):continue
            for m in([None]if a is None and b is None else[4,5]):patterns.append([a,b,2+dx,2+dy,mu,m])
    assert len(patterns)==22 and sorted(patterns,key=str)==sorted(d['literal_patterns'],key=str)
    canonical=set()
    for a,b,dx,dy,mu,m in patterns:
        if b is not None:a,b,dx,dy=b,a,dy,dx
        canonical.add((None if a is None else a-3,dx-2,dy-2,mu,m))
    assert len(canonical)==14 and canonical==set(map(tuple,d['canonical_patterns']))
    # Exact cosine7 bounds by a separate rational polynomial recurrence.
    def ra(a,b):
        r=[F(0)]*max(len(a),len(b))
        for i,x in enumerate(a):r[i]+=x
        for i,x in enumerate(b):r[i]+=x
        return r
    def rm(a,b):
        r=[F(0)]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):r[i+j]+=x*y
        return r
    T0,T1=[F(1)],[F(0),F(1)]
    for _ in range(2,8):T0,T1=T1,ra(rm([F(0),F(2)],T1),[-x for x in T0])
    C=list(map(F,[-1,-4,4,8]));assert ra(T1,[-F(1)])==rm([-F(1),F(1)],rm(C,C))
    lo=F(31,50);hi=2*F(901,1000)**2-1
    vals=[sum(v*u**i for i,v in enumerate(C))for u in(lo,hi)]
    assert vals[0]<0<vals[1]and 24*F(1,2)**2+8*F(1,2)-4>0 and 48*F(1,2)+8>0
    assert(1+lo)/2==F(81,100)
    assert R.sign(R.sub(R.scale(R.add(k(2),(F(0),F(1),F(0),F(0))),F(1,4)),k(F(23,25)**2)))>0
    assert F(3,4)>F(6,7)**2
    one={(0,0,0):O};beta={(1,0,0):O};c={(0,1,0):O};t={(0,0,1):O}
    q=pa(one,ps(pp(c,2),k(-1)));v=pa(one,ps(pp(t,2),k(-1)));qv=pm(q,v)
    checks=[];count=Counter();keys=set()
    for row in d['continuous_cases']:
        foreign=row['foreign_owner_root'];dx,dy=map(F,row['third_anchor_cosines']);mu=row['mutual_anchor_label'];m=row['owner_internal'];branch=row['branch']
        pattern=(foreign,int(2*dx),int(2*dy),mu,m);assert pattern in canonical
        key=pattern+(branch,);assert key not in keys;keys.add(key)
        aa=pa(qv,ps(q,k(F(-1,4))));ab=pa(ps(qv,R.neg(R.COS[mu])),ps(pm(t,q),k(F(-1,4))))
        at=ps(qv,k(-dx));bt=ps(qv,k(-dy));cc=pm(pa(q,ps(pp(beta,2),k(-1))),v)
        aaf=aa
        if foreign is not None:
            aaf=pa(aa,ps(v,k(F(-1,4))));at=pa(at,ps(pm(pm(beta,c)if foreign==0 else beta,v),k(F(-1,2))))
        G=[[aaf,ab,at],[ab,aa,bt],[at,bt,cc]]
        polynomial=divide_square(divide_square(P.detpoly(G),1),2)
        assert polynomial==parse(row['scaled_determinant'],{'beta':beta,'c':c,'t':t})
        assert max(ex[2]for ex in polynomial)<=2
        exceptional=pattern==(0,0,1,2,4)
        expectedbranches={'third_ge6','third5_last_ge8','third5_last7'}if exceptional else{'uniform'}
        assert branch in expectedbranches
        lower,upper,target=F(9,10),F(1),1
        bs=pa({(0,0,0):R.COS[5]},beta)
        if branch=='third_ge6':bs=pa({(0,0,0):k(F(6,7))},beta)
        if branch=='third5_last_ge8':bs={(0,0,0):R.COS[5]};lower=F(23,25)
        if branch=='third5_last7':bs={(0,0,0):R.COS[5]};upper=F(901,1000);target=-1
        cs=pa({(0,0,0):R.COS[4]},c)if m is None else{(0,0,0):R.COS[m]}
        ts=pa({(0,0,0):k(lower)},ps(t,k(upper-lower)))
        result=compose(polynomial,[bs,cs,ts]);groups={}
        for(a,b,j),coef in result.items():groups.setdefault((a,b),{})[j]=coef
        assert(0,0)in groups
        bern={}
        for ex,terms in groups.items():
            for j in range(3):
                val=R.total([R.scale(terms.get(i,Z),F(comb(j,i),comb(2,i)))for i in range(j+1)])
                assert R.sign(val)==target;bern[ex+(j,)]=val
        saved={tuple(r['nonnegative_power'])+(r['bernstein_index'],):dec(r['coefficient'])for r in row['coefficients']}
        assert len(saved)==len(row['coefficients'])and saved==bern
        assert row['target_sign']==target and list(map(F,row['tail_interval']))==[lower,upper]
        count[target]+=len(bern)
        checks.append(dict(pattern=list(pattern),branch=branch,target_sign=target,
          polynomial=[dict(exponents=list(ex),coefficient=enc(vv))for ex,vv in sorted(polynomial.items())],
          bernstein=[dict(index=list(ex),coefficient=enc(vv))for ex,vv in sorted(bern.items())]))
    expectedkeys={p+(b,)for p in canonical for b in(['third_ge6','third5_last_ge8','third5_last7']if p==(0,0,1,2,4)else['uniform'])}
    assert keys==expectedkeys and len(keys)==16 and count=={1:165,-1:3}
    fixed=[[3,14,2],[11,12,2],[11,13,3],[12,17,3],[13,17,2]]
    bounds=[[3,4,4],[6,7,7],[9,10,7],[12,13,7],[4,14,5]]
    assert len(d['families'])==len(d['literal_orientations'])==4
    for i,(owner,last)in enumerate(product((False,True),repeat=2)):
        perm=list(range(18))
        if owner:perm[3],perm[4]=4,3
        if last:perm[12],perm[13]=13,12
        assert perm==d['literal_orientations'][i]
        assert{frozenset(perm[x]for x in z)for z in V}==set(map(frozenset,V))
        assert{frozenset(perm[x]for x in z)for z in MN}==MN
        for name,rows in(('fixed',fixed),('bounds',bounds)):
            assert sorted([sorted((perm[a],perm[b]))+[m]for a,b,m in rows])==d['families'][i][name]
    batch=read('results/fresh_h9_star_single_third_20261014_search_batch_013.json')
    labels={tuple(r[:2]):r[2]for r in batch['solution']['pair_labels']}
    assert any(all(labels[a,b]==m for a,b,m in f['fixed'])and all(labels[a,b]>=m for a,b,m in f['bounds'])for f in d['families'])
    assert batch['cumulative_additional_nodes']==d['covered_saved_cursor']==26951
    maps_path='results/fresh_h9_star_core_transfers_20261014.json'
    map_review_path='results/fresh_h9_star_core_transfers_20261014_review.json'
    old_path='results/fresh_h9_star_single_third_20261014_gates.json'
    maps,mr,old=read(maps_path),read(map_review_path),read(old_path)
    assert mr['sources'][maps_path]==sha(maps_path)
    assert maps['source_core']==CORE and maps['source_maximal_faces']==d['complete_maximal_core_faces']
    def family_key(f):return tuple(tuple(sorted(map(tuple,f[n])))for n in('fixed','bounds'))
    keys={family_key(f)for f in old['families']};assert len(keys)==349
    transfers=0
    for placement in maps['placements']:
        mp=dict(placement['core_map']);assert set(mp)==set(CORE)and len(set(mp.values()))==11
        target=set(mp.values());faces={frozenset(v&target)for v in V}
        mx={v for v in faces if not any(v<w for w in faces)}
        assert mx=={frozenset(mp[i]for i in z)for z in maximal}
        for row in d['positive_kernel_guards']+d['projection_budgets']:
            image={mp[x]for x in row['support']};assert any(image<=v for v in V);transfers+=1
        for f in d['families']:
            image={n:sorted([sorted((mp[a],mp[b]))+[m]for a,b,m in f[n]])for n in('fixed','bounds')}
            key=family_key(image);assert key not in keys;keys.add(key)
    assert len(maps['placements'])==12 and transfers==684 and len(keys)==397
    paths=[f'scripts/{TAG}.py',f'results/{TAG}.json',f'scripts/{TAG}_review.py',
           'scripts/fresh_h9_star_lastblock_20261014_review.py','scripts/fresh_h9_star_single_third_20261014_review.py',
           maps_path,map_review_path,old_path]
    return dict(status='PASS_INDEPENDENT_OPPOSITE_SINGLE_THIRD_DOMAINS_AND_POLYNOMIALS',sources={p:sha(p)for p in paths},
      arithmetic='No producer or SymPy imports. Separate Fraction Qsqrt2sqrt5 arithmetic, permutation determinants and polynomial long division; Qsqrt2/Qsqrt3 null-kernel checks.',
      recursive_manifests=len(seen),actual_vertices=378,core_maximal_faces=9,literal_guards=54,guard_counts=dict(guardcounts),projection_budgets=3,
      literal_patterns=22,canonical_patterns=14,polynomial_regions=16,positive_coefficients=165,negative_coefficients=3,
      rational_tail=dict(cubic_coefficients=list(map(str,C)),tests=[str(lo),str(hi)],values=list(map(str,vals)),cos7_interval=['9/10','901/1000'],cos8_lower='23/25'),
      independent_cases=checks,literal_families=4,covered_cursor=26951,
      transfer_check=dict(induced_core_maps=12,actual_witness_images=684,new_predicates=48,total_predicates=397),
      scope='Owner internal>=4; other three internals>=7; third own(2,k),k>=5; opposite singlepi3 anchors; all noncore contacts unrestricted. Every literal domain guard and all new polynomial arithmetic checked. No whole-incidence exclusion.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}_review.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['polynomial_regions'],d['positive_coefficients'],d['negative_coefficients'])
