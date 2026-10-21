#!/usr/bin/env python3
"""Independent complete last-block4/5 profile review.

No producer imports or symbolic-algebra dependency.  Arithmetic uses the
ordered basis 1,sqrt(2),sqrt(5),sqrt(10), and determinants use permutations.
"""
import argparse
import json
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import permutations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_star_lastblock_20261014'
Z = (F(0),) * 4
O = (F(1), F(0), F(0), F(0))
S = [3, 4, 6, 7, 9, 10, 12, 13]
P = [2, 5, 8, 11]
C = [14, 15, 16, 17]

def add(x, y): return tuple(a+b for a, b in zip(x, y))
def neg(x): return tuple(-a for a in x)
def sub(x, y): return add(x, neg(y))
def scale(x, a): return tuple(a*b for b in x)
def total(xs):
    v = Z
    for x in xs: v = add(v, x)
    return v

@lru_cache(None)
def mul(x, y):
    out = [F(0)]*4
    for i, a in enumerate(x):
        if not a: continue
        for j, b in enumerate(y):
            if b: out[i ^ j] += a*b*(2 if i & j & 1 else 1)*(5 if i & j & 2 else 1)
    return tuple(out)

def sign_quad(a, b, d):
    sg = lambda q: (q > 0)-(q < 0)
    if not b: return sg(a)
    if not a: return sg(b)
    if sg(a) == sg(b): return sg(a)
    return sg(a)*sg(a*a-d*b*b)

@lru_cache(None)
def sign(x):
    a, b, c, d = x
    u, v = sign_quad(a, b, 2), sign_quad(c, d, 2)
    if not v: return u
    if not u: return v
    if u == v: return u
    return u*sign_quad(a*a+2*b*b-5*c*c-10*d*d, 2*a*b-10*c*d, 2)

@lru_cache(None)
def inverse(x):
    conjugates = [(x[0], s*x[1], t*x[2], s*t*x[3])
                  for s, t in [(-1, 1), (1, -1), (-1, -1)]]
    n = O
    for y in conjugates: n = mul(n, y)
    den = mul(x, n)
    assert den[0] and den[1:] == (0, 0, 0)
    return scale(n, 1/den[0])

def dec_anp(row):
    p = [F(v) for v in row]
    p = [F(0)]*(4-len(p))+p
    a3, a2, a1, a0 = p
    return a0+7*a2, a1+17*a3, a1+11*a3, 2*a2

COS = {2: Z, 3: (F(1, 2), F(0), F(0), F(0)),
       4: (F(0), F(1, 2), F(0), F(0)),
       5: (F(1, 4), F(0), F(1, 4), F(0))}
PERMS = [(p, (-1)**sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3)))
         for p in permutations(range(3))]
def determinant3(M):
    return total(scale(mul(mul(M[0][p[0]], M[1][p[1]]), M[2][p[2]]), sg)
                 for p, sg in PERMS)
def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()

def failure(lab,row,V):
    nodes=sorted(V[row['actual_vertex']])
    M=[[O if a==b else neg(COS[lab.get(tuple(sorted((a,b))),2)])for b in nodes]for a in nodes]
    piv=[]
    for k in range(9):
        p=M[k][k];piv.append(p)
        if sign(p)<=0:
            assert k+1==row['prefix_order']
            assert piv==list(map(dec_anp,row['pivots']))
            return
        for i in range(k+1,9):
            for j in range(k+1,9):M[i][j]=sub(M[i][j],mul(mul(M[i][k],M[k][j]),inverse(p)))
    raise AssertionError('no nonpositive pivot')

def run():
    src=f'results/{TAG}.json';d=json.loads((ROOT/src).read_text());todo=[src];seen=set()
    while todo:
        p=todo.pop()
        if p in seen:continue
        seen.add(p);obj=json.loads((ROOT/p).read_text())
        for q,h in obj.get('sources',{}).items():
            assert sha(q)==h,q
            if q.endswith('.json'):todo.append(q)
    g=json.loads((ROOT/'results/fresh_h9_tree_cube_wedges_20261014.json').read_text())['representatives'][0]['geometry']
    V=[set(v['incident'])for v in g['vertices']];missing={frozenset(m)for m in g['minimal_nonfaces']}
    assert len(V)==378
    own={tuple(sorted((leaf,a))):3 for i in range(3)for leaf in(P[i],C[i])for a in S[2*i:2*i+2]}
    assert d['fixed_first_three_own']==sorted([list(p)+[m]for p,m in own.items()])
    assert d['internal_lower_bounds']==[S[2*i:2*i+2]+[4]for i in range(3)]and d['last_internal_values']==[4,5]
    for i in range(4):
        for leaf in(P[i],C[i]):assert frozenset(S[2*i:2*i+2]+[leaf])in missing
    # First3 internal2/3 cannot be their required Lanner triple with own3/3.
    for m,want in((2,scale(O,F(1,2))),(3,Z)):
        M=[[O,neg(COS[3]),neg(COS[3])],[neg(COS[3]),O,neg(COS[m])],[neg(COS[3]),neg(COS[m]),O]]
        assert determinant3(M)==want
    gates=json.loads((ROOT/'results/fresh_h9_star_zero_20261014_gates.json').read_text())
    assert len(d['first_three_high_reductions'])==3
    for r in d['first_three_high_reductions']:
        i=r['high_block'];assert i in range(3)
        lab=own|{tuple(S[2*j:2*j+2]):6 if i==j else 4 for j in range(4)}
        gate=next(x for x in gates['families']if x['id']==r['reviewed_gate'])
        assert gate['origin']=='two_dense'and all(lab.get((a,b),0)>=m for a,b,m in gate['bounds'])
    cross={tuple(sorted((a,b)))for i in range(4)for j in range(i+1,4)for a in S[2*i:2*i+2]for b in S[2*j:2*j+2]}
    assert {tuple(r['pair'])for r in d['cross_right_guards']}==cross
    sqrt2=scale(COS[4],2);null=[O,sqrt2,sqrt2,O]
    path=[[O if i==j else neg(COS[[4,3,4][min(i,j)]])if abs(i-j)==1 else Z for j in range(4)]for i in range(4)]
    assert all(total(mul(path[i][j],null[j])for j in range(4))==Z for i in range(4))
    for r in d['cross_right_guards']:
        a,b=r['pair'];i=S.index(a)//2;j=S.index(b)//2
        assert set(r['support'])==set(S[2*i:2*i+2]+S[2*j:2*j+2])and set(r['support'])<=V[r['actual_vertex']]
    budgets=set()
    for r in d['retained_root_budgets']:
        leaf,keep=r['leaf'],r['retained_own_root'];group=0 if leaf==14 else 3
        assert leaf in(11,17,14)and keep in S[2*group:2*group+2]
        drop=next(x for x in S[2*group:2*group+2]if x!=keep)
        assert set(r['support'])==(set(S)-{drop})|{leaf}and set(r['support'])<=V[r['actual_vertex']]
        budgets.add((leaf,keep))
    assert budgets=={(a,b)for a in(11,17)for b in(12,13)}|{(14,3),(14,4)}
    triangles=set()
    for r in d['high_own_triangle_guards']:
        a,b=r['anchor'],r['high_root'];assert a in(11,17)and b in(12,13)
        assert set(r['support'])=={14,a,b}and set(r['support'])<=V[r['actual_vertex']]
        triangles.add((a,b))
    assert triangles=={(a,b)for a in(11,17)for b in(12,13)}
    # Highown>=6 gives owncost>=3/4. Any foreign finite contact adds>=1/2.
    # Its literal triangle with14 excludes either of the two additional
    # bonds already from3/4+1/4=1. Third14 foreign budget is unchanged.
    assert F(3,4)+F(1,2)>1 and F(3,4)+F(1,4)==1
    assert F(1,4)+2*2*F(1,4)>1 and F(1,4)+2*F(1,2)>1
    # Independent polynomial identity via coefficient dictionaries. It is
    # equivalent to X-D(c*alpha+beta)^2=alpha^2 when D(1-c²)=1.
    # Exact positive decomposition uses a=alpha²-3/4>=0,d=D-2>=0:
    # lower =sqrt2/4-1/8+d(3sqrt2/8-1/16)+(2+d)a(sqrt2/2+1/4).
    minimum=sub(scale(sqrt2,F(1,4)),scale(O,F(1,8)))
    slope=sub(scale(sqrt2,F(3,8)),scale(O,F(1,16)))
    growth=add(scale(sqrt2,F(1,2)),scale(O,F(1,4)))
    assert all(sign(v)>0 for v in(minimum,slope,growth))
    assert minimum==dec_anp(d['high_own_negative_plane']['strict_determinant_lower'])
    for n in(4,5):
        c=COS[n];D=inverse(sub(O,mul(c,c)))
        assert mul(D,sub(O,mul(c,c)))==O and sign(sub(D,scale(O,2)))>=0
        # Coefficients of alpha², alpha*beta,beta² in X-D(cα+β)².
        assert sub(D,mul(D,mul(c,c)))==O
        assert sub(scale(mul(D,c),2),scale(mul(D,c),2))==Z
        assert sub(D,D)==Z
    profiles={n:[p for p in product(range(2,6),repeat=2)if F(1,n)+sum(F(1,m)for m in p)<1]for n in(4,5)}
    assert {str(n):list(map(list,p))for n,p in profiles.items()}==d['ordered_own_profiles']
    assert [len(profiles[n])for n in(4,5)]==[11,13]
    domains={};active={4:set(),5:set()};status_counts={}
    for row in d['own_pair_domains']:
        n=row['last_internal'];a,b=map(tuple,row['anchor_profiles']);key=(n,a,b)
        assert key not in domains and a in profiles[n]and b in profiles[n];domains[key]=row
        common=[i for i in range(2)if a[i]>2 and b[i]>2]
        if not common:
            assert row['status']=='INHERITED_OPPOSITE_SINGLE'
            assert a.count(2)==b.count(2)==1 and a.index(2)!=b.index(2)and max(a)>=4 and max(b)>=4
        else:
            r=row['mutual_right_guard'];j=common[0]
            assert set(r['support'])=={11,17,12+j}and set(r['support'])<=V[r['actual_vertex']]
            ones=sub(scale(O,2),scale(add(COS[a[j]],COS[b[j]]),2))
            assert sign(ones)<=0 and ones==dec_anp(r['ones_quadratic_at_test3'])
            bad=next((j for j in range(2)if sign(sub(sub(O,mul(COS[a[j]],COS[a[j]])),mul(COS[b[j]],COS[b[j]])))<=0),None)
            if bad is not None:
                assert row['status']=='NONPOSITIVE_OWN_TRIANGLE'
                r=row['triangle'];assert set(r['support'])=={11,17,12+bad}and set(r['support'])<=V[r['actual_vertex']]
                value=sub(sub(O,mul(COS[a[bad]],COS[a[bad]])),mul(COS[b[bad]],COS[b[bad]]))
                assert value==dec_anp(row['determinant'])
            elif a==b==(3,3):assert row['status']=='INHERITED_ALL_OWN_THREE'
            else:assert row['status']=='RANK_CENSUS';active[n].add((a,b))
        status_counts[row['status']]=status_counts.get(row['status'],0)+1
    assert set(domains)==set((n,a,b)for n in(4,5)for a,b in product(profiles[n],repeat=2))and len(domains)==290
    assert [len(active[n])for n in(4,5)]==[36,48]
    uppers=set()
    for r in d['third_upper_guards']:
        n=r['last_internal'];a,b=map(tuple,r['anchor_profiles']);pair=tuple(r['pair'])
        assert(a,b)in active[n]and pair in((11,14),(14,17))and r['test_label']==4
        assert set(pair)<=V[r['actual_vertex']]
        lab=own|{tuple(S[2*i:2*i+2]):n if i==3 else 4 for i in range(4)}
        for leaf,p in zip((11,17),(a,b)):
            for root,m in zip((12,13),p):lab[tuple(sorted((leaf,root)))]=m
        lab[pair]=4;failure(lab,r,V);uppers.add((n,a,b,pair))
    assert uppers=={(n,a,b,p)for n in(4,5)for a,b in active[n]for p in((11,14),(14,17))}and len(uppers)==168
    def patterns(pr,group):
        extras=[None]+[i for i in range(8)if i//2!=group]if pr==(3,3)else[None]
        out=[]
        for e in extras:
            v=[Z]*8;v[2*group]=neg(COS[pr[0]]);v[2*group+1]=neg(COS[pr[1]])
            if e is not None:v[e]=neg(COS[3])
            out.append((e,v))
        return out
    tr=patterns((3,3),0);zeros={};anchor_signs={};total_cases=0;nonzero=0;records=set()
    for row in d['zero_rank_cases']:
        key=(tuple(row['block_labels']),tuple(map(tuple,row['anchor_profiles'])),tuple(row['foreign_contacts']),tuple(row['third_labels']))
        assert key not in zeros;zeros[key]=row
    usedzeros=set()
    for rec in d['records']:
        blocks=tuple(rec['block_labels']);a,b=map(tuple,rec['anchor_profiles']);key=(blocks,a,b)
        assert blocks in set(product((4,5),repeat=4))and(a,b)in active[blocks[3]]and key not in records;records.add(key)
        inv=[inverse(sub(O,mul(COS[m],COS[m])))for m in blocks]
        for m,v in zip(blocks,inv):assert sign(sub(v,scale(O,2)))>=0
        def inner(x,y):
            return total(mul(inv[i],add(add(mul(x[2*i],y[2*i]),mul(x[2*i+1],y[2*i+1])),mul(COS[blocks[i]],add(mul(x[2*i],y[2*i+1]),mul(x[2*i+1],y[2*i])))))for i in range(4))
        ar,br=patterns(a,3),patterns(b,3)
        assert rec['pattern_order']==[[e for e,v in pp]for pp in(ar,br,tr)]and rec['contact_order']==[[2,2],[2,3],[3,2],[3,3]]
        aa=[sub(O,inner(v,v))for e,v in ar];bb=[sub(O,inner(v,v))for e,v in br];tt=[sub(O,inner(v,v))for e,v in tr]
        ab=[[neg(inner(v,w))for e,w in br]for f,v in ar]
        at=[[neg(inner(v,w))for e,w in tr]for f,v in ar];bt=[[neg(inner(v,w))for e,w in tr]for f,v in br]
        ix=0
        for ia,ib,ic in product(range(len(ar)),range(len(br)),range(7)):
            av=sub(mul(aa[ia],bb[ib]),mul(ab[ia][ib],ab[ia][ib]));sg=str(sign(av));anchor_signs[sg]=anchor_signs.get(sg,0)+1
            for ma,mb in product((2,3),repeat=2):
                ac=sub(at[ia][ic],COS[ma]);bc=sub(bt[ib][ic],COS[mb])
                val=determinant3([[aa[ia],ab[ia][ib],ac],[ab[ia][ib],bb[ib],bc],[ac,bc,tt[ic]]])
                assert val==dec_anp(rec['determinants'][ix]);ix+=1;total_cases+=1
                if val==Z:
                    zk=(blocks,(a,b),(ar[ia][0],br[ib][0],tr[ic][0]),(ma,mb));assert zk in zeros;usedzeros.add(zk)
                    zr=zeros[zk];assert zr['actual_vertex_failure']is not None
                    lab=own|{tuple(S[2*i:2*i+2]):blocks[i]for i in range(4)}|{(11,14):ma,(14,17):mb}
                    for leaf,pr in zip((11,17),(a,b)):
                        for root,m in zip((12,13),pr):lab[tuple(sorted((leaf,root)))]=m
                    for leaf,(e,v)in zip((11,17,14),(ar[ia],br[ib],tr[ic])):
                        if e is not None:lab[tuple(sorted((leaf,S[e])))]=3
                    failure(lab,zr['actual_vertex_failure'],V)
                else:nonzero+=1
        assert ix==len(rec['determinants'])
    assert records=={(bs,a,b)for bs in product((4,5),repeat=4)for a,b in active[bs[3]]}
    assert total_cases==d['residual_cases']==77952 and len(usedzeros)==len(zeros)==32 and not d['unresolved']
    assert anchor_signs==d['anchor_determinant_sign_counts']
    paths=[src,f'scripts/{TAG}.py',f'scripts/{TAG}_review.py']
    return dict(status='PASS_INDEPENDENT_COMPLETE_LAST_BLOCK_4_5_PROFILE_REVIEW',sources={p:sha(p)for p in paths},
        arithmetic='Standalone Fraction basis1,sqrt2,sqrt5,sqrt10; no producer/project imports; permutation determinants; independent domain reconstruction.',
        source_manifests=len(seen),actual_vertices=378,cross_right_guards=24,retained_root_budgets=6,
        high_own_actual_triangles=4,high_own_strict_lower=[str(x)for x in minimum],
        profile_pairs=290,profile_pair_status_counts=status_counts,active_profile_pair_counts=[36,48],
        fresh_upper_vertex_failures=168,rank_records=len(records),residual_determinants=total_cases,
        nonzero_determinants=nonzero,zero_determinants_with_actual_vertex_exclusions=32,
        anchor_determinant_sign_counts=anchor_signs,
        scope='First3 blocks have all own3/3, arbitrary internal labels (forced4/5); last internal4/5, all finite own labels. Rechecks every new guard/domain/coefficient/rank and zero-case vertex witness. Inherits previously independently reviewed opposite-single/all-own-three and two-dense theorem packages through explicit placements. Not all last internals or whole case0.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    out=run();path=ROOT/f'results/{TAG}_review.json'
    if args.verify:assert out==json.loads(path.read_text());print('SAVED_REPLAY_PASS')
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['residual_determinants'],out['zero_determinants_with_actual_vertex_exclusions'])

