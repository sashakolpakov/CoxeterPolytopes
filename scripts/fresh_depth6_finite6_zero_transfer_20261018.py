#!/usr/bin/env python3
"""Universal deep finite-six zero transfer with an accounted A1/A12 unit.

The report supplies arbitrary-support geometry. Exact profiles, exceptional
seals, protected triangles, literal interfaces and typed costs are checked.
"""
import argparse,json
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import sympy as S
import fresh_depth6_finite6_adaptive_pilot_20261018 as Z
import fresh_order15_d16_d17_ambient59_frontier_20261018 as C
import fresh_finite7_ideal5_triangle_payments_20261017 as TR

ROOT=Z.ROOT;P=Z.P;TAG='fresh_depth6_finite6_zero_transfer_20261018'
CYCLE='fresh_cycle11_half_pair_transfer_20261018'
STEMS=tuple(dict.fromkeys((*C.STEMS,C.TAG,Z.TAG,Z.OLD.TAG,Z.CAT,TR.TAG,CYCLE)))
enc=Z.enc


def count(t):return max(0,t[1]-{'A':18,'B':19,'D':20}[t[0]]) if t[0] in 'ABD' else 0
def actual_base(t):return C.actual_base(t)-16*count(t)
def upper_base(t):return C.upper_base(t)-16*count(t)
def actual_cross(a,b):return C.actual_cross(a,b)+(32 if {a,b}=={('A',1),('A',12)} else 0)


def graph(g,ids):return {i:[j for j in ids if i!=j and g[i,j]!=0] for i in ids}


def distances(g,a):
    d={a:0};q=deque([a])
    while q:
        i=q.popleft()
        for j in g[i]:
            if j not in d:d[j]=d[i]+1;q.append(j)
    return d


def seal(g,mark,tail):
    h=g.inv();b=h.extract(tail,tail);piv=[b[:i,:i].det() for i in range(1,len(tail)+1)]
    assert all(x>0 for x in piv)
    lower=S.simplify(h[mark,mark]-(h.extract([mark],tail)*b.inv()*h.extract(tail,[mark]))[0]);assert lower>2
    return dict(mark=mark,free_tail=tail,tail_gram=[[str(x) for x in r] for r in b.tolist()],
        positive_leading_determinants=list(map(str,piv)),projection_lower=str(lower))


def exceptions(pilot):
    ps=[{int(i):S.sympify(v) for i,v in x.items()} for x in pilot['profiles']['profiles']]
    def core(left,right):
        g=S.diag(P.path(19),S.diag(2,2))
        for i,v in left.items():g[i,19]=g[19,i]=-v
        for i,v in right.items():g[i,20]=g[20,i]=-v
        return g
    gs={
        'mixed':core({6:1},{11:1}),
        'mixed_cycle':core({6:1},{11:1,0:1}),
        'central':core({7:1},{11:1}),
        'central_cycle':core({7:1},{11:1,0:1})}
    rev=[18-i for i in range(19)]+[20,19];mapping=[]
    for row in pilot['schur']['remaining_branches']:
        assert row['c_patch']==row['d_patch']=='none'
        gap=None
        if row['t'] or row['u']:
            index=0 if row['t'] else 1
            value=sum((S.Rational(r['value']['rational'])+S.Rational(r['value']['sqrt2'])*S.sqrt(2))*(S.sqrt(2)-1)**r['powers'][index]
                      for r in row['coefficients'])
            gap=P.field(value);assert gap['sign']>0
            assert all(r['value']['sign']>=0 for r in row['coefficients'] if sum(r['powers']))
        c={**ps[row['c']],**({6:1} if row['t'] else {})}
        d={18-i:v for i,v in ps[row['d']].items()}
        if row['u']:d[12]=1
        g=core(c,d);matches=[(name,ref) for name,h in gs.items() for ref in (False,True) if g==(h.extract(rev,rev) if ref else h)]
        assert matches
        mapping.append(dict(branch=[row[k] for k in ('c','d','t','u')],native_gap=gap,canonical=matches[0]))
    assert len(mapping)==7
    proofs=[]
    for name,g in gs.items():
        h=g[:19,:19];v=g[:19,19:];sch=g[19:,19:]-v.T*h.inv()*v
        assert sch.det()<0
        assert C.Q.inertia([[F(x) for x in row] for row in g.tolist()])==[20,1,0]
        white=[i for i in range(19) if i not in(6,12)]
        for b in (6,12,19,20):assert C.Q.inertia([[F(g[i,j]) for j in white+[b]] for i in white+[b]])==[18,0,0]
        proof=dict(kind=name,gram=[[str(x) for x in r] for r in g.tolist()],inertia=[20,1,0],seals=[])
        if name=='mixed':
            proof['seals']=[seal(g,20,[0,18])];assert proof['seals'][0]['projection_lower']=='9/4'
            proof['paying_corner']=[6,20];proof['retained_unit']='marked leaf20 with fixed degree-two neighbor11'
        if name.endswith('cycle'):
            ids=list(range(12))+[20];cycle=g.extract(ids,ids)
            assert cycle*S.ones(13,1)==S.zeros(13,1) and cycle.rank()==12
            proof['whole_affine_cycle']=dict(nodes=ids,kernel=[1]*13,corner=[6,20])
            proof['seals']=[seal(g,19,[18])]
            assert proof['seals'][0]['projection_lower']==('78/7' if name=='mixed_cycle' else '78/31')
            corner=white+[19,20];adj=graph(g,corner)
            if name=='mixed_cycle':
                assert adj[19]==[]
                path=list(range(5,-1,-1))+[20]+list(range(11,6,-1))
                assert g.extract(path,path)==P.path(12) and path.index(20)==6
                proof['whole_finite_alternative']=dict(types=['A1','A12'],A1=[19],A12=path,marks=[19,20],A12_canonical_marks=[5,6],payment='1')
            else:
                path=list(range(5,-1,-1))+[20]+list(range(11,6,-1))+[19]
                assert g.extract(path,path)==P.path(13) and distances(adj,19)[20]==6 and len(adj[19])==1
                proof['whole_finite_alternative']=dict(type='A13',nodes=path,marks=[19,20],distance=6,payment='1')
        if name=='central':
            inv=g.inv();candidates=[];surv={19:[],20:[]}
            for mark in (19,20):
                for contacts in ({mark:1},{mark:S.sqrt(2)},{mark:1,0:1},{mark:1,18:1}):
                    q=S.expand(sum(v*w*inv[i,j] for i,v in contacts.items() for j,w in contacts.items()))
                    row=dict(mark=mark,contacts={i:str(v) for i,v in contacts.items()},projection=P.field(q),excluded=P.field(q-2)['sign']>0)
                    candidates.append(row)
                    if not row['excluded']:surv[mark].append(contacts)
            assert surv=={19:[{19:1},{19:1,18:1}],20:[{20:1},{20:1,0:1}]}
            pairs=[]
            for x,y in product(surv[19],surv[20]):
                inner=lambda a,b:sum(v*w*inv[i,j] for i,v in a.items() for j,w in b.items())
                b=S.Matrix([[2-inner(x,x),-inner(x,y)],[-inner(x,y),2-inner(y,y)]])
                assert b.det()<0
                pairs.append(dict(c_contacts=x,d_contacts=y,schur=[[str(z) for z in r] for r in b.tolist()],determinant=str(b.det())))
            assert [r['determinant'] for r in pairs]==['-13/12','-7/8','-7/8','-5/6']
            proof['extra_neighbor_domain']=candidates;proof['incompatible_neighbor_pairs']=pairs
            proof['retained_alternative']='At least one of19,20 is a whole marked simple leaf with fixed degree-two neighbor; finite distance6 or exceptional affine order<=9 unit.'
        proofs.append(proof)
    return dict(branch_reduction=mapping,proofs=proofs,
        extra_white_contacts_to_selected_path_only=[0,18],selected_finite_forest_forbids_contact_to_both_ends=True)


def protected_interfaces(docs):
    tri=TR.run();assert tri==docs[TR.TAG];protected=[]
    for row in tri['rows']:
        for c in row['corners']:
            if c['elementary_payment'] and c['length']==6:
                w=c['feature_witness'];assert w is not None and len(w['path'])-1<=5
                protected.append(dict(diagram=row['diagram'],marks=row['marks'],pair=c['pair'],feature_witness=w))
    assert len(protected)==54
    core=next(r for r in docs[Z.CAT]['diagrams'] if r['index']==25)
    graph0={i:[] for i in range(core['order'])}
    for i,j,m in core['edges']:graph0[i].append(j);graph0[j].append(i)
    special=[]
    for r in tri['exceptional_resolution']['rows']:
        if len(r['alternate_path'])-1!=6:continue
        omit=next(i for i in r['marks'] if i not in r['alternate_pair'])
        g={i:[j for j in js if j!=omit] for i,js in graph0.items() if i!=omit}
        dist=distances(g,r['fixed_simple_leaf']);depth=min(dist[m] for m in r['alternate_pair']);assert depth<=2
        special.append(dict(**r,marked_distance_to_protected_leaf=depth))
    assert len(special)==3
    cyc=docs[CYCLE]['profiles']['single_profiles'];long=[r for r in cyc if r['marked_distance']==6]
    assert [r['profile'] for r in long]==['end3_1','end4_1'] and all(r['fixed_leaf']==1 for r in long)
    assert 5-1==4<6
    return dict(fresh_high_core_triangle_triples=1025,protected_elementary_distance6=protected,
        protected_exceptional_distance6=special,cycle11_distance6=long,cycle11_old_mark_distance_to_fixed_leaf=4,
        low_triangle_core_bound='Two paths of length<=5 by the unchanged all-label Steiner theorem.',
        short5_scope='All current short5 ideal units except the independent cycle11 halves; their full proof uses only these protected finite6 outputs, shorter paths, or distinct retained halves.')


def counts():
    rows=[]
    for kind in 'ABD':
        for n in range(1 if kind=='A' else 2 if kind=='B' else 4,81):
            edges=[(i,i+1,3) for i in range(n-1)] if kind!='D' else [(i,i+1,3) for i in range(n-3)]+[(n-3,n-2,3),(n-3,n-1,3)]
            if kind=='B':edges[-1]=(n-2,n-1,4)
            adj={i:[] for i in range(n)}
            for i,j,m in edges:adj[i].append(j);adj[j].append(i)
            features={i for i in adj if len(adj[i])!=2}|{i for a,b,m in edges if m!=3 for i in(a,b)}
            ds=[distances(adj,i) for i in range(n)];depth=[min(ds[i][f] for f in features) for i in range(n)]
            pairs=[(i,j) for i,j in combinations(range(n),2) if ds[i][j]==6 and min(depth[i],depth[j])>=6]
            assert len(pairs)==count((kind,n));rows.append(dict(kind=kind,rank=n,deleted_pairs=pairs))
    return dict(formulas=dict(A='max(0,n-18)',B='max(0,n-19)',D='max(0,n-20)'),rows=rows)


def coupon(g,blocks,marks):
    a,b=[next(b for b in blocks if m in b['nodes']) for m in marks]
    if a is b:return []
    out=[]
    for sm,la,m,n in ((a,b,*marks),(b,a,marks[1],marks[0])):
        if sm['type']=='A1' and la['type']=='A12':
            images=C.T.mark_pairs(g,la['nodes'],[n],12,'A')
            if images==[(5,),(6,)]:out.append(dict(small=sm,large=la,marks=[m,n],images=images,weight='1'))
    assert len(out)<=1;return out


def literal(docs):
    old={C.R.key(r['source'],r['case']):r for r in docs[C.TAG]['combined_literal_interface']['rows']};rows=[];changed=[]
    for stem,index,g,white,cs in C.T.literal_inputs(docs):
        before=old[C.R.key(stem,index)];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==before['corners'][j]['marks'];ids=sorted(set(white)|set(marks));sig=C.Q.inertia(C.A.sub(g,ids));debit=F();cp=[]
            if sig==[len(ids),0,0]:
                bs=C.D.finite_blocks(g,ids);cp=coupon(g,bs,marks)
                b=next((b for b in bs if set(marks)<=set(b['nodes'])),None)
                if b:
                    ns=b['nodes'];adj={i:[k for k in ns if i!=k and g[i][k]!=0] for i in ns}
                    feature={i for i in ns if len(adj[i])!=2 or any(g[i][k]**2!=1 for k in adj[i])}
                    dd=[distances(adj,m) for m in marks]
                    if dd[0][marks[1]]==6 and min(min(x[i] for i in feature) for x in dd)>=6:debit=F(1)
            inc=F(len(cp));pay=F(before['corners'][j]['retained_payment_lower_bound'])-debit+inc;assert pay>=0
            corners.append(dict(marks=marks,debit=debit,coupons=cp,increment=inc,retained_payment_lower_bound=pay))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(c['debit'] or c['increment'] for c in corners):changed.append(row)
    assert len(rows)==216
    return dict(configurations=216,changed=changed,rows=rows)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    fresh=Z.run();assert fresh==docs[Z.TAG]
    assert fresh['profiles']['accepted_count']==19 and fresh['schur']['branch_count']==829 and fresh['schur']['strict_branches']==822
    rep=C.P.Z.B.B.representative;types=[t for t in C.P.Z.B.B.types() if t[1]<=59]
    for a in types:
        assert upper_base(rep(a))>=actual_base(a)>=0
        for b in types:assert actual_cross(a,b)==actual_cross(b,a)<=C.upper_cross(rep(a),rep(b))
    assert C.actual_cross(('A',1),('A',12))==0 and actual_cross(('A',1),('A',12))==32 and C.upper_cross(('A',1),('D',12))==40
    paths={Path(__file__),Path(Z.__file__),Path(P.__file__),Path(C.__file__),Path(TR.__file__)}
    paths.update(ROOT/f'{f}/{s}.{e}' for s in STEMS for f,e in (('scripts','py'),('results','json')))
    paths.update(ROOT/f'results/{s}.md' for s in STEMS if s!=Z.TAG)
    return enc(dict(status='PASS_DEPTH6_FINITE6_ZERO_LOCAL_TRANSFER',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        feature_depth=6,finite_distance=6,new_weight='0',profiles=fresh['profiles'],minimal_core_screen=fresh['minimal_core_screen'],
        adaptive_schur=fresh['schur'],exceptional_payment_proofs=exceptions(fresh),protected_interfaces=protected_interfaces(docs),
        deleted_pair_counts=counts(),current_literal_interface=literal(docs),simultaneous_actual_type_pairs=len(types)**2,
        new_coupon=dict(small='whole_A1',large='whole_A12',large_marks=[5,6],weight='1',component_pair_increment='2',
            new_actual_cost='2',unchanged_representative_allowance='5/2'),
        all_ideal_weights_unchanged=True,all_representative_cross_costs_unchanged=True,triangle_minimum=2,quadrilateral_minimum=1,
        finite_forest_count_performed=False,new_dimension_bound=None,
        scope='Universal arbitrary-support local rule; complete retained profiles and adaptive minimal-core reductions, seven unbounded/native branches reduced to four exact seals, independent short5 and protected triangle composition. All added actual payments included. No new finite forest maximum or ambient bound.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();p=ROOT/f'results/{TAG}.json';out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'literal_changed',len(out['current_literal_interface']['changed']),flush=True)
