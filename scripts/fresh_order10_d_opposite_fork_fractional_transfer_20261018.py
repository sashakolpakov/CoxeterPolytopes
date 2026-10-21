#!/usr/bin/env python3
"""Order10 opposite-fork fractional transfer with explicit finite charges.

Exact local identities and complete bounded consequences of all-rank
inverse inequalities. No finite forest DP or cusp convolution is run.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import sympy as s
import fresh_order10_cycle_half_transfer_20261018 as A10
import fresh_order10_b_terminal_half_transfer_20261018 as B10
import fresh_order11_bcd_half_pair_transfer_20261018 as K
import fresh_depth6_finite8_ambient58_frontier_20261018 as BASE

ROOT=BASE.ROOT;TAG='fresh_order10_d_opposite_fork_fractional_transfer_20261018'
STEMS=(A10.TAG,B10.TAG,B10.CZERO,K.TAG,BASE.TAG,BASE.LOCAL,
 'fresh_spectator_energy_below_one_20261017','fresh_spectator_five_quarters_transfer_20261017',
 'fresh_a17_five_pair_alternate_transfer_20261018','fresh_finite_d15_a4_a11_half_zero_transfer_20261018',
 'fresh_finite_d16_alternate_half_zero_transfer_20261018','fresh_finite_d17_leaf_a16_half_transfer_20261018',
 'fresh_a15_d15_dimension58_20261018')
DELTA=F(5,32);WEIGHT=1-DELTA;N,R,A,B=10,9,6,8
CHANGES=((6,8),(6,9),(7,8),(7,9))


def affine():
    g=K.gram(10,[(i,i+1,1)for i in range(5)]+[(0,6,1),(0,7,1),(5,8,1),(5,9,1)])
    nu=s.Matrix([1]*6+[s.Rational(1,2)]*4)
    assert g*nu==s.zeros(10,1)and g.rank()==9
    return g,nu


def profiles(mark):
    mate,center=(7,0)if mark==A else(9,5);t=s.Symbol('t',nonnegative=True)
    return [dict(name=name,contacts={mark:t,**extra},moment=l,c0=c0,height=h,t=t)
        for name,extra,l,c0,h in zip(('isolated','end3','end4','next3'),
        ({},{mate:1},{mate:s.sqrt(2)},{center:1}),(0,1,s.sqrt(2),2),(2,1,0,0),
        (t/2,(t+1)/2,(t+s.sqrt(2))/2,1+t/2))]


def physical():
    g,nu=affine();W=[i for i in range(10)if i not in(A,B)];J=g.extract(W,W).inv()
    vs=[-nu[i]*g.extract([i],W).T for i in(A,B)]
    sc=s.diag(nu[A],nu[B]);assert sc*g.extract([A,B],[A,B])*sc-s.Matrix.hstack(*vs).T*J*s.Matrix.hstack(*vs)==s.Matrix([[1,-1],[-1,1]])/9
    rows=[];invrows=[];immediate=[];cross=[];u=s.Symbol('u',nonnegative=True)
    for mark,other,far in((A,B,vs[1]),(B,A,vs[0])):
        ids=[i for i in range(10)if i!=mark];deleted=g.extract(ids,ids);K.pd(deleted);inv=deleted.inv()
        assert all(z.is_nonnegative for z in inv)
        allow=[ids[i]for i in range(9)if inv[i,i]<=2];expected=[0,7]if mark==A else[5,9];assert allow==expected
        ii=[ids.index(i)for i in allow];two=sum(inv[i,j]for i in ii for j in ii);assert two>2
        assert sorted(inv[i,i]for i in ii)==[1,2]and inv[ids.index(other),ids.index(other)]==s.Rational(9,4)
        invrows.append(dict(mark=mark,nodes=ids,diagonal=list(map(str,inv.diagonal())),allowed=allow,two_contact_minimum=str(two)))
        for p in profiles(mark):
            x=s.Matrix([p['contacts'].get(i,0)for i in W]);h,l,c0,t=[p[k]for k in('height','moment','c0','t')]
            assert s.simplify(sum(nu[i]*q for i,q in p['contacts'].items())-h)==0
            assert s.simplify((far.T*J*x)[0]-l/9)==0 and s.simplify(2-(x.T*J*x)[0]-c0-l*l/9)==0
            cost=s.factor(2*l/h+c0/h**2);minimum=1 if p['name']=='isolated'else 0;maximum=s.simplify(cost.subs(t,minimum))
            num,den=s.fraction(s.factor((maximum-cost).subs(t,u+minimum)))
            assert all(z.is_nonnegative for z in s.Poly(num,u).all_coeffs()+s.Poly(den,u).all_coeffs())
            row=dict(mark=mark,profile=p['name'],height=str(h),moment=str(l),c0=str(c0),cost=str(cost),maximum=str(maximum))
            if p['name']in('isolated','end3'):
                cutoff=s.sqrt(2)if p['name']=='isolated'else 1
                assert s.simplify(cost.subs(t,cutoff))<=4
                row['nonextreme_maximum']=str(s.simplify(cost.subs(t,cutoff)))
            rows.append(row)
            if p['name']in('end4','next3'):
                gg=2*s.eye(10);gg[:9,:9]=deleted
                for i,q in p['contacts'].items():
                    if i!=mark:gg[ids.index(i),9]=gg[9,ids.index(i)]=-q
                null=gg.nullspace();assert len(null)==1;v=null[0]if null[0][0].is_positive else-null[0];assert all(z.is_positive for z in v)
                kind=B10.type_of_affine(gg,(9,ids.index(other)));assert kind==('B'if p['name']=='end4'else'D')
                pay=F(1,2)if kind=='B'else WEIGHT;assert WEIGHT+pay>=1
                immediate.append(dict(mark=mark,profile=p['name'],target=kind,kernel=list(map(str,v)),distinct_payment=str(pay),total=str(WEIGHT+pay)))
    for p,q in product(profiles(A),profiles(B)):
        x=s.Matrix([p['contacts'].get(i,0)for i in W]);y=s.Matrix([q['contacts'].get(i,0)for i in W]);z=s.simplify((x.T*J*y)[0])
        assert s.simplify(z-p['moment']*q['moment']/9)==0
        cross.append(dict(left=p['name'],right=q['name'],value=str(z)))
    h,k,l,m,c0,d0,p,q,r=s.symbols('h k l m c0 d0 p q r',nonzero=True);alpha=s.Rational(1,9)
    core=s.Matrix([[alpha,-alpha,-h+l*alpha,-m*alpha],[-alpha,alpha,-l*alpha,-k+m*alpha],[-h+l*alpha,-l*alpha,c0+l*l*alpha-p,-l*m*alpha-q],[-m*alpha,-k+m*alpha,-l*m*alpha-q,d0+m*m*alpha-r]])
    formula=1-(2*l/h+c0/h**2+2*m/k+d0/k**2)/9+(p/h**2-2*q/(h*k)+r/k**2)/9
    assert s.factor(core.det()/(h*h*k*k)-formula)==0
    assert 1-s.Rational(8,9)==s.Rational(1,9)and s.Rational(16-9,4)==s.Rational(7,4)
    gap=s.factor(2+1/(2*k*k)-(2/(k*k)-1));assert s.simplify(gap-3+3/(2*k*k))==0
    return dict(gram=[list(map(str,row))for row in g.tolist()],positive_kernel=list(map(str,nu)),inverses=invrows,profiles=rows,cross_identities=cross,
        immediate=immediate,determinant=str(formula),nonextreme_gap='1/9',one_extreme_energy_maximum='3',isolated_two_coordinate_gap=str(gap),two_extreme_energy='7/4')


def spectrum():
    cases=[];bounds=[]
    assert F(30,17)>F(7,4)
    for d in range(4,15):
        bound=F(4*d*d,4*d-7);bounds.append(dict(distance=d,N_at_least=d+2,N_at_most=str(bound)))
        for total in range(d+2,bound.numerator//bound.denominator+1):
            for i in range(1,total-d):
                j=i+d;a=F(i*(total-i),total);b=F(j*(total-j),total);c=F(i*(total-j),total)
                if max(a,b)>=2:continue
                assert a+b-2*c==F(d*(total-d),total)<=F(7,4)
                cases.append(dict(kind=f'A{total-1}',marks=[i-1,j-1],distance=d,block=[a,b,c]))
    for rank in(6,7):
        cases.append(dict(kind=f'D{rank}',marks=[0,rank-2],distance=rank-2,block=[F(1),F(rank,4),F(1,2)]))
    cases.append(dict(kind='E6',marks='eligible endpoints',distance=4,block=[F(4,3),F(4,3),F(2,3)]))
    exact=[]
    for row in cases:
        a,b,c=map(s.Rational,row.pop('block'));energy=a+b-2*c
        assert a-c>0 and b-c>0 and energy+s.Rational(1,2)>s.Rational(7,4)
        gradients=(a*s.sqrt(2)-c,b-c*s.sqrt(2),a-c*s.sqrt(2),b*s.sqrt(2)-c)
        assert all(s.simplify(z).is_positive for z in gradients)
        nonunit=(2*a+b-2*s.sqrt(2)*c,a+2*b-2*s.sqrt(2)*c)
        assert all(s.simplify(z-s.Rational(7,4)).is_positive for z in nonunit)
        row.update(inverse_block=list(map(str,(a,b,c))),unit_energy=str(energy),nonunit_minima=list(map(str,nonunit)))
        if energy==s.Rational(7,4):exact.append(row['kind'])
    assert exact==['A15','D7']
    return dict(A_far_bounds=bounds,complete_far_inverse_blocks=cases,exact_total_components=exact,
        minimum_other_active_component='1/2',no_mixed_one_side_maximum='7/8',no_mixed_whole_A_rank_maximum=7,
        all_rank_reason='A distance>=15 has lower bound30/17>7/4; explicit remaining distance/N ranges exhaust A. B has no far allowable pair, D onlyD6/D7, E onlyE6; F/H/I2 diameter<=3. Native nonunit minima all exceed7/4.')


def finite_coupons():
    f,nu=affine();rows=[]
    for attachment,r in product((A,7),range(8)):
        c=10;g=2*s.eye(11+r);g[:10,:10]=f;g[c,attachment]=g[attachment,c]=-1
        if r:
            for i,j in[(c,11)]+[(i,i+1)for i in range(11,10+r)]:g[i,j]=g[j,i]=-1
        ids=[i for i in range(g.rows)if i!=A];K.pd(g.extract(ids,ids));parts=K.components(g,ids)
        cp=next(z for z in parts if c in z);bp=next(z for z in parts if B in z);ad=K.graph(g)
        if attachment==A:
            assert len(cp)==r+1 and len(bp)==9 and cp!=bp and len(parts)==2
            assert len(ad[c]&set(cp))==(0 if r==0 else 1)
            assert len(ad[B]&set(bp))==1 and len(ad[next(iter(ad[B]&set(bp)))]&set(bp))==3
            pay=DELTA;role=f'wholeA{r+1} endpoint/wholeD9 shortleaf'
        else:
            assert cp==bp and len(cp)==10+r
            adj=K.graph(g.extract(cp,cp));dist=K.distances(adj,cp.index(c))[cp.index(B)]
            assert dist==8 and len(adj[cp.index(B)])==1
            pay=F(1);role=f'wholeD{10+r} marked shortleaf distance8'
        assert WEIGHT+pay>=1
        rows.append(dict(spectator_rank=r,old_attachment=attachment,whole_components=parts,marks=[c,B],role=role,payment=str(pay),total=str(WEIGHT+pay)))
    return rows


def mixed_corners():
    f,nu=affine();a,b,c,d=A,B,10,11;rows=[]
    for kind,r,edges,contacts in(('A15',15,[(i,i+1)for i in range(14)],(0,14)),('D7',7,[(i,i+1)for i in range(4)]+[(4,5),(4,6)],(0,5))):
        for ori in range(1 if kind=='A15'else 2):
            for la,rb in product((6,7),(8,9)):
                g=2*s.eye(12+r);g[:10,:10]=f
                for i,j in[(c,la),(d,rb)]+[(12+i,12+j)for i,j in edges]:g[i,j]=g[j,i]=-1
                cp=contacts if not ori else contacts[::-1]
                for v,z in zip((c,d),cp):g[v,12+z]=g[12+z,v]=-1
                W=[i for i in range(g.rows)if i not in(a,b,c,d)];edgesigs=[K.inertia(g.extract(W+[v],W+[v]))for v in(a,b,c,d)]
                corners=[]
                for marks in((a,b),(b,c),(c,d),(a,d)):
                    ids=W+list(marks);sig=K.inertia(g.extract(ids,ids));parts=K.components(g,ids);whole=next(z for z in parts if marks[0]in z)
                    row=dict(marks=list(marks),inertia=sig,components=parts)
                    if marks[1]in whole:
                        sub=g.extract(whole,whole);ad=K.graph(sub);row.update(whole_order=len(whole),marked_distance=K.distances(ad,whole.index(marks[0]))[whole.index(marks[1])],marked_degrees=[len(ad[whole.index(v)])for v in marks])
                    corners.append(row)
                viable=all(sig[1:]==[0,0]for sig in edgesigs)and all(z['inertia'][1]==0 for z in corners)
                payment=None
                if viable:
                    if kind=='A15'and(la,rb)==(A,B):
                        out=corners[2];assert out['whole_order']==17 and out['marked_distance']==16 and out['marked_degrees']==[1,1]
                        payment=dict(corner=2,role='wholeA17 endpoints restored5/32',weight=str(DELTA))
                    elif kind=='A15':
                        index=1 if la==7 else 3;out=corners[index]
                        assert out['whole_order']==25 and out['marked_distance']==8 and sorted(out['marked_degrees'])==[1,2]and out['inertia'][1:]==[0,0]
                        payment=dict(corner=index,role='wholeD25 shortleaf distance8',weight='1')
                    else:
                        assert(la,rb)==(A,B);out=corners[2]
                        assert out['whole_order']==9 and out['marked_distance']==7 and out['marked_degrees']==[1,1]and out['inertia'][1:]==[0,1]
                        sub=g.extract(next(z for z in out['components']if c in z),next(z for z in out['components']if c in z));null=sub.nullspace();assert len(null)==1
                        v=null[0]if null[0][0].is_positive else-null[0];assert all(z.is_positive for z in v)
                        payment=dict(corner=2,role='whole affineE8 order9 unit',weight='1',positive_kernel=list(map(str,v)))
                    assert WEIGHT+F(payment['weight'])>=1
                rows.append(dict(spectator=kind,orientation=ori,old_attachments=[la,rb],edge_inertias=edgesigs,corners=corners,viable=viable,payment=payment,gram=[list(map(str,row))for row in g.tolist()]))
    assert len(rows)==12 and sum(z['viable']for z in rows)==6
    return rows


def actual_base(t):return 2*BASE.actual_base(t)+(5 if t==('A',17)else 0)
def upper_base(t):return 2*BASE.upper_base(t)+(5 if t==('A',17)else 0)
def actual_cross(a,b):
    extra=next((10 if r==1 else 20 for r in range(1,9)if frozenset((a,b))==frozenset((('A',r),('D',9)))),0)
    return 2*A10.actual_cross(a,b)+extra

def upper_cross(a,b):
    extra=20 if frozenset((a,b))in{frozenset((t,('B',9)))for t in(('A',5),('E',6),('A',7),('A',8))}else 0
    return 2*B10.upper_cross(a,b)+extra


def run():
    docs={tag:json.loads((ROOT/f'results/{tag}.json').read_text())for tag in STEMS}
    for tag,doc in docs.items():
        assert doc['status'].startswith('PASS')
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    ph=physical();sp=spectrum();fc=finite_coupons();mc=mixed_corners()
    f,nu=affine();ad=K.graph(f);weights={p:F(1)for p in combinations(range(10),2)}
    for pair in CHANGES:assert K.distances(ad,pair[0])[pair[1]]==7;weights[pair]=WEIGHT
    assert sum(weights.values())==F(355,8)
    old_a17={tuple(r['pair']):F(r['weight'])for r in docs['fresh_a17_five_pair_alternate_transfer_20261018']['weighted_pairs']}
    assert old_a17[0,16]==0 and sum(old_a17.values())==F(BASE.actual_base(('A',17)),16)
    new_a17=old_a17|{(0,16):DELTA};assert sum(new_a17.values())==F(actual_base(('A',17)),32)
    retained=[]
    for rank,tag,field in ((15,'fresh_finite_d15_a4_a11_half_zero_transfer_20261018','weighted_pairs'),(16,'fresh_finite_d16_alternate_half_zero_transfer_20261018','D16_weighted_pairs'),(17,'fresh_finite_d17_leaf_a16_half_transfer_20261018','weighted_pairs')):
        gg=K.gram(rank,[(i,i+1,1)for i in range(rank-3)]+[(rank-3,rank-2,1),(rank-3,rank-1,1)]);adj=K.graph(gg)
        for row in docs[tag][field]:
            i,j=row['pair']
            if K.distances(adj,i)[j]<=8:assert F(row['weight'])==1
        retained.append(dict(rank=rank,all_distance_at_most8_unit=True))
    allowed_a={0,7};allowed_b={5,9};assert not allowed_a&allowed_b and A not in allowed_b and B not in allowed_a
    engine=B10.C.C.C.P.Z.B.B;types=engine.types();rep=engine.representative
    for a in types:
        assert upper_base(rep(a))>=actual_base(a)>=0
        for b in types:assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)==actual_cross(b,a)>=0
    paths={Path(__file__),Path(A10.__file__),Path(B10.__file__),Path(K.__file__),Path(BASE.__file__)}
    paths.update(ROOT/f'{folder}/{tag}.{ext}'for tag in STEMS for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    paths.add(ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md')
    return dict(status='PASS_ORDER10_D_OPPOSITE_FORK_FRACTIONAL_TRANSFER',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        delta=str(DELTA),physical=ph,spectator_spectrum=sp,whole_finite_coupons=fc,mixed_whole_corners=mc,
        weighted_pairs=[dict(pair=p,weight=str(w))for p,w in sorted(weights.items())],old_cost='45',new_cost='355/8',
        actual_component_increments={f'A{r}/D9':str(2*DELTA if r==1 else 4*DELTA)for r in range(1,9)},
        representative_increments_beyond_B10={f'{t}/B9':'5/8'for t in('A5','E6','A7','A8')},
        restored_A17_pair=[0,16],A17_increment=str(DELTA),A17_new_cost=str(F(actual_base(('A',17)),32)),
        A17_weighted_pairs=[dict(pair=p,weight=str(w))for p,w in sorted(new_a17.items())],
        retained_finite_short8_guards=retained,ordered_type_pairs=len(types)**2,scale=32,triangle=dict(allowed_a=sorted(allowed_a),allowed_b=sorted(allowed_b),affected_triangle_impossible=True),
        current_order10_family_costs=dict(A='85/2',B='44',C='44',D='355/8'),new_order10_envelope='71/16',
        triangle_minimum=2,quadrilateral_minimum=1,new_forest_count=False,new_dimension_bound=None,
        scope='Universal local fractional theorem with full arbitrary-support spectrum, explicit ordinary coupons and A17 restoration. No new forest budget or concentration theorem is asserted.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();path=ROOT/f'results/{TAG}.json'
    if not args.verify:assert not path.exists()
    out=json.loads(json.dumps(run()))
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_cost'],out['A17_new_cost'],flush=True)
