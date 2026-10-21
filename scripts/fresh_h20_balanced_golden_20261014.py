#!/usr/bin/env python3
"""Exact coupled radial caps in a golden A6 x A6 x A7 cusp model."""
import argparse,hashlib,json
from collections import Counter
from fractions import Fraction as F
from itertools import product,combinations
from pathlib import Path
from fresh_h20_balanced_20261014 import ROOT,exact_data

TAG='fresh_h20_balanced_golden_20261014'
Z=(F(0),F(0));O=(F(1),F(0));P=(F(0),F(1))
def qa(a,b):return a[0]+b[0],a[1]+b[1]
def qs(a,b):return a[0]-b[0],a[1]-b[1]
def qm(a,b):return a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1]
def qt(a,t):return a[0]*t,a[1]*t
def q(x):return F(x),F(0)
def qsg(a):
    x,y=a[0]+a[1]/2,a[1]/2
    if not x:return (y>0)-(y<0)
    if not y:return (x>0)-(x<0)
    if (x>0)==(y>0):return (x>0)-(x<0)
    d=x*x-5*y*y
    if not d:return 0
    return ((x>0)-(x<0))if d>0 else((y>0)-(y<0))
def qi(a):
    b=(a[0]+a[1],-a[1]);d=qm(a,b);assert d[1]==0 and d[0]
    return qt(b,1/d[0])
N=(F(196),F(84));D=qt(qs(N,O),12)
ZERO=(Z,Z);ONE=(O,Z)
def c(x):return q(x),Z
def lift(x):return x,Z
def add(a,b):return qa(a[0],b[0]),qa(a[1],b[1])
def sub(a,b):return qs(a[0],b[0]),qs(a[1],b[1])
def mul(a,b):return qa(qm(a[0],b[0]),qm(D,qm(a[1],b[1]))),qa(qm(a[0],b[1]),qm(a[1],b[0]))
def scale(a,t):return qt(a[0],t),qt(a[1],t)
def sign(a):
    s,t=qsg(a[0]),qsg(a[1])
    if not s:return t
    if not t or s==t:return s
    u=qsg(qs(qm(a[0],a[0]),qm(D,qm(a[1],a[1]))))
    return s*u
def square(a):return mul(a,a)
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return [enc(y)for y in x]
    if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
    return x
def ratdot(a,b):
    return sum(x*y for x,y in zip(a,b,strict=True))-sum((a[i]*b[i+1]+a[i+1]*b[i])/2 for lo,hi in[(0,6),(6,12),(12,19)]for i in range(lo,hi-1))
def dist(a,b):
    x=abs(a[0]-b[0]);y=abs(a[1]-b[1]);z=abs(a[2]-b[2])
    return qa(q(14*x*(7-x)+14*y*(7-y)),qt(qm(P,P),16*z*(8-z)))
def run():
    base=exact_data();fv=base['factor_vertices'];assert qsg(D)>0
    basis=[[F(i==j)for j in range(19)]for i in range(19)]
    assert base['spatial_Gram']==[[ratdot(a,b)for b in basis]for a in basis]
    assert all(ratdot(v,v)==1 for v in base['seed_vectors'])
    # Reconstruct all three factor metric formulas from the rational source.
    for f,h in enumerate([7,7,8]):
        for i,j in product(range(h),repeat=2):
            v=[x-y for x,y in zip(fv[f][i],fv[f][j],strict=True)]
            a=abs(i-j);assert ratdot(v,v)==2*h*a*(h-a)
        assert all(ratdot(v,v)==[56,56,84][f]for v in fv[f])
    assert N==qa(q(112),qt(qm(P,P),84))
    den=qi(qs(N,q(4)));k=(qm(qa(N,q(2)),den),qt(den,-1))
    assert sign(k)>0 and sign(sub(ONE,k))>0
    assert sub(add(mul(lift(qs(N,q(4))),square(k)),lift(qs(N,q(4)))),scale(mul(lift(qa(N,q(2))),k),2))==ZERO
    powers=[ONE]
    for j in range(8):powers.append(mul(powers[-1],k))
    vertices=list(product(range(7),range(7),range(8)))
    distance_counts=Counter(dist(a,b)for a,b in combinations(vertices,2))
    assert sum(distance_counts.values())==392*391//2
    for dd in distance_counts:assert dd==q(84)or qsg(qs(dd,q(112)))>0
    ray_code=[(a%7,a//7,a)for a in range(8)]
    assert len(set(ray_code))==8
    ray_dists=[dist(a,b)for a,b in combinations(ray_code,2)]
    min_ray=ray_dists[0]
    for dd in ray_dists:
        if qsg(qs(dd,min_ray))<0:min_ray=dd
    def gate_dotted(num,r2,s2):
        assert sign(num)>0 and sign(sub(square(num),scale(mul(r2,s2),4)))>0
    def cap(j,om):
        de=powers[j];ga=sub(ONE,de);rr=scale(de,2)
        return dict(level=j,omitted=om,delta=de,gamma=ga,radius=rr,radius_square=square(rr))
    configurations=[]
    for name,rr2,specs in[
        ('two_step_retaining_large_endpoint_layer',c(28),[((0,0,0),j)for j in[1,2]]),
        ('eight_rays_linked_to_endpoint_layer',mul(lift(qs(N,q(4))),square(powers[8])),[(om,j)for om in ray_code for j in range(1,9)]),
        ('full_second_layer_with_eight_inner_bridges',c(28),[(om,1)for om in ray_code]+[(om,2)for om in vertices])]:
        if name.startswith('eight'):assert sign(sub(c(21),rr2))>0
        counts=Counter();newcounts=Counter();caps=[cap(j,om)for om,j in specs]
        old_gate_cache={};new_gate_cache={}
        # Literal original vertical-seed Gram, and its central contacts.
        for a,b in combinations(range(22),2):
            v=ratdot(base['seed_vectors'][a],base['seed_vectors'][b]);assert v in[0,F(-1,2)]
            counts['right'if not v else'pi3']+=1
        counts['pi3']+=14;counts['pi5']+=8
        # Endpoint layer: all omitted ratios are dotted; finite contacts right.
        for sf,A in [(O,7),(O,7),(P,8)]:assert sign(sub(lift(qt(qm(sf,sf),A*A)),rr2))>0
        counts['right']+=392*19;counts['dotted']+=392*3
        cp=sub(sub(lift(N),rr2),c(4));gate_dotted(cp,rr2,c(4));counts['dotted']+=392
        for dd,multiplicity in distance_counts.items():
            num=sub(lift(dd),scale(rr2,2))
            if num==rr2:counts['pi3']+=multiplicity
            else:gate_dotted(num,rr2,rr2);counts['dotted']+=multiplicity
        # Exact new seed and spherical pair gates.
        chain_contacts=[]
        for i,row in enumerate(caps):
            ga,de,R,r2=row['gamma'],row['delta'],row['radius'],row['radius_square']
            assert sign(ga)>0 and sign(sub(ONE,ga))>0
            for sf,h in [(O,7),(O,7),(P,8)]:
                assert mul(lift(sf),de)==mul(R,lift(qt(sf,F(1,2))))
                omitted=mul(lift(sf),add(de,scale(ga,h)))
                assert sign(sub(omitted,R))>0
            newcounts['pi3']+=12;newcounts['pi5']+=7;newcounts['dotted']+=3
            center2=mul(lift(N),square(ga));num=sub(sub(center2,r2),c(4))
            if row['level']==1:
                assert num==scale(R,2);newcounts['pi3']+=1
                chain_contacts.append(['central',i,'pi3'])
            else:gate_dotted(num,r2,c(4));newcounts['dotted']+=1
            for om in vertices:
                dd=dist(row['omitted'],om);key=(row['level'],dd)
                if key not in old_gate_cache:
                    distance=add(mul(lift(N),square(de)),mul(ga,lift(dd)))
                    num=sub(sub(distance,r2),rr2)
                    if row['level']==8 and dd==Z:
                        assert num==ZERO;typ='right'
                    else:gate_dotted(num,r2,rr2);typ='dotted'
                    # Existing finite endpoint has height²=rr2 and stays strict.
                    assert sign(add(sub(distance,r2),rr2))>0
                    old_gate_cache[key]=typ
                typ=old_gate_cache[key];newcounts[typ]+=1
                if typ=='right':
                    chain_contacts.append([i,'endpoint_'+str(vertices.index(om)),'right'])
            for j,old in enumerate(caps[:i]):
                dd=dist(row['omitted'],old['omitted']);key=(row['level'],old['level'],dd)
                if key not in new_gate_cache:
                    distance=add(mul(lift(N),square(sub(ga,old['gamma']))),mul(mul(ga,old['gamma']),lift(dd)))
                    num=sub(sub(distance,r2),old['radius_square'])
                    if dd==Z and abs(row['level']-old['level'])==1:
                        assert num==mul(R,old['radius']);typ='pi3'
                    else:gate_dotted(num,r2,old['radius_square']);typ='dotted'
                    new_gate_cache[key]=typ
                typ=new_gate_cache[key];newcounts[typ]+=1
                if typ=='pi3':chain_contacts.append([j,i,'pi3'])
            # H=(0,200,1), H²=-400, is strict for every new sphere.
            assert sign(add(sub(center2,r2),c(400)))>0
        allcounts=counts+newcounts;wallcount=415+len(caps)
        assert sum(allcounts.values())==wallcount*(wallcount-1)//2
        # Actual old endpoints: seedPD19 + perpendicular endpoint normal.
        assert all(F(j+1,2*j)>0 for n in[6,6,7]for j in range(1,n+1))
        assert sign(rr2)>0 and sign(add(sub(lift(N),rr2),c(400)))>0
        # The selected radial point was strict on the central+endpoint core.
        lam=F(7,20);x2=qt(N,lam*lam)
        target_power=[]
        for row in caps:
            distance=add(mul(lift(N),square(sub(c(lam),row['gamma']))),
                scale(mul(row['gamma'],lift(dist((0,0,0),row['omitted']))),lam))
            target_power.append(sub(distance,row['radius_square']))
        cutters=[i for i,v in enumerate(target_power)if sign(v)<0];assert cutters
        assert sign(sub(lift(x2),c(4)))>0
        for om in vertices:
            basepower=sub(add(lift(qt(N,(1-lam)**2)),lift(qt(dist((0,0,0),om),lam))),rr2)
            assert sign(sub(basepower,c(1)))>0
        # A different product direction gives a strict final escape.
        same_ray_power=sub(mul(lift(N),square(sub(c(lam),sub(ONE,powers[2])))),scale(square(powers[2]),4))
        assert sign(same_ray_power)<0
        assert qsg(qs(qt(N,lam*lam),q(36)))>0
        assert sign(sub(mul(lift(N),square(sub(ONE,powers[2]))),square(sub(c(6),scale(powers[2],2)))))>0
        target_coverage=len({tuple(row['omitted'])for row in caps if row['level']==2})
        if name.startswith('full'):
            escape_om=(3,3,0);escape_lam=F(1,8)
            assert {tuple(row['omitted'])for row in caps if row['level']==2}==set(vertices)
            assert target_coverage==392
        else:escape_om=(0,0,1);escape_lam=lam
        escape_x2=qt(N,escape_lam**2);powers_at_escape=[sub(lift(escape_x2),c(4))]
        for om in vertices:
            powers_at_escape.append(sub(add(lift(qt(N,(1-escape_lam)**2)),lift(qt(dist(escape_om,om),escape_lam))),rr2))
        for row in caps:
            distance=add(mul(lift(N),square(sub(c(escape_lam),row['gamma']))),
                scale(mul(row['gamma'],lift(dist(escape_om,row['omitted']))),escape_lam))
            powers_at_escape.append(sub(distance,row['radius_square']))
        minimum=powers_at_escape[0]
        for v in powers_at_escape:
            assert sign(sub(v,c(1)))>0
            if sign(sub(v,minimum))<0:minimum=v
        escapeD=qt(qs(escape_x2,O),F(1,2));assert qsg(qs(q(-200),escapeD))<0
        # Every new sphere makes ten explicit affineA2 triples with adjacent
        # retained A6 seed nodes. This is a noncompact construction.
        for row in caps:
            old_A6=[i for i in range(14)if i not in [row['omitted'][0],7+row['omitted'][1]]]
            adjacent=[(i,j)for i,j in combinations(old_A6,2)
                if ratdot(base['seed_vectors'][i],base['seed_vectors'][j])==F(-1,2)]
            assert len(adjacent)==10
            row['affine_A2_seed_pairs']=adjacent
        tri=[[F(1 if i==j else -F(1,2))for j in range(3)]for i in range(3)]
        assert all(sum(row)==0 for row in tri)and tri[0][0]*tri[1][1]-tri[0][1]**2==F(3,4)
        # Finite-contact graph of spherical caps only (vertical seeds excluded).
        parent=list(range(393+len(caps)))
        def find(a):
            while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
            return a
        def join(a,b):parent[find(a)]=find(b)
        if rr2==c(28):
            for ia,ib in combinations(range(392),2):
                if dist(vertices[ia],vertices[ib])==q(84):join(1+ia,1+ib)
        def contact_index(v):
            if v=='central':return 0
            if isinstance(v,str):return 1+int(v.split('_')[1])
            return 393+v
        for a,b,label in chain_contacts:join(contact_index(a),contact_index(b))
        components=Counter(find(a)for a in range(len(parent)))
        # The complete pair gates leave no further incident wall at any of
        # these null rays, except the right endpoint mate of a level8 cap.
        # Indeed a wall outside the triple is orthogonal to its positive null
        # sum iff it is orthogonal to all three roots. All new-cap seed pairs
        # are strictly negative; all cap/cap pairs are strictly negative,
        # except the recorded level8/own-endpoint right relation.
        assert newcounts.get('right',0)==sum(row['level']==8 for row in caps)
        deficient=Counter({2:10*sum(row['level']!=8 for row in caps),
                           3:10*sum(row['level']==8 for row in caps)})
        configurations.append(dict(name=name,wall_count=wallcount,endpoint_radius_square=rr2,
            cap_recipes=caps,old_pair_counts=dict(counts),new_pair_counts=dict(newcounts),
            all_pair_counts=dict(allcounts),finite_contact_paths=chain_contacts,
            spherical_cap_contact_component_sizes=sorted(components.values(),reverse=True),
            actual_facets=wallcount,preserved_old_finite_vertices=392,
            target_probe=dict(multiple=lam,vertex=(0,0,0),cutters=cutters,powers=target_power),
            transported_target_probes_cut=target_coverage,
            strict_escape=dict(multiple=escape_lam,vertex=escape_om,D=escapeD,U=1,norm=1,
                minimum_ball_power=minimum,seed_gap_minimum=1-escape_lam),
            explicit_affine_A2_triples_per_new_cap=10,
            deficient_cusp_counts_by_quotient_rank={rank:n for rank,n in deficient.items()if n},
            deficient_cusp_scope='These explicit positive-null triples only; complete pair gates prove each incident set is the triple plus its level8 right endpoint mate when present. Rank2 or3 is below required19. No census of other cusp rays.',
            scope='Exact pair-compatible actual facets and preserved392finite endpoints. Contact paths are actual codimension-two faces, not a complete vertex/edge census. The saved strict escape proves infinite volume.'))
    paths=['scripts/fresh_h20_balanced_20261014.py',str(Path(__file__).resolve().relative_to(ROOT))]
    return enc(dict(status='PASS_EXACT_GOLDEN_RADIAL_CONTACT_CHAINS_WITH_STRICT_ESCAPES',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        field='Q(phi)(sqrt(12N-12)); phi²=phi+1',vertex_square=N,radicand=D,kappa=k,
        kappa_equation='(N-4)k²-2(N+2)k+(N-4)=0, 0<k<1',
        ray_code=ray_code,minimum_ray_separation_square=min_ray,
        vertex_distance_distribution=[dict(square=d,multiplicity=n)for d,n in sorted(distance_counts.items())],
        signature=[20,1],strict_center=dict(spatial_zero=True,D=200,U=1,norm=-400),
        configurations=configurations,
        central_contact_interior_center_bound=dict(radius_upper_bound=2,ball_containment_radius=6,
            proof='The seven A6 gaps sum to7, so one positive gap is<=1. A positive Coxeter ratio is>=1/2, hence R<=2. Finite/parallel centralR2 contact gives center_length<=R+2, hence outer_radius<=6. Boundary centers are excluded.'),
        scope='Three exactly coherent golden balanced constructions. The longer chain changes the endpoint layer radius and links eight endpoint facets to central by finite contacts. The full second layer cuts all392 transported7z/20 probes. Every configuration still has a strict escape; arbitrary future completion remains open.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status']);print([(r['wall_count'],r['new_pair_counts'])for r in out['configurations']])
