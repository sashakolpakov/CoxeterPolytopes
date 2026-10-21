#!/usr/bin/env python3
"""Exact physical remote U extension of the original 63-wall chamber.

Directed rational intervals enclose the original physical roots. The new
root is defined by a positive-definite three-frame and its timelike center
complement, not by the numerical N/S/T construction.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
import math
from itertools import combinations
from pathlib import Path
import fresh_oblique_low_s_continuous_n_interval_20261018 as E

P=E.P; I=P.I; ROOT=P.ROOT
TAG='fresh_remote60_exact_original_extension_20261018'
PRIMARY='fresh_nonarithmetic_release49_62_replacement_20261017'
DISCOVERY='fresh_remote60_compact_face_bridge_20261018'
V,add,sub,mul,div,neg,sqrt,dot=P.V,P.add,P.sub,P.mul,P.div,P.neg,P.sqrt,I.dot
WHITE=I.WHITE; FRAME=[60,58,59]; NEW=75


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v)for v in x]
    return x


def refine(z,i,j):
    identity=None; label=None
    if i==j:identity=V(2);label='norm2'
    elif -1<z[0]<=0<=z[1]<1:identity=V(0);label=2
    elif -2<z[0]<=z[1]<0:
        m=round(math.pi/math.acos(-float(sum(z)/2)/2))
        if m>=3 and z[0]>-E.theta(m+1)[0]and z[1]<-E.theta(m-1)[1]:
            identity=neg(E.theta(m));label=m
    if identity is None:return z,None
    intersection=max(z[0],identity[0]),min(z[1],identity[1])
    assert intersection[0]<=intersection[1]
    return intersection,dict(pair=[i,j],raw=z,identity=label,identity_interval=identity,intersection=intersection)


def vector_sub(a,b):return [sub(x,y)for x,y in zip(a,b)]
def linear(coeff,vecs):return [I.total(mul(c,v[k])for c,v in zip(coeff,vecs))for k in range(len(vecs[0]))]
def row_times(v,g):return [I.total(mul(v[k],g[k][j])for k in range(len(v)))for j in range(len(g))]


def run():
    primary=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    discovery=json.loads((ROOT/f'results/{DISCOVERY}.json').read_text())
    for d in(primary,discovery):
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert primary['total_wall_count']==63 and primary['complete_finite_vertex_count']==199
    metric=[list(map(V,r))for r in primary['computational_gram']]
    center=list(map(V,primary['common_center']))
    roots={r['wall']:[div(mul(V(2),sub(x,y)),V(primary['center_products'][str(r['wall'])]))
             for x,y in zip(center,map(V,r['point']))]for r in primary['facet_witnesses']}
    # The original last root has a direct coordinate enclosure much tighter
    # than reconstruction through its coarse center/facet-witness products.
    assert primary['new_wall']==66
    roots[66]=list(map(V,primary['new_normal']))
    ids=sorted(roots);assert len(ids)==63 and len(WHITE)==6 and len(set(ids)-set(WHITE))==57
    hc={i:V(primary['center_products'][str(i)])for i in ids}
    assert all(z[1]<0 for z in hc.values())
    h2=V(primary['center_norm']);assert h2[1]<0
    native={};identities=[]
    for i in WHITE+FRAME:
        native[i]={}
        for j in ids:
            raw=dot(roots[i],metric,roots[j]);z,proof=refine(raw,i,j)
            native[i][j]=z
            if proof:identities.append(proof)
    wi,wpiv=I.inverse([[native[i][j]for j in WHITE]for i in WHITE]);assert all(z[0]>0 for z in wpiv)
    wc={i:[native[w][i]for w in WHITE]for i in ids}
    hw=[hc[w]for w in WHITE]
    hp={i:sub(hc[i],dot(hw,wi,wc[i]))for i in ids}
    hp2=sub(h2,dot(hw,wi,hw));assert hp2[1]<0
    projected={i:{j:sub(native[i][j],dot(wc[i],wi,wc[j]))for j in ids}for i in FRAME}
    gf=[[projected[i][j]for j in FRAME]for i in FRAME]
    inv,piv=I.inverse(gf);assert all(z[0]>0 for z in piv)
    assert next(r['identity']for r in identities if r['pair']==[58,59])==9
    assert all(wc[i][k]==V(0)for i in(58,59)for k in range(6))
    rhs=list(map(V,[0,0,-1]));hf=[hp[i]for i in FRAME]
    remainder=sub(V(2),dot(rhs,inv,rhs));q=sub(hp2,dot(hf,inv,hf))
    assert remainder[1]<0 and q[1]<0
    factor=sqrt(div(remainder,q));coeff=row_times(rhs,inv);hcoeff=row_times(hf,inv)
    products={}
    for j in ids:
        cj=[projected[i][j]for i in FRAME]
        products[j]=add(dot(rhs,inv,cj),mul(factor,sub(hp[j],dot(hf,inv,cj))))
    # The white orthogonality and frame contacts are defining exact equations.
    for j in WHITE:products[j]=V(0)
    for j,z in zip(FRAME,rhs):products[j]=z
    center_product=add(dot(rhs,inv,hf),mul(factor,q))
    opposite_center=sub(dot(rhs,inv,hf),mul(factor,q))
    assert center_product[1]<0 and opposite_center[0]>0
    dotted=set(ids)-set(WHITE+FRAME)
    assert len(dotted)==54
    wcoeff={j:row_times(wc[j],wi)for j in ids}
    projroots={j:vector_sub(roots[j],linear(wcoeff[j],[roots[w]for w in WHITE]))for j in ids}
    ph=vector_sub(center,linear(row_times(hw,wi),[roots[w]for w in WHITE]))
    residual=vector_sub(ph,linear(hcoeff,[projroots[j]for j in FRAME]))
    u=[add(a,mul(factor,b))for a,b in zip(linear(coeff,[projroots[j]for j in FRAME]),residual)]
    unorm=dot(u,metric,u);assert unorm[0]<=2<=unorm[1]
    direct_products={j:dot(u,metric,roots[j])for j in ids}
    assert all(max(z[0],direct_products[j][0])<=min(z[1],direct_products[j][1])for j,z in products.items())
    schur_products=products.copy()
    products={j:(max(z[0],direct_products[j][0]),min(z[1],direct_products[j][1]))for j,z in products.items()}
    assert all(products[j][1]<-2 for j in dotted),[(j,enc(products[j]))for j in dotted if products[j][1]>=-2]
    facet_point=vector_sub(center,[mul(div(center_product,V(2)),z)for z in u])
    facet_square=sub(h2,div(mul(center_product,center_product),V(2)))
    facet_products={i:sub(hc[i],mul(div(center_product,V(2)),products[i]))for i in ids}
    old_witness_products={i:sub(center_product,mul(div(hc[i],V(2)),products[i]))for i in ids}
    assert facet_square[1]<0 and all(z[1]<0 for z in facet_products.values())and all(z[1]<0 for z in old_witness_products.values())
    # The three vertices of the old60 face in the fixed-white H3 section.
    allroots={**projroots,NEW:u};g={}
    for i in FRAME+[NEW]:
        g[i]={}
        for j in ids+[NEW]:
            g[i][j]=V(2)if i==j==NEW else products[j]if i==NEW else products[i]if j==NEW else projected[i][j]
    hp[NEW]=center_product
    vertices=[]
    for sidepair in([58,59],[58,NEW],[59,NEW]):
        support=[60]+sidepair;gv=[[g[i][j]for j in support]for i in support]
        vi,vp=I.inverse(gv);assert all(z[0]>0 for z in vp)
        cv=[hp[i]for i in support];vc=row_times(cv,vi)
        vv=vector_sub(ph,linear(vc,[allroots[i]for i in support]))
        square=sub(hp2,dot(cv,vi,cv));assert square[1]<0
        vpds={}
        for j in ids+[NEW]:
            z=sub(hp[j],dot(cv,vi,[g[i][j]for i in support]))
            if j in WHITE+support:z=V(0)
            else:assert z[1]<0,(support,j,enc(z))
            vpds[j]=z
        vertices.append(dict(section_support=support,ambient_support=WHITE+support,Gram=gv,pivots=vp,inverse=vi,
            coefficients=vc,point=vv,square=square,products=vpds))
    assert len(vertices)==3
    # A finite vertex containing U can use only the nine nondotted old
    # neighbors. Any trial containing both58 and59 has this indefinite
    # principal submatrix, and the two other trials are the vertices above.
    triangle_gram=[[g[i][j]for j in(58,59,NEW)]for i in(58,59,NEW)]
    triangle_determinant=P.determinant(triangle_gram);assert triangle_determinant[1]<0
    trials=[]
    new_supports=[sorted(v['ambient_support'])for v in vertices if NEW in v['ambient_support']]
    for subset in combinations(sorted(WHITE+FRAME),8):
        support=sorted([NEW,*subset])
        if 58 in subset and 59 in subset:
            trials.append(dict(support=support,indefinite_principal=[58,59,NEW]))
        else:
            assert support in new_supports
            trials.append(dict(support=support,certified_finite=True))
    assert len(trials)==9 and sum('certified_finite'in r for r in trials)==2
    escape=primary['strict_escape'];ep=list(map(V,escape['point']))
    escape_product=dot(u,metric,ep)
    assert V(escape['norm'])[0]>0 and V(escape['future'])[1]<0 and escape_product[1]<0
    assert len(escape['products'])==63 and all(V(z)[1]<0 for z in escape['products'].values())
    print('EXACT_U_PASS',len(ids),'original walls;',len(dotted),'strict dotted;',len(vertices),'finite section vertices',flush=True)
    sources={f'{folder}/{stem}.{suffix}':sha256((ROOT/f'{folder}/{stem}.{suffix}').read_bytes()).hexdigest()
        for stem,folders in [(TAG,[('scripts','py')]),(PRIMARY,[('scripts','py'),('results','json'),('results','md')]),
            (DISCOVERY,[('scripts','py'),('results','json'),('results','md')]),
            (E.TAG,[('scripts','py'),('results','md')]),(P.TAG,[('scripts','py'),('results','md')]),
            (I.TAG,[('scripts','py'),('results','md')])]
        for folder,suffix in folders}
    return dict(status='PASS_EXACT_REMOTE60_ORIGINAL_EXTENSION',sources=sources,old_walls=ids,new_wall=NEW,
        deleted_original_walls=[],fixed_white_walls=WHITE,frame=FRAME,profile=0,labels=[2,2,3],
        physical_input=dict(metric=metric,center=center,roots=roots,center_products=hc,center_square=h2),
        native_isolation=identities,white_inverse=wi,white_pivots=wpiv,white_covectors=wc,
        projected_center_products=hp,projected_center_square=hp2,projected_frame_products=projected,
        frame_Gram=gf,frame_inverse=inv,frame_pivots=piv,rhs=rhs,remainder=remainder,
        center_complement_square=q,positive_factor=factor,normal=u,direct_norm=unorm,
        direct_products=direct_products,schur_products=schur_products,native_products=products,center_product=center_product,
        opposite_sheet_center_product=opposite_center,
        facet_witness=dict(point=facet_point,square=facet_square,products=facet_products),
        retained_old_facet_witness_products=old_witness_products,section_vertices=vertices,
        new_finite_vertex_trials=trials,triangle_Gram=triangle_gram,triangle_determinant=triangle_determinant,
        retained_strict_escape=dict(inherited=escape,new_product=escape_product),
        summary=dict(original_walls=63,total_walls=64,strict_dotted_original_contacts=54,
            right_original_contacts=8,simple_original_contacts=1,retained_old_finite_vertices=199,
            added_certified_finite_vertices=2,complete_finite_vertices=201,
            compact_section_face=60,section_face_angles=[9,2,3]),
        scope=dict(exact_physical_extension=True,numerical_chain_required=False,
            complete_new_finite_vertex_census=True,ambient_facet60_compact=False,
            infinite_volume_proved=True,finite_volume_proved=False,new_dimension_example_proved=False))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';out=enc(run())
    if args.write:
        assert not path.exists(),'refusing to overwrite an existing certificate'
        path.write_text(json.dumps(out,indent=2)+'\n')
    else:assert json.loads(path.read_text())==out
    print('PASS',TAG,flush=True)

if __name__=='__main__':main()
