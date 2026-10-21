#!/usr/bin/env python3
"""Stdlib native K arithmetic and saved one-pair all89-cover accounting.

Independent arithmetic on the saved literal frame; historical owner/native
embedding completeness remain the bound producer and predecessor inputs.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_raw4_all89_retained_pair_barrier_20261018'
BARRIER='fresh_new5121_raw4_all89_cover_barrier_20261018'
RAW='fresh_new5121_raw4_next_four_domain_20261018'
BASE='fresh_new5121_raw4_two_face_local_covers_20261018'
OLD='fresh_cross_vertex_current290_local_covers_20261018'
MODEL='fresh_cross_vertex_current_full_model_20261018'


def k(a=0,b=0):return F(a),F(b)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def sub(x,y):return add(x,neg(y))
def mul(x,y):return x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def div(x,y):
    norm=y[0]*y[0]-5*y[1]*y[1];assert norm
    a,b=mul(x,(y[0],-y[1]));return a/norm,b/norm

def sign(x):
    a,b=x
    if not b:return (a>0)-(a<0)
    if not a:return (b>0)-(b<0)
    if (a>0)==(b>0):return (a>0)-(a<0)
    c=a*a-5*b*b;assert c
    return ((a>0)-(a<0))*((c>0)-(c<0))

def kout(x):return list(map(str,x))
def integral(x):return (2*x[0]).denominator==1 and (x[0]*x[0]-5*x[1]*x[1]).denominator==1
@lru_cache(None)
def read(tag):return json.loads((ROOT/f'results/{tag}.json').read_text())
def jhash(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def run():
    d,b,r,c=map(read,(TAG,BARRIER,RAW,BASE))
    assert d['status']=='PASS_ALL89_CURRENT_PAIR_GATES_RETAIN_ABOVE12_WITNESS'
    for path,h in d['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    assert (d['selected_vertex'],d['selected_face'])==(495,14368)
    raw=next(x for x in r['forest_records'] if x['vertex']==495)
    inv=next(x for x in r['cover_inventory']['forests'] if x['vertex']==495)
    assert d['target']==raw['target'] and d['forest']==raw['forest']
    assert d['baseline']==inv['baseline'] and d['complete_selected_occurrences']==inv['actual_occurrences']
    assert d['raw_star_indices']==[0,9] and d['merged_star_indices']==[306,315]
    assert len(d['selected_native_stars'])==len(d['complete_selected_orbits'])==2
    for rawi,merged,saved in zip(d['raw_star_indices'],d['merged_star_indices'],d['selected_native_stars']):
        assert saved==dict(index=merged,**raw['literal_stars'][rawi]['native'])
        assert inv['raw4_to_merged_star_indices'][rawi]==merged
    assert [o['orbit'] for o in d['complete_selected_occurrences']]==[13,14]
    for gi,z,o in zip((13,14),d['complete_selected_orbits'],d['complete_selected_occurrences']):
        assert z['index']==gi and z['literal_stars']==o['literal_star_indices']
        assert z['origins']==[o['origin']] and z['positive_face_coefficients']==[[14368,1]]
        assert len(z['literal_stars'])==8
    assert not set(d['complete_selected_orbits'][0]['literal_stars'])&set(d['complete_selected_orbits'][1]['literal_stars'])
    assert all(i>=306 for z in d['complete_selected_orbits'] for i in z['literal_stars'])
    stars=d['selected_native_stars']
    assert not set(map(tuple,stars[0]['occupied_pairs']))&set(map(tuple,stars[1]['occupied_pairs']))
    assert stars[0]['contacts'][0]==stars[1]['contacts'][0]
    walls=[stars[0]['contacts'][0],stars[0]['contacts'][1],stars[1]['contacts'][1]]
    assert len({json.dumps(x) for x in walls})==3 and walls==d['distinct_external_contacts']
    # Rebuild all13 normalized Gram entries and invert using independent K pairs.
    tau=k(F(7,2),F(1,2));phi=k(F(1,2),F(1,2));n=13;target={int(i):{int(j):e for j,e in v.items()} for i,v in d['target'].items()}
    assert set(target)==set(range(n))
    G=[[k(2) if i==j else neg(k(1) if target[i].get(j)==3 else phi) if j in target[i] else k() for j in range(n)] for i in range(n)]
    assert all(e in (3,5) for row in target.values() for e in row.values())
    assert G==list(map(list,zip(*G)))
    A=[G[i][:]+[k(int(i==j)) for j in range(n)] for i in range(n)]
    for i in range(n):
        pivot=A[i][i];assert sign(pivot)>0
        A[i]=[div(x,pivot) for x in A[i]]
        for j in range(n):
            if j!=i:
                q=A[j][i];A[j]=[sub(x,mul(q,y)) for x,y in zip(A[j],A[i])]
    J=[row[n:] for row in A]
    assert all(A[i][j]==k(int(i==j)) for i in range(n) for j in range(n))
    def pmul(x,y):return add(mul(x[0],y[0]),mul(tau,mul(x[1],y[1]))),add(mul(x[0],y[1]),mul(x[1],y[0]))
    def contact(edge):
        label,a,z=edge
        if label==0:assert k(a,z)==tau;return k(),k(-1)
        assert label in (3,5);return (k(-1) if label==3 else neg(phi)),k()
    vectors=[]
    for wall in walls:
        vec=[(k(),k()) for i in range(n)]
        assert len({i for i,e in wall})==len(wall)
        for i,e in wall:vec[i]=contact(e)
        vectors.append(vec)
    def projection(a,b):
        x,y=k(),k()
        for i in range(n):
            for j in range(n):
                p,q=pmul(vectors[a][i],vectors[b][j]);x=add(x,mul(p,J[i][j]));y=add(y,mul(q,J[i][j]))
        return x,y
    hs=[]
    for i in range(3):
        p,q=projection(i,i);assert q==k();h=sub(p,k(2));assert sign(h)>0;hs.append(h)
    assert list(map(kout,hs))==d['squared_external_heights']
    assert hs==[k(-F(1,4),F(1,4)),k(F(3,2),1),k(1,1)]
    for star,ids in zip(stars,((0,1),(0,2))):assert star['squared_timelike_heights']==[kout(hs[i]) for i in ids]
    S=k(F(13,2),F(5,2));products={(0,1):k(),(0,2):k(-1)}
    assert projection(0,1)==(k(),k(F(1,2))) and mul(hs[0],hs[1])==div(tau,k(4))
    assert projection(0,2)==(k(),k()) and mul(hs[0],hs[2])==k(1)
    assert projection(1,2)==(k(),k()) and mul(hs[1],hs[2])==S and sign(sub(S,k(4)))>0
    assert integral(S) and 2*S[0]==13 and S[0]*S[0]-5*S[1]*S[1]==11
    assert mul(tau,S)==mul(k(3,2),k(3,2)) and sign(k(3,2))>0
    assert d['explicit_products']==dict(AB='0',AC='-1',BC_negative_square=kout(S),tau=kout(tau),cycle_positive_K_root=kout(k(3,2)))
    assert [x['walls'] for x in d['complete_distinct_wall_products']]==[[0,1],[0,2],[1,2]]
    for row,sg,sg1 in zip(d['complete_distinct_wall_products'],[0,-1,-1],[1,0,-1]):
        i,j=row['walls'];p,q=projection(i,j)
        assert row['projection']==[kout(p),kout(q),kout(k()),kout(k())]
        assert row['pair_field_radicals'][0]==kout(tau)
        assert row['height_product']==kout(mul(hs[i],hs[j]))
        assert (row['sign'],row['sign_plus_one'],row['arithmetic_gate_code'],row['arithmetic_gate_reason'])==(sg,sg1,None,None)
    # Thus AB=0, AC=-1 and BC<-2 with integral K square: all stated pair gates retain.
    assert d['pair_outcome']['newly_owned_pair_domain']==[[306,315,None]]
    assert d['pair_outcome']['new_exclusions']==0 and d['pair_outcome']['new_unknown_pairs']==1
    old=next(x for x in read(OLD)['forest_records'] if x['vertex']==495)
    origin=next(x for x in read(old['source'])['forest_records'] if x['vertex']==495)
    vv=origin['maximal_count_vectors'];assert jhash(vv)==old['count_vectors_sha256']
    cert=old['certificate'];weights=list(map(F,cert['weights']));dem=list(map(F,old['demands']))
    assert len(weights)==len(dem)==13 and all(w>=0 for w in weights)
    assert len(vv)==old['complete_count_vectors']==19 and all(len(v)==13 for v in vv)
    for v in vv:assert sum((w*x for w,x in zip(weights,v)),F())<=1
    covered=[F() for i in range(13)];cost=F()
    for index,m in cert['mixture']:
        m=F(m);assert m>=0;cost+=m
        for i,x in enumerate(vv[index]):covered[i]+=m*x
    assert all(x>=y for x,y in zip(covered,dem))
    assert cost==sum((x*y for x,y in zip(weights,dem)),F())==F(old['exact_local_optimum'])
    assert jhash(old)==d['baseline']['record_sha256']
    alpha={f:F(x) for f,x in c['face_multipliers']};mu={v:F(x) for v,x in c['vertex_multipliers']}
    mass=alpha[14368];assert str(mass)==d['pair_cover_mass']==d['saving_after_multiplier_slack']
    assert d['face_multipliers']==c['face_multipliers'] and d['initial_vertex_multipliers']==c['vertex_multipliers']
    before={x['vertex']:x for x in b['all89']['forests']};after=d['all89_cover_upper_records']
    assert len(after)==len(before)==238 and {x['vertex'] for x in after}==set(before)
    gain=F();upper=dict(mu)
    for row in after:
        v=row['vertex'];oldrow=before[v];save=mass if v==495 else F();assert F(row['pair_cover_mass'])==save
        for key in oldrow:
            if key not in ('possible_vertex_mass_upper','maximum_possible_mu_increase'):assert row[key]==oldrow[key]
        m,t,e=map(F,(oldrow['baseline_vertex_mass'],oldrow['baseline_local_cover'],oldrow['actual_corner_demand']))
        assert e>=2*save;u=max(m,t+e-save);gain+=u-m;upper[v]=u
        assert F(row['possible_vertex_mass_upper'])==u and F(row['maximum_possible_mu_increase'])==u-m
        assert F(row['saving_after_multiplier_slack'])==F(oldrow['maximum_possible_mu_increase'])-(u-m)
    assert F(b['all89']['maximum_possible_mu_increase'])-gain==mass
    N,M=F(b['objective_numerator']),F(b['current_normalizer']);den=sum(upper.values(),F())
    assert den==M+gain and all(upper.get(v,F())>=x for v,x in mu.items())
    assert d['upper_vertex_multipliers']==[[v,str(x)] for v,x in sorted(upper.items()) if x]
    assert F(d['maximum_possible_mu_increase'])==gain and F(d['normalizer_upper'])==den
    assert F(d['ratio_lower_bound'])==N/den==F(8222765954018112,684963831123457)>12
    assert F(d['strict_numerator_gap_over12'])==N-12*den>0
    two=0
    for col in read(MODEL)['prepared']['new_columns']:
        load=sum((alpha.get(f,F())*n for f,n in col['integer_face_coefficients']),F())
        oldcost=sum((mu.get(v,F())*n for v,n in col['integer_vertex_capacities']),F())
        newcost=sum((upper.get(v,F())*n for v,n in col['integer_vertex_capacities']),F())
        assert load<=oldcost<=newcost
        two+=len(col['integer_vertex_capacities'])>1
    assert two==14
    assert d['pair_tests']==1 and all(d[x]==0 for x in ('new_clique_enumerations','covering_LPs','global_model_optimizations'))
    assert (d['authoritative_classes'],d['authoritative_forests'],d['new_dimension_bound'])==(324,600,None)
    print('PASS_RETAINED_PAIR_SAVED_CHECK',dict(frame=13,distinct_external_walls=3,wall_products=3,old_charge_inequalities=19,old_cover_terms=len(cert['mixture']),all89_forests=238,two_owner_columns=two,ratio=d['ratio_lower_bound']))

if __name__=='__main__':run()
