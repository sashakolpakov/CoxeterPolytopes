#!/usr/bin/env python3
"""Exact persistent null frontier and all847 commonescape; bounded graph scout.

The exact conclusions inherit algebraic recipe interval boxes. Numerical
compatibility is reported separately and is never an exact angle certificate.
"""
import argparse,json,math
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from collections import Counter
import numpy as np
import networkx as nx
from fresh_h20_coupled_roles_20261014_field import Field,parse,B
from fresh_h20_s61_multicap_20261014_exact import FILES
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_complementary_frontier_20261016'
BASE='results/fresh_h20_s61_multicap_20261014_exact.json'
BANK='results/fresh_h20_nonright_s61_20261014_scout.json'
GATE='results/fresh_h20_nonright_s61_20261014_bank_gate.json'
ARC='results/fresh_h20_dotted_arc_cap_20261016.json'
CAP='results/fresh_h20_complementary_cap_20261016.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def I(x):return (F(x),F(x))
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def scale(x,t):return (x[0]*t,x[1]*t)if t>=0 else(x[1]*t,x[0]*t)
def mul(x,y):
    v=[a*b for a in x for b in y];return(min(v),max(v))
def root(n):
    den=1<<256;m=math.isqrt(n*den*den)
    return(F(m,den),F(m+1,den))
def enc(v):return list(map(str,v))
def midpoint(v):return float((F(v[0])+F(v[1]))/2)
def run():
    base,bank,gate,arc,cap=[json.loads((ROOT/p).read_text())for p in [BASE,BANK,GATE,ARC,CAP]]
    for doc in [base,bank,gate,arc,cap]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    S=base['support'];J=[[F(x)for x in row]for row in base['support_inverse']]
    ef=list(map(F,cap['strict_escape']['frame_coefficients']));eb=list(map(F,cap['strict_escape']['dual_coefficients']))
    assert ef==[F(1),F(-1,100),F(-49,50)]
    assert eb==[F(-1)if i==17 else F(-1,100000)for i in S]
    wt=[sum(J[i][j]*eb[j]for j in range(18))for i in range(18)]
    norm=-ef[0]**2+ef[1]**2+ef[2]**2+sum(eb[i]*wt[i]for i in range(18))
    assert norm==F(cap['strict_escape']['norm'])==F(1842000049,4000000000)>0
    assert all(F(v[1])<0 for v in cap['strict_escape']['all135_pair_bounds'])
    assert gate['common_center']['inherited_recipe']==arc['common_center']['inherited_recipe']
    assert F(gate['common_center']['maximum_bank_pair_upper'])<0
    profiles={r['id']:r for r in bank['profiles']+bank['skipped']}
    mags={3:I(1),4:root(2),5:scale(add(I(1),root(5)),F(1,2)),6:root(3)}
    exact=[];Y=[];D=[]
    assert len(gate['records'])==847
    for k,r in enumerate(gate['records']):
        assert r['index']==k and r['profile']==bank['proposals'][k]['profile']
        pr=profiles[r['profile']];d=[I(0)]*18
        for i,m in zip(pr['support'],pr['labels']):d[i]=mags[m]
        y=[tuple(map(F,b))for b in r['residual_frame']]
        val=I(0)
        for a,b in zip(y,[-ef[0],ef[1],ef[2]]):val=add(val,scale(a,b))
        for a,b in zip(d,wt):val=add(val,scale(a,-b))
        assert val[1]<0,(k,val)
        assert F(r['common_center_pair'][1])<0
        exact.append(dict(index=k,E1_pair_bounds=enc(val)))
        Y.append([midpoint(v)for v in y]);D.append([midpoint(v)for v in d])
    # Full incident set of nu=sqrt3*r134+2*n0+n1, using positivity of its
    # coefficients and acuteness for every root outside the affine triple.
    old,bridge,edge,polar,cur=[json.loads((ROOT/p).read_text())for p in FILES]
    A=Field(old,bridge,edge)
    M=[tuple(parse(x,1)for x in v)for v in old['normals']]
    n95=tuple(parse(x,2)for x in bridge['new_normal'])
    right=[r['wall']for r in cap['pairs']if r['relation']=='right'];assert len(right)==18
    white=[];excluded=[];rows=[]
    for i in right:
        if i==1:continue
        if i==95:
            g0=A.dot(A.vlift(M[0],1,2),n95,2);g1=A.dot(A.vlift(M[1],1,2),n95,2);level=2
        else:g0=A.dot(M[0],M[i],1);g1=A.dot(M[1],M[i],1);level=1
        s0=A.sign(g0,level);s1=A.sign(g1,level);assert s0<=0 and s1<=0
        if s0==s1==0:white.append(i)
        else:excluded.append(i)
        rows.append(dict(wall=i,pair_with_0_sign=s0,pair_with_1_sign=s1,incident=(s0==s1==0)))
    expected=[3,4,6,7,8,9,10,12,13,14,15,17,31]
    assert white==expected and excluded==[2,61,87,95]
    G=[[A.rational(A.dot(M[i],M[j],1),1)for j in white]for i in white]
    components=[];todo=set(range(len(white)))
    while todo:
        seed=min(todo);comp={seed};stack=[seed];todo.remove(seed)
        while stack:
            i=stack.pop()
            for j in sorted(todo.copy()):
                if G[i][j]!=0:todo.remove(j);comp.add(j);stack.append(j)
        inds=sorted(comp);nodes=[white[i]for i in inds]
        block=[[G[i][j]for j in inds]for i in inds]
        assert block==[[F(2)if i==j else F(-1)if abs(i-j)==1 else F(0)for j in range(len(inds))]for i in range(len(inds))]
        components.append(dict(type='A'+str(len(inds)),walls=nodes,Gram=[[str(x)for x in row]for row in block]))
    assert [len(v['walls'])for v in components]==[2,5,4,1,1]
    assert sum(len(v['walls'])for v in components)+2==15
    assert A.dot(M[0],M[1],1)==A.c(-1,1)
    assert cap['finite_obstructions']['positive_kernel']==['sqrt3','2','1']
    # Bounded selection on the frozen847 recipes; this part uses double precision.
    Y=np.array(Y);D=np.array(D);JJ=np.array(J,float)
    Q=D@JJ@D.T+Y[:,1:]@Y[:,1:].T-Y[:,0,None]*Y[None,:,0]
    y0=np.array([midpoint(v)for v in arc['definition']['frame_bounds']]);d0=np.zeros(18);d0[S.index(31)]=midpoint(arc['definition']['dual_coefficient_bounds'])
    g0=D@JJ@d0-Y[:,0]*y0[0]+Y[:,1:]@y0[1:]
    cut0=-Y[:,0]-.01*Y[:,1]+D@JJ[:,S.index(5)]
    cut1=np.array([midpoint(r['E1_pair_bounds'])for r in exact])
    tol=1e-10
    def good(q):
        if q<-2-tol:return True
        if abs(q)<tol or abs(q+2)<tol:return True
        if q>0 or q<-2:return False
        m=round(np.pi/np.arccos(-q/2))
        return m>=3 and abs(q+2*np.cos(np.pi/m))<tol
    valid=[i for i in range(847)if good(g0[i])];graph=nx.Graph();graph.add_nodes_from(valid)
    for ii,i in enumerate(valid):
        for j in valid[:ii]:
            if good(Q[i,j]):graph.add_edge(i,j)
    cliques=sorted([sorted(v)for v in nx.find_cliques(graph)])
    compatible_both=[i for i in valid if i!=44 and good(Q[i,44])]
    assert len(valid)==765 and not compatible_both
    assert Counter(map(len,cliques))=={1:173,2:3975,3:1}
    assert [v for v in cliques if len(v)==3]==[[100,805,810]]
    priority=sorted([dict(index=i,degree=graph.degree(i),E0_pair=float(cut0[i]),E1_pair=float(cut1[i]),
        profile=bank['proposals'][i]['profile'],anchors=bank['proposals'][i]['anchors'],
        anchor_labels=bank['proposals'][i]['anchor_labels'])for i in valid if cut0[i]>1e-8],
        key=lambda r:(-r['degree'],-r['E0_pair'],r['index']))
    assert len(priority)==476 and max(cut1)<0
    paths=[BASE,BANK,GATE,ARC,CAP,f'scripts/{TAG}.py']+FILES+[
        'scripts/fresh_h20_coupled_roles_20261014_field.py','scripts/fresh_h20_s61_multicap_20261014_exact.py']
    return dict(status='PASS_EXACT_PERSISTENT_NULL_FRONTIER_AND_ALL847_ESCAPE_WITH_BOUNDED_SELECTION',
        sources={p:sha(p)for p in paths},
        exact_bank_escape=dict(recipe_count=847,escape_dependency=CAP,norm=str(norm),records=exact,
            maximum_upper_bound=str(max(F(r['E1_pair_bounds'][1])for r in exact)),
            maximizer=max(range(847),key=lambda i:F(exact[i]['E1_pair_bounds'][1])),
            common_center_preserved=True,
            conclusion='The relaxed intersection of all847exactrecipehalfspaces andold134retainedwalls has a strict spacelike escape and commonstrictcenter. Everycompatible banksubset witharc133 remains infinitevolume.'),
        exact_null_frontier=dict(positive_null_combination={'134':'sqrt3','0':'2','1':'1'},
            incident_walls=[0,1]+white+[134],incident_count=16,other_right_wall_tests=rows,
            affine_component=dict(type='affineG2',walls=[134,0,1],positive_kernel=['sqrt3','2','1'],rank=2),
            finite_components=components,Euclidean_link_dimension=19,normal_rank=15,lineality_dimension=4,
            minimum_additional_incident_walls_for_bounded_link=5,
            future_null_proof='The exactaffineGram has the stated positivekernel; its weighted pairing withcommonfuturecenter isnegative. Alloutsidewallpairings are nonpositive, with equality exactlyon the recordedincidentset.',
            persistence='Any coherentacute addedwall has nonpositive pairing withnu. Itcannotcutnu. Ifincidentatnu,itmustberighttoallthree affineG2roots. Boundingthe4freedirections requiresatleast5 newincidentnormals.'),
        numerical_selection=dict(tolerance=tol,bank_size=847,compatible_with_arc133=valid,
            E0_cutters_among_compatible=len(priority),compatible_with_both_newcaps=compatible_both,
            induced_compatible_edge_count=graph.number_of_edges(),maximal_clique_size_counts={str(k):v for k,v in sorted(Counter(map(len,cliques)).items())},
            maximal_cliques=cliques,E0_cutters_by_compatible_degree=priority,
            scope='Onlynumericalangleclassificationat1e-10; noexactacceptance/exclusion oralllabelseparationclaim. Theexactall847E1gate above doesnotdependonthisgraph.'),
        scope='Exactall847recipeescape andcompleteincidentsetofonefutureidealpoint; notanambientidealvertexcensus. Numericalcompatiblecliquediagnosticseparatelymarked. Noarithmeticity assumption ornewcapsearch.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:
        old=json.loads(p.read_text())
        # Exact fields anddiscrete graphmustreplay; floatingrecordvalues arediagnostic.
        for q in [d,old]:
            for r in q['numerical_selection']['E0_cutters_by_compatible_degree']:
                r['E0_pair']=round(r['E0_pair'],10);r['E1_pair']=round(r['E1_pair'],10)
        assert d==old
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['exact_bank_escape']['recipe_count'],d['exact_null_frontier']['incident_count'],flush=True)
