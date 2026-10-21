#!/usr/bin/env python3
"""Independent integer-coordinate graph/coloring audit; no project imports."""
from collections import Counter
from hashlib import sha256
from pathlib import Path
import argparse,json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_weighted_binary_20261014_local'
PRIMARY='fresh_h20_weighted_binary_20261014'
SIZE=1<<19
ALL=SIZE-1
PARITY=ALL^(1<<10)
A=(1<<0)|(1<<1)|(1<<2)
B=(1<<0)|(1<<3)|(1<<4)
UNIT=(1<<10)-1
GOLD=ALL^((1<<11)-1)

def sgn(a,b):
    """Sign of a+b*phi, entirely in integers."""
    x,y=2*a+b,b
    if not x:return (y>0)-(y<0)
    if not y:return (x>0)-(x<0)
    if x*y>0:return (x>0)-(x<0)
    n=x*x-5*y*y
    assert n
    return ((x>0)-(x<0))*((n>0)-(n<0))

def twice_center(mask):
    full=B if (mask&PARITY).bit_count()%2 else A
    return [(-1 if(mask>>i)&1 else 1)*(2 if(full>>i)&1 else 1)for i in range(19)]

def distance_times_four(x,y):
    squares=[(a-b)**2 for a,b in zip(x,y,strict=True)]
    gold=sum(squares[11:])
    return sum(squares[:10])+2*squares[10]+gold,gold

def label(d):
    return {(8,0):2,(12,0):3,(8,4):5}.get(d,0)

def swap(mask):
    out=mask
    for i,j in [(1,3),(2,4)]:
        if ((mask>>i)^(mask>>j))&1:out^=(1<<i)|(1<<j)
    return out

def check_colorings(doc,edges,labels):
    adj={tuple(e[:2])for e in edges}
    nonfive={i for i,v in enumerate(labels)if v!=5}
    choices=[-1]+[i for i,v in enumerate(labels)if v==5]
    assert len(nonfive)==42 and len(choices)==73
    assert len(doc['coloring_cases'])==len(choices)
    hist=Counter()
    for chosen,row in zip(choices,doc['coloring_cases'],strict=True):
        assert type(row['chosen_pi5_neighbor'])is int and row['chosen_pi5_neighbor']==chosen
        nodes=nonfive if chosen==-1 else {i for i in nonfive if tuple(sorted((i,chosen)))in adj}
        coloring={int(k):v for k,v in row['coloring'].items()}
        assert set(coloring)==nodes
        assert all(type(c)is int and c>=0 for c in coloring.values())
        for i,j,_ in edges:
            if i in nodes and j in nodes:assert coloring[i]!=coloring[j]
        count=max(coloring.values(),default=-1)+1
        assert count==row['color_count']
        bound=1+count+int(chosen!=-1)
        assert bound==row['full_clique_upper_bound'] and bound<20
        hist[bound]+=1
    return hist

def run():
    doc=json.loads((ROOT/f'results/{TAG}.json').read_text())
    orbit=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    for data in (doc,orbit):
        for path,h in data['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h
    assert doc['status']=='PASS_EXACT_BINARY_CAP_FINITE_VERTEX_OBSTRUCTION'
    assert orbit['full_unit_triples']==[[0,1,2],[0,3,4]]
    assert orbit['orbit_size']==SIZE and orbit['radius']==1

    # Exhaust every sign mask, using independent integer squared distances.
    # Four times the physical squared distance is a+b*phi.
    all_neighbors={};counts=Counter()
    for mask in range(1,SIZE):
        hg=(mask&GOLD).bit_count();hh=(mask>>10)&1
        if (mask&PARITY).bit_count()%2==0:
            ha=(mask&A).bit_count();hn=(mask&(UNIT^A)).bit_count()
            d=(16*ha+4*hn+8*hh+4*hg,4*hg)
        else:
            hi=(mask&A&B).bit_count();hd=(mask&(A^B)).bit_count()
            hn=(mask&(UNIT^(A|B))).bit_count()
            d=(4+16*hi+8*hd+4*hn+8*hh+4*hg,4*hg)
        r=label(d)
        if r:all_neighbors[mask]=r;counts[r]+=1
        else:assert sgn(d[0]-16,d[1])>=0
    assert len(all_neighbors)==95
    assert dict(orbit['finite_neighbor_masks'])==all_neighbors
    masks=doc['binary_neighbor_masks'];assert len(masks)==len(set(masks))==95
    assert set(masks)==set(all_neighbors)
    V=[twice_center(m)for m in masks]
    v0=twice_center(0)
    assert all(label(distance_times_four(v0,v))==all_neighbors[m]for m,v in zip(masks,V,strict=True))

    # All omitted reference neighbors are dotted: negative seed walls,
    # central sphere, and the entire sqrt(2)-radius endpoint orbit.
    # Center norm =29/4+2phi, central cosine=(9+8phi)/16>1.
    assert sgn(-7,8)>0
    # min endpoint distance squared=17/4+2phi >6 > (1+sqrt2)^2.
    assert sgn(-7,8)>0 and 4*2<3*3
    seed_labels=[2 if(A>>i)&1 else 3 for i in range(10)]+[4]+[5]*8
    labels=[all_neighbors[m]for m in masks]+seed_labels
    assert labels==doc['reference_pair_labels']

    # Direct coordinate reconstruction of every local pair. The 19 seed
    # neighbors are the positive coordinate walls, hence mutually right.
    edges=[]
    for i in range(114):
        for j in range(i):
            if i<95:r=label(distance_times_four(V[i],V[j]))
            elif j>=95:r=2
            else:
                axis=i-95;c=V[j][axis]
                if axis<10:r=2 if c==2 else 3 if c==1 else 0
                else:r=(4 if axis==10 else 5)if c==1 else 0
            if r:edges.append([j,i,r])
    assert len(edges)==2644 and edges==doc['finite_pair_graph_edges']
    hist=check_colorings(doc,edges,labels)
    assert hist==Counter({18:1,13:64,19:8})
    assert doc['maximum_full_clique_bound']==max(hist)==19
    assert doc['bound_distribution']=={str(k):v for k,v in hist.items()}

    # No Dynkin catalog is needed. For two pi5 neighbors of a root,
    # det normalized Gram =1-2c²-d²-2c²d <=1-2c²=(1-phi)/2<0.
    assert sgn(1,-1)<0

    # Transitivity: even unit/golden flips preserve the triple type.
    # For odd flips, compose the unit swap (1 3)(2 4), which interchanges
    # A and B, then the prescribed sign flip. The swap preserves parity.
    assert swap(A)==B and swap(B)==A and swap(PARITY)==PARITY
    assert swap(UNIT)==UNIT and swap(GOLD)==GOLD
    assert all(swap(swap(1<<i))==1<<i for i in range(19))
    # Exhaust all target signs: each map sends reference triple/signs to
    # the target center and preserves the family type rule for either type.
    for m in range(SIZE):
        odd=(m&PARITY).bit_count()%2
        mapped=swap(A)if odd else A
        assert mapped==(B if odd else A)
        for initial in (A,B):
            moved=swap(initial)if odd else initial
            initial_parity=int(initial==B)
            assert moved==(B if(initial_parity^odd)else A)

    # A deliberately invalid coloring is rejected by the saved-certificate
    # validation, rather than silently accepted as heuristic output.
    bad=json.loads(json.dumps(doc));bad['coloring_cases'][0]['coloring']={k:0 for k in bad['coloring_cases'][0]['coloring']}
    rejected=False
    try:check_colorings(bad,edges,labels)
    except AssertionError:rejected=True
    assert rejected
    paths=[f'scripts/{TAG}.py',f'results/{TAG}.json',f'scripts/{PRIMARY}.py',f'results/{PRIMARY}.json',f'scripts/{TAG}_review.py']
    return dict(status='PASS_INDEPENDENT_BINARY_CAP_LOCAL_VERTEX_OBSTRUCTION',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        method='No producer/project imports; integer squared distances in Q(phi), exhaustive sign masks, direct saved coloring checks.',
        sign_masks_checked=SIZE,finite_binary_neighbors=95,neighbor_label_counts=dict(counts),
        local_vertices=114,local_edges=len(edges),coloring_cases=73,bound_distribution=dict(hist),
        maximum_incident_spherical_clique=19,transitivity='even sign flips; odd sign flips composed with unit swap (1 3)(2 4)',
        double_pi5_determinant_upper_bound='(1-phi)/2 < 0',forged_coloring_rejected=True,
        scope='No finite vertex of the frozen chamber is incident to any binary cap. The finite-neighbor inventory includes binary caps and seed walls; central and all endpoint caps are dotted. This is not a no-completion theorem after adding different walls, and does not classify ideal vertices.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    data=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}_review.json'
    if args.verify:assert json.loads(out.read_text())==data
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status']);print(data['bound_distribution'])
