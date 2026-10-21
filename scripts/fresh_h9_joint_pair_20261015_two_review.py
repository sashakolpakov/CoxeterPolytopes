#!/usr/bin/env python3
"""Independent two-cap reconstruction and two-step complete 30-to-32 cone clipping.

The initial complete ray incidence list is the previously independently reviewed
30-wall cone. Every starting vector is reconstructed from its nullspace. New-ray
adjacency uses the complete face ray set, not the producer's rank criterion.
"""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_h9_quartet_polar_extension_20261014_review as R

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_joint_pair_20261015_two_review'
BASE='results/fresh_h9_nine_contact_extension_20261015.json'
CAP='results/fresh_h9_joint_pair_20261015_twocaps.json'
PRIMARY='results/fresh_h9_joint_pair_20261015_two_extension.json'
OLD_REVIEW='results/fresh_h9_nine_contact_extension_20261015_review.json'
K,Z,O,C=R.K,R.Z,R.O,R.C

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(v):
    if isinstance(v,(tuple,list)):return [enc(x)for x in v]
    return list(map(str,v.to_list()))


def clip_saved_cone(b,N,Q,H):
    functionals=[C.mv(Q,n)for n in N]
    def normalized(v):
        first=next(C.dot(v,row)for row in functionals[:30]if C.dot(v,row))
        assert R.sg(first)<0
        return tuple(x/-first for x in v)
    rays=[];old_counts=Counter()
    for kind,sets in b['complete_census'].items():
        for ids in sets:
            v=C.null_one([functionals[i]for i in ids])
            vals=[C.dot(v,row)for row in functionals[:30]]
            if R.sg(next(g for g in vals if g))>0:v=[-x for x in v];vals=[-g for g in vals]
            assert all(R.sg(g)<=0 for g in vals)
            assert ids==[i for i,g in enumerate(vals)if not g]
            q=C.pair(v,v,Q);time=C.pair(v,H,Q)
            actual='hyperideal_rays'if R.sg(q)>0 else'finite_vertices'if R.sg(q)<0 else'cusps'
            if R.sg(q)<=0 and R.sg(time)>=0:actual='opposite_sheet_rays'
            assert actual==kind
            rays.append((normalized(v),sum(1<<i for i in ids)))
            old_counts[kind]+=1
        print('reconstructed old',kind,len(sets),flush=True)
    assert len(rays)==894 and len({v for v,_ in rays})==894
    stages=[]
    for wall in range(30,len(N)):
        values=[C.dot(v,functionals[wall])for v,_ in rays]
        neg=[i for i,g in enumerate(values)if R.sg(g)<0]
        pos=[i for i,g in enumerate(values)if R.sg(g)>0]
        zero=[i for i,g in enumerate(values)if not g]
        new={rays[i][0]:rays[i][1]for i in neg}
        new.update({rays[i][0]:rays[i][1]|(1<<wall)for i in zero})
        crossing=0
        for i in neg:
            for j in pos:
                common=rays[i][1]&rays[j][1]
                if common.bit_count()<8:continue
                if any(k not in(i,j)and mask&common==common for k,(_,mask)in enumerate(rays)):continue
                crossing+=1
                v=normalized([values[j]*a-values[i]*b for a,b in zip(rays[i][0],rays[j][0])])
                active=common|(1<<wall)
                if v in new:assert new[v]==active
                new[v]=active
        stages.append(dict(wall=wall,negative=len(neg),positive=len(pos),zero=len(zero),
                           crossing_edges=crossing,resulting_rays=len(new)))
        print('clip',stages[-1],flush=True)
        rays=list(new.items())
    result={kind:{}for kind in ('finite_vertices','cusps','hyperideal_rays','opposite_sheet_rays')}
    for v,active in new.items():
        vals=[C.dot(v,row)for row in functionals]
        assert all(R.sg(g)<=0 for g in vals)
        ids=tuple(i for i,g in enumerate(vals)if not g)
        assert active==sum(1<<i for i in ids)
        assert len(C.eliminate([functionals[i]for i in ids])[1])==9
        q=C.pair(v,v,Q);time=C.pair(v,H,Q)
        kind='hyperideal_rays'if R.sg(q)>0 else'finite_vertices'if R.sg(q)<0 else'cusps'
        if R.sg(q)<=0 and R.sg(time)>=0:kind='opposite_sheet_rays'
        assert ids not in result[kind];result[kind][ids]=v
    return result,dict(old_counts),stages


def run():
    paths=[BASE,CAP,PRIMARY,OLD_REVIEW]
    docs=[json.loads((ROOT/p).read_text())for p in paths]
    for doc in docs:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    b,p,d,old_review=docs
    assert old_review['status']=='PASS_INDEPENDENT_FULL_30_WALL_CONE'
    assert old_review['sources'][BASE]==sha(BASE)
    assert old_review['counts']==b['counts']
    Q,N,H=R.dm(d['form']),R.dm(d['normals']),R.dv(d['center'])
    assert Q==R.dm(b['form'])and H==R.dv(b['center'])and len(N)==32
    assert len(p['caps'])==2
    assert N==R.dm(b['normals'])+[R.dv(row['normal'])for row in p['caps']]
    pair=lambda a,b:C.pair(a,b,Q)
    cov=[C.mv(Q,a)for a in N];norms=[pair(a,a)for a in N]
    assert norms==R.dv(d['norms'])and all(R.sg(q)>0 for q in norms)
    assert norms[-2:]==[2*O,2*O]
    G=[[C.dot(a,row)for row in cov]for a in N];assert G==R.dm(d['gram'])
    ids=[0,1,2,3,5,6,7,10,29]
    kernel=C.null_one([cov[i]for i in ids]);kernel=[x/kernel[9]for x in kernel]
    assert len(C.eliminate([cov[i]for i in ids])[1])==9
    expected_contacts=[{7:5,29:3},{3:5,5:3,7:5}]
    beta=[['30471/18764','7413/18764'],['7236/4691','4059/4691']]
    delta=[['3611919213/88021924','8201201589/440109620'],
           ['1748292174/22005481','3918049056/110027405']]
    phi=R.dec([['1/2','1/2'],['0','0']]);oldE=R.dv(b['strict_escape']['point'])
    accepted_branches=[]
    for j,row in enumerate(p['caps']):
        rec=row['reconstruction'];n=N[30+j];part=R.dv(rec['particular'])
        v=R.dv(rec['kernel']);root=R.dec(rec['positive_square_root'])
        assert v==kernel and pair(part,v)==Z
        assert root==R.dec([['0','0'],beta[j]])and R.sg(root)>0
        assert R.dec(rec['delta'])==R.dec([delta[j],['0','0']])
        assert root*root==R.dec(rec['delta'])==(2*O-pair(part,part))/pair(v,v)
        contacts=expected_contacts[j]
        assert row['contacts']==[[i,contacts.get(i,2)]for i in ids]
        for i in ids:
            value=-phi if contacts.get(i)==5 else -O if contacts.get(i)==3 else Z
            assert pair(part,N[i])==value==G[30+j][i]
        options=[[a+sign*root*b for a,b in zip(part,v)]for sign in(-1,1)]
        accepted=[]
        for sign,option in zip((-1,1),options):
            assert pair(option,option)==2*O
            if R.sg(pair(option,H))>=0 or R.sg(pair(option,oldE))<=0:continue
            if all(R.relation(pair(option,N[i]),2*O,norms[i])[0]!='forbidden'for i in range(30)):
                accepted.append((sign,option))
        assert len(accepted)==1 and accepted[0][1]==n
        accepted_branches.append(accepted[0][0])
        assert G[30+j][:30]==R.dv(row['pairings'])
        assert norms[30+j]==R.dec(row['norm'])
        assert pair(n,H)==R.dec(row['center_pairing'])and pair(n,oldE)==R.dec(row['old_escape_pairing'])
    assert G[30][31]==-O==R.dec(p['mutual_pairing'])and p['mutual_label']==3
    counts=Counter();new_counts=Counter();saved_new={tuple(x['pair']):x for x in d['new_pairs']}
    for i,j in combinations(range(32),2):
        typ,label,u=R.relation(G[i][j],norms[i],norms[j]);assert typ!='forbidden'
        key=typ+(str(label)if label else'');counts[key]+=1
        if j>=30:
            row=saved_new.pop((i,j))
            assert row['kind']==typ and row['label']==label
            records=[row]
            if i<30:records.append(p['caps'][j-30]['relations'][i])
            for rr in records:
                assert rr['kind']==typ and rr['label']==label
                assert u==(None if rr['invariant']is None else R.dec(rr['invariant']))
            new_counts[key]+=1
    assert not saved_new and dict(counts)==d['pair_counts']
    assert sum(new_counts.values())==61
    assert R.sg(pair(H,H))<0 and all(R.sg(C.dot(H,row))<0 for row in cov)
    assert len(d['facet_points'])==32
    for i,row in enumerate(d['facet_points']):
        P=[h-C.dot(H,cov[i])/norms[i]*a for h,a in zip(H,N[i])]
        assert P==R.dv(row['point'])and pair(P,P)==R.dec(row['norm'])
        assert R.sg(pair(P,P))<0 and R.sg(pair(P,H))<0
        vals=[C.dot(P,cc)for cc in cov]
        assert all(R.sg(g)==(0 if j==i else-1)for j,g in enumerate(vals))
        if i>=30:
            witness=p['caps'][i-30]['facet_witness']
            assert P==R.dv(witness['point'])and pair(P,P)==R.dec(witness['norm'])
            assert vals==R.dv(witness['pairings'])
    cn,old_counts,stage=clip_saved_cone(b,N,Q,H)
    cc={kind:len(rows)for kind,rows in cn.items()}
    assert cc==d['counts']
    assert {kind:sorted(map(list,rows))for kind,rows in cn.items()}==d['complete_census']
    saved={tuple(row['incident']):row for row in d['finite_vertices']}
    assert set(saved)==set(cn['finite_vertices'])
    for ids,v in cn['finite_vertices'].items():
        row=saved[ids];sv=R.dv(row['point'])
        assert len(ids)==9 and R.canon(v)==R.canon(sv)
        assert R.pd([[G[i][j]for j in ids]for i in ids])
        assert pair(sv,sv)==R.dec(row['norm'])and R.sg(pair(sv,sv))<0 and R.sg(pair(sv,H))<0
        vals=[C.dot(sv,row)for row in cov]
        assert vals==R.dv(row['pairings'])and all(R.sg(g)==(0 if i in ids else-1)for i,g in enumerate(vals))
    old={tuple(row['incident'])for row in b['finite_vertices']}
    assert old<=set(saved)and len(saved)-len(old)==d['new_finite_vertices']
    vertices=d['finite_vertices'];edges=[]
    for i,j in combinations(range(len(vertices)),2):
        support=sorted(set(vertices[i]['incident'])&set(vertices[j]['incident']))
        if len(support)==8:
            assert R.pd([[G[a][bb]for bb in support]for a in support])
            edges.append(dict(vertices=[i,j],support=support))
    assert edges==d['finite_edges']
    closed={(i,tuple(row['support']))for row in edges for i in row['vertices']}
    germs=[dict(vertex=i,support=[j for j in row['incident']if j!=release],released_wall=release)
           for i,row in enumerate(vertices)for release in row['incident']
           if(i,tuple(j for j in row['incident']if j!=release))not in closed]
    assert germs==d['open_finite_germs']and 9*len(saved)==2*len(edges)+len(germs)
    oldE=R.dv(b['strict_escape']['point']);vals=[C.dot(oldE,row)for row in cov]
    assert vals==R.dv(d['old_escape_pairings'])
    assert all(R.sg(g)<0 for g in vals[:30])and all(R.sg(g)>0 for g in vals[30:])
    escape=d['strict_escape'];E,V=R.dv(escape['point']),R.dv(escape['hyperideal_ray'])
    ids=tuple(escape['active']);eps=K.convert(Fraction(escape['center_epsilon']))
    assert ids in cn['hyperideal_rays']and R.canon(V)==R.canon(cn['hyperideal_rays'][ids])
    assert R.sg(eps)>0 and E==[v+eps*h for v,h in zip(V,H)]
    vals=[C.dot(E,row)for row in cov]
    assert vals==R.dv(escape['pairings'])and all(R.sg(g)<0 for g in vals)
    assert pair(E,E)==R.dec(escape['norm'])and R.sg(pair(E,E))>0
    assert pair(E,H)==R.dec(escape['center_pairing'])and R.sg(pair(E,H))<0
    paths += [f'scripts/{TAG}.py','scripts/fresh_h9_joint_pair_20261015_twocaps.py',
              'scripts/fresh_h9_joint_pair_20261015_two_extension.py',
              'scripts/fresh_h9_nine_contact_extension_20261015_review.py',
              'scripts/fresh_h9_quartet_polar_extension_20261014_review.py',
              'scripts/fresh_h9_h8_all_high_20261014_review.py',
              'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    joint=[list(ids)for ids in sorted(saved)if 30 in ids and 31 in ids]
    assert joint==d['joint_finite_vertices']==[[0,1,2,5,6,10,29,30,31]]
    joint_indices=[i for i,row in enumerate(vertices)if row['incident']in joint]
    assert not any(i in joint_indices for edge in edges for i in edge['vertices'])
    joint_germs=[row for row in germs if row['vertex']in joint_indices]
    assert len(joint_germs)==9
    polar=R.dv(b['strict_escape']['hyperideal_ray'])
    polar_pairings=[pair(n,polar)for n in N[30:]]
    assert all(R.sg(g)>0 for g in polar_pairings)
    return dict(status='PASS_INDEPENDENT_TWO_CAPS_AND_COMPLETE_32_WALL_CONE',
                sources={p:sha(p)for p in paths},
                arithmetic='Separate SymPy ANP quartic field and192-bit outward rational sign bounds. No producer arithmetic or chamber engine imported.',
                inherited_complete_starting_cone=OLD_REVIEW,
                reconstructed_starting_rays=sum(old_counts.values()),starting_counts=old_counts,
                clipping_stages=stage,counts=cc,complete_census=d['complete_census'],
                cap_norms=enc(norms[-2:]),contact_support=[0,1,2,3,5,6,7,10,29],accepted_branches=accepted_branches,mutual_label=3,new_pair_counts=dict(new_counts),
                pair_counts=dict(counts),all_pair_checks=496,actual_facets=32,facet_sign_checks=1024,
                finite_vertices=len(saved),new_finite_vertices=len(saved)-len(old),
                new_finite_vertex_incidents=[list(ids)for ids in sorted(set(saved)-old)],
                joint_finite_vertex_incidents=joint,joint_vertex_open_germs=joint_germs,
                polar_cap_pairings=enc(polar_pairings),
                finite_vertex_sign_checks=32*len(saved),finite_edges=len(edges),open_finite_germs=len(germs),
                old_escape_cut=True,strict_escape_norm=enc(pair(E,E)),
                strict_escape_center_pairing=enc(pair(E,H)),strict_escape_pairings=enc(vals),
                scope='Complete two-step clipping of the previously independently certified30-wall cone. All894 starting rays are reconstructed exactly; every new crossing edge is detected from its complete face ray set. Both exact quartic recipes, all32facets, complete finite/exterior ray census, joint finite vertices, finite edges, germs and a new strict spacelike escape are checked. This32-wall chamber remains infinite-volume. No arbitrary-cap or finite-volume completion exclusion.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['counts'],data['finite_edges'],data['open_finite_germs'])
