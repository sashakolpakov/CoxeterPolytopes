#!/usr/bin/env python3
"""Independent full induced-subset replay of twelve local core transfers.

Standard-library only; no producer, geometric-engine or algebra imports.
The two source obstruction theorems remain explicit dependencies.
"""
import argparse,json
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_star_core_transfers_20261014_review'
PRIMARY='results/fresh_h9_star_core_transfers_20261014.json'
PRODUCER='scripts/fresh_h9_star_core_transfers_20261014.py'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
THEOREMS=['results/fresh_h9_star_mixed_singles_20261014.json','results/fresh_h9_star_opposite_core_20261014.json']
REVIEWS=['results/fresh_h9_star_mixed_singles_20261014_review.json','results/fresh_h9_star_opposite_core_20261014_review.json']
BLOCKS=[(3,4),(6,7),(9,10),(12,13)]
P=[2,5,8,11];C=[14,15,16,17]
CORE=sorted([x for b in BLOCKS for x in b]+[11,14,17])
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def bitset(v):return sum(1<<x for x in v)
def elements(mask):return [i for i in range(mask.bit_length())if mask>>i&1]
def run():
    seen={}
    def read(p):
        if p in seen:return seen[p]
        d=json.loads((ROOT/p).read_text());seen[p]=d
        for f,h in d.get('sources',{}).items():
            assert sha(f)==h,f
            if f.endswith('.json'):read(f)
        return d
    d=read(PRIMARY);g=read(GEOM)['representatives'][0]['geometry']
    th=[read(p)for p in THEOREMS]
    reviews=[read(p)for p in REVIEWS]
    for p,r in zip(THEOREMS,reviews):
        assert r['status'].startswith('PASS_INDEPENDENT_')and r['sources'][p]==sha(p)
    assert d['status']=='PASS_TWELVE_LITERAL_CORE_TRANSFERS'
    assert [t['status']for t in th]==['PASS_CORE_SHARED_SINGLE_ALL_INTERNAL_LABELS','PASS_CORE_OPPOSITE_SINGLE_ALL_INTERNAL_LABELS']
    vertices=[bitset(v['incident'])for v in g['vertices']]
    allmissing={bitset(m)for m in g['minimal_nonfaces']}
    assert all(v.bit_count()==9 for v in vertices)
    def incidence(labels):
        # Unlike the producer's maximal-intersection calculation, test all
        # subsets directly against the full saved actual-vertex list.
        masks=[bitset(labels[i]for i in range(11)if b>>i&1)for b in range(2048)]
        faces=[any(m&v==m for v in vertices)for m in masks]
        maximal=[masks[b]for b in range(2048)if faces[b]and all(not faces[b|(1<<i)]for i in range(11)if not b>>i&1)]
        minimal=[masks[b]for b in range(1,2048)if not faces[b]and all(faces[b^(1<<i)]for i in range(11)if b>>i&1)]
        return masks,faces,sorted(maximal),sorted(minimal)
    sm,sf,smax,smin=incidence(CORE)
    assert sum(sf)==1456 and len(smax)==9 and len(smin)==3
    assert set(smin)=={bitset([3,4,14]),bitset([11,12,13]),bitset([12,13,17])}
    assert all(sf[(1<<i)|(1<<j)]for i in range(11)for j in range(i))
    assert set(map(bitset,d['source_maximal_faces']))==set(smax)
    assert set(map(bitset,d['source_induced_minimal_nonfaces']))==set(smin)
    for t in th:
        assert t['rank_core']==CORE
        assert set(map(bitset,t['complete_maximal_core_faces']))==set(smax)
        assert len(t['families'])==2
    # Independently check the literal incidence premises that the source
    # proofs advertise as core-only; their algebra remains inherited.
    source_witness_counts=[]
    for t in th:
        witness_rows=[]
        for k in ('core_positive_kernel_guards','core_positive_quadratic_guards','projection_budgets'):
            witness_rows+=t.get(k,[])
        for k in ('high_last_third_right','mutual_upper_six'):
            if k in t:witness_rows.append(t[k])
        witness_rows += [r['witness']for r in t['new_core_fallbacks']]
        for w in witness_rows:
            support=bitset(w['support']);assert support&bitset(CORE)==support
            assert support&vertices[w['actual_vertex']]==support
            if 'component_order'in w:assert set(w['component_order'])<=set(CORE)
        for r in t['new_core_fallbacks']:
            assert all(a in CORE and b in CORE for a,b,_ in r['relaxed_core_labels'])
        source_witness_counts.append(len(witness_rows))
    assert source_witness_counts==[479,3644]
    source_families={(path,f['id']):f for path,t in zip(THEOREMS,th)for f in t['families']}
    target_families={f['id']:f for f in d['families']}
    assert len(target_families)==48
    placements=[];pairs=set();familykeys=set();used=set();witnesses=0
    for row in d['placements']:
        anchor,owner=row['anchor_block'],row['owner_block']
        assert anchor!=owner and 0<=anchor<4 and 0<=owner<4
        assert (anchor,owner)not in pairs;pairs.add((anchor,owner))
        mp=dict(row['core_map']);assert set(mp)==set(CORE)and len(set(mp.values()))==11
        assert mp[11]==P[anchor]and mp[17]==C[anchor]
        assert mp[14]==(C[owner]if owner<anchor else P[owner])
        assert row['anchors']==[P[anchor],C[anchor]]and row['third_leaf']==mp[14]
        blockmap=[owner]+[i for i in range(4)if i not in (anchor,owner)]+[anchor]
        assert row['block_map']==blockmap
        for i,j in enumerate(blockmap):assert [mp[x]for x in BLOCKS[i]]==list(BLOCKS[j])
        ordered=[mp[x]for x in CORE]
        tm,tf,tmax,tmin=incidence(ordered)
        assert tf==sf # all2048 face/nonface decisions preserved
        target=bitset(ordered)
        assert sorted(ordered)==row['target_core']
        assert set(tmin)=={m for m in allmissing if m&target==m}
        assert set(tmin)==set(map(bitset,row['induced_minimal_nonfaces']))
        ws=row['maximal_face_witnesses'];assert len(ws)==9
        assert set(tmax)=={bitset(x['support'])for x in ws}
        for x in ws:
            mask=bitset(x['support']);v=vertices[x['actual_vertex']]
            assert mask&v==mask;witnesses+=1
        assert len(row['family_ids'])==4
        for fid in row['family_ids']:
            assert fid not in used;used.add(fid);f=target_families[fid]
            src=source_families[(f['source_theorem'],f['source_family'])]
            for kind,n in [('fixed',6),('bounds',4)]:
                assert len(src[kind])==len(f[kind])==n
                expected=sorted([min(mp[a],mp[b]),max(mp[a],mp[b]),label]for a,b,label in src[kind])
                assert f[kind]==expected
                assert all((1<<a|1<<b)&target==(1<<a|1<<b)for a,b,_ in f[kind])
            fk=tuple(tuple(tuple(x)for x in f[k])for k in ('fixed','bounds'))
            assert fk not in familykeys;familykeys.add(fk)
        placements.append(dict(anchor=anchor,owner=owner,all_induced_subsets=2048,
            face_count=sum(tf),maximal_faces=[elements(m)for m in tmax],
            minimal_nonfaces=[elements(m)for m in tmin],literal_vertex_witnesses=9,
            transported_family_ids=row['family_ids']))
    assert pairs=={(i,j)for i in range(4)for j in range(4)if i!=j}
    assert used==set(target_families)and len(familykeys)==48 and witnesses==108
    assert d['placement_count']==12 and d['predicate_count']==48
    sources={p:sha(p)for p in [f'scripts/{TAG}.py',PRODUCER,PRIMARY,GEOM,*THEOREMS,*REVIEWS]}
    return dict(status='PASS_INDEPENDENT_TWELVE_INDUCED_CORE_TRANSFERS',sources=sources,
        arithmetic='Standard-library integer bit sets; no producer imports.',
        source_core=CORE,source_induced_subsets=2048,source_face_count=sum(sf),
        placements=placements,total_target_subsets=24576,total_face_decisions=26624,
        transported_predicates=48,literal_maximal_face_witnesses=108,
        inherited_theorem_literal_core_witnesses=source_witness_counts,
        proof_review='Every face/nonface in each eleven-node induced core is preserved, not merely the displayed witnesses. The source theorems use only these core faces, four internal bounds, six fixed contacts and Lorentz rank/index. Relabeling therefore transports each theorem. No whole-nerve automorphism is required.',
        dependencies='This replay does not reprove the two obstruction theorems or their determinant banks. Both existing independent source-theorem reviews are explicitly hash-bound.',
        scope='Forty-eight conditional label-family exclusions on star0; not an exclusion of the complete incidence or of any unlisted family.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS independent saved twelve-core transfer replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print('PASS independent twelve-core transfer generation')
