#!/usr/bin/env python3
"""Final source-bound completeness, old-support containment and quotient audit.

No full DFS repeated. Inherits the separately frozen independent protocol/two-
block audit; independently checks final accounting, every retained literal
corner/provenance, all5440old supports with exact Fraction principal inertia,
and every explicit quotient bijection. Opens only a frozen read-only SQLite.
"""
import argparse,json,sqlite3,sys
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from math import comb,perm
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import fresh_sqrt5_all_quad_unions_20261015_review as R
Q=R.Q
BASE='fresh_sqrt5_all_quad_unions_20261015';TAG=BASE+'_review_complete'

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    data=json.loads((ROOT/f'results/{BASE}.json').read_text());state=data['state'];counts=state['counts']
    assert data['status']=='COMPLETE_ALL_QUAD_UNIONS'
    assert [state[k]for k in('left','right','mapping','stack','complete')]==[153,153,0,None,True]
    bank=R.complete_pieces();white=[r['order']-2 for r in bank]
    expected=sum(sum(comb(b,k)*perm(a,k)for k in range(min(a,b)+1))for i,a in enumerate(white)for b in white[i:])
    assert expected==counts['base_identifications']==127875 and counts['completed_left_pieces']==len(bank)==153
    assert counts['dfs_nodes']==counts['base_corner_and_field_pass']+counts['label_trials']-counts['corner_PD_prunes']-counts['finite_cycle_field_prunes']==977040
    last=0
    for batch in state['batches']:
        assert batch['initial_nodes']==last and batch['final_nodes']>=last
        last=batch['final_nodes']
    assert last==counts['dfs_nodes']and state['batches'][-1]['reason']=='COMPLETE'
    path=ROOT/data['database'];assert digest(data['database'])==data['database_sha256']
    con=sqlite3.connect(f'file:{path}?mode=ro',uri=True)
    assert con.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
    saved=json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0]);assert saved==state
    actual=dict(con.execute('SELECT key,status FROM examined'));status_counts=Counter(actual.values())
    assert dict(status_counts)==data['distinct_status_counts']
    assert set(status_counts)=={'PASS_DOTTED_PARAMETER_UNFILTERED','PASS_FINITE_PRINCIPAL_NECESSARY_TESTS','REJECT_CONJUGATE_NOT_PSD','REJECT_ORTHOGONAL_LANNER_PAIR','REJECT_TWO_NEGATIVE_DIRECTIONS'}
    assert all(counts[k]==v for k,v in status_counts.items())
    assert len(actual)+counts['duplicate_complete_union']==counts['labelled_completions']==100796
    rows=data['labelled_union_supports'];database=list(con.execute("SELECT key,status,payload FROM examined WHERE status LIKE 'PASS_%' ORDER BY key"));assert len(rows)==len(database)==38755
    corner_checks=0
    for row,(key,status,payload)in zip(rows,database):
        assert row==json.loads(payload);n=row['nodes'];g=Q.edges(row)
        assert key==R.key(n,g)and row['white']==list(range(4,n))and R.corners(n,g)and R.cycle_ok(n,g)
        assert (0 in row['labels'])==(status=='PASS_DOTTED_PARAMETER_UNFILTERED')
        L,U=row['left_nodes'],row['right_nodes'];assert L[:2]==[0,2]and U[:2]==[1,3]
        assert sorted(set(L+U))==list(range(n))
        for piece,mp in ((bank[row['left_piece']],L),(bank[row['right_piece']],U)):
            assert len(mp)==len(set(mp))==piece['order']
            assert all(g[Q.pair(mp[i],mp[j])]==m for(i,j),m in piece['edges'].items())
        assert set(L)&set(U)or any(g[Q.pair(a,b)]!=2 for a,b in product(L,U))
        corner_checks+=4
    old=json.loads((ROOT/'results/fresh_sqrt5_quad_exhaustion_20261015.json').read_text());assert len(old['labelled_union_supports'])==5440
    old_status=Counter();old_records=[]
    for index,row in enumerate(old['labelled_union_supports']):
        n=row['nodes'];g=Q.edges(row);assert R.corners(n,g)
        if not R.cycle_ok(n,g):
            assert R.key(n,g)not in actual;status='EARLY_FIXED_FINITE_CYCLE_REJECTION'
        else:
            status=R.status(n,g,row['left_nodes'],row['right_nodes']);assert actual[R.key(n,g)]==status
        old_status[status]+=1;old_records.append(dict(old_index=index,status=status))
    con.close()
    classes=json.loads((ROOT/f'results/{BASE}_classes.json').read_text());assert classes['literal_rows']==38755 and classes['marked_classes']==8298
    seen=set();pairchecks=0;types=Counter();opposites={frozenset((0,2)),frozenset((1,3))}
    for record in classes['marked_isomorphism_classes']:
        ri=record['representative'];target=rows[ri];n=target['nodes'];tg=Q.edges(target);member_ids=[]
        types['dotted'if 0 in target['labels']else'finite']+=1
        for member in record['members']:
            i=member['index'];assert i not in seen;seen.add(i);member_ids.append(i);source=rows[i];mp=member['to_representative']
            assert source['nodes']==n and sorted(mp)==list(range(n))and all((j<4)==(mp[j]<4)for j in range(n))
            assert {frozenset((mp[0],mp[2])),frozenset((mp[1],mp[3]))}==opposites
            for(a,b),label in Q.edges(source).items():assert tg[Q.pair(mp[a],mp[b])]==label;pairchecks+=1
        assert ri in member_ids
    assert seen==set(range(38755))and pairchecks==classes['literal_pair_checks']==1015927
    assert dict(types)==classes['classes_by_parameter_type']=={'finite':5904,'dotted':2394}
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',data['database'],f'scripts/{BASE}_classes.py',f'results/{BASE}_classes.json',
        f'scripts/{BASE}_review.py',f'results/{BASE}_review.json','scripts/fresh_sqrt5_quad_contacts_20261015_review.py','results/fresh_sqrt5_quad_contacts_20261015_review.json',
        'results/fresh_sqrt5_quad_exhaustion_20261015.json','results/fresh_sqrt5_quad_filters_20261015_review.json']
    Q.manifest_check(paths)
    return dict(status='PASS_COMPLETE_ENUMERATION_ACCOUNTING_ALL_5440_PRIOR_SUPPORTS_AND_MARKED_QUOTIENT_REVIEW',sources={p:digest(p)for p in paths},
        complete_left_blocks=153,base_identifications=expected,dfs_nodes=977040,completed_distinct_graphs=len(actual),terminal_status_counts=dict(status_counts),
        retained_rows=38755,independent_retained_corner_checks=corner_checks,prior_support_count=5440,prior_status_counts=dict(old_status),prior_support_records=old_records,
        quotient=dict(classes=8298,types=dict(types),member_rows=38755,literal_pair_checks=pairchecks),
        scope='Full terminal-accounting/old-census-containment audit inherits frozen independent traversal source/protocol proof. Every retained corner and piece provenance rechecked; old5440 exact principal statuses independently recomputed. Full new finite terminal inertia replay is producer-owned, not duplicated here. Quotient bijections independently checked; canonical completeness source-reviewed. No external-white/rank-completion/realization/global dimension claim.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();data=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert data==json.loads(out.read_text());print('PASS saved completion review')
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['prior_status_counts'],data['quotient'])
