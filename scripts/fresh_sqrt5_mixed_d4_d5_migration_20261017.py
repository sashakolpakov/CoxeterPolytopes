#!/usr/bin/env python3
"""Complete finite same-type endpoint attachments and two-special triangles.
No framing-distance cutoff, quadrilateral-completeness or dimension claim.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_triangle_extension_complete_20261016 as T
Q=T.Q;ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_mixed_d4_d5_migration_20261017'
TARGETS=[('3D4',['D4']*3),('D5_2D4',['D5','D4','D4']),('2D5_D4',['D5','D5','D4']),('3D5_A2',['D5']*3+['A2'])]
ALPHABET=(2,3,4,5,6,10)

def forest(names):
    g={};components=[]
    for name in names:
        n=int(name[1:]);start=len(g);ids=list(range(start,start+n))
        for i in ids:g[i]=set()
        es=[(0,1)]if name=='A2'else[(0,1),(0,2),(0,3)]+([(3,4)]if name=='D5'else[])
        for a,b in es:g[start+a].add(start+b);g[start+b].add(start+a)
        components.append(dict(type=name,nodes=ids))
    return g,components

def classify(g):
    unseen=set(g);out=[]
    while unseen:
        found={min(unseen)};stack=list(found)
        while stack:
            i=stack.pop()
            for j in g[i]:
                if j not in found:found.add(j);stack.append(j)
        unseen-=found;degrees=sorted(len(g[i])for i in found)
        name={(1,1):'A2',(1,1,1,3):'D4',(1,1,1,2,3):'D5'}.get(tuple(degrees))
        if name is None:return None
        assert sum(degrees)==2*(len(found)-1)
        out.append(dict(type=name,nodes=sorted(found)))
    return out

def endpoints(names):
    g,comps=forest(names);n=len(g);target=Counter(names);out=[]
    for port in range(n):
        old=next(c for c in comps if port in c['nodes']);remaining=[i for i in g if i!=port]
        matches=[]
        for size in range(1,4):
            for contacts in combinations(remaining,size):
                h={i:g[i]-{port}for i in remaining};h[n]=set(contacts)
                for i in contacts:h[i]=h[i]|{n}
                result=classify(h)
                if result is None or Counter(c['type']for c in result)!=target:continue
                connected=next(c for c in result if n in c['nodes'])
                matches.append(dict(contacts=list(contacts),copy=set(contacts)==g[port],
                    destination_component=connected,
                    old_component_indices=[i for i,c in enumerate(comps)if set(c['nodes'])&set(connected['nodes'])]))
        assert matches
        if old['type']=='D4':assert len(matches)==1 and matches[0]['copy']
        out.append(dict(port=port,old_type=old['type'],old_local_port=old['nodes'].index(port),attachments=matches))
    return g,comps,out

def row_from(g,port,other,contacts,label):
    n=len(g);order=[port,other,n]+[i for i in g if i not in(port,other)];mp={i:j for j,i in enumerate(order)}
    es=[(mp[i],mp[j],3)for i in g for j in g[i]if i<j]
    es +=[(mp[i],mp[n],3)for i in contacts]
    if label!=2:es.append((mp[port],mp[n],label))
    return T.from_edges(n+1,es),order

def triangle_arithmetic(row):
    col=Q.coloring(row)
    if not col['passes']:return dict(status='REJECT_CYCLE_FIELD')
    m=Q.scaled_matrix(row,col['colors']);ids=list(range(1,row['nodes']))
    s,_=Q.schur(m,ids);sign=s[0][0].sign()
    if sign!=-1:return dict(status='REJECT_NONLORENTZ_TRIANGLE',physical_last_sign=sign)
    # The frozen arithmetic routine supplies independent conjugate pivots.
    return T.arithmetic(row)

def one(name,names):
    g,components,ports=endpoints(names);n=len(g);counts=Counter();accepted=[];pairs=set();rows=[];buckets={}
    for rec in ports:
        i=rec['port']
        for ai,attach in enumerate(rec['attachments']):
            for j in range(n):
                if i==j:continue
                for label in ALPHABET:
                    row,order=row_from(g,i,j,attach['contacts'],label);info=T.corner_info(row)
                    if info is None:counts['REJECT_NONSPHERICAL_THIRD_CORNER']+=1;continue
                    check=triangle_arithmetic(row);counts[check['status']]+=1
                    if check['status']!='ACCEPT':continue
                    corner=info['corners'][1]
                    types=sorted(c['type']for c in corner['components'])
                    assert types!=sorted(names),'An all-special triangle would have a rational indefinite Gram'
                    pair=tuple(sorted((i,j)));pairs.add(pair)
                    key=T.bucket(row);found=None
                    for idx in buckets.get(key,[]):
                        mp=T.isomorphism(row,rows[idx]['row'])
                        if mp is not None:found=(idx,mp);break
                    if found is None:
                        idx=len(rows);rows.append(dict(index=idx,row=row,arithmetic=check,third_corner_types=types))
                        buckets.setdefault(key,[]).append(idx);found=(idx,list(range(n+1)))
                    accepted.append(dict(port=i,other=j,attachment=ai,label=label,third_corner_types=types,
                        full_original_order=order,target=found[0],map=found[1],migrating=not attach['copy']))
        print('MIXED_PORT',name,i,'accepted',len(accepted),'classes',len(rows),flush=True)
    return dict(name=name,dimension=n,forest=names,components=components,ports=ports,
        endpoint_attachments=sum(len(p['attachments'])for p in ports),endpoint_copies=sum(a['copy']for p in ports for a in p['attachments']),
        triangle_tests=dict(counts),accepted_triangle_occurrences=accepted,triangle_classes=rows,
        possible_two_special_triangle_corners=[list(p)for p in sorted(pairs)],
        corner_count=len(pairs),corner_rank_pair_counts=dict(Counter('+'.join(sorted(next(c['type']for c in components if i in c['nodes'])for i in p))for p in pairs)))

def run():
    results=[one(name,names)for name,names in TARGETS]
    calibration=results[0];assert calibration['endpoint_attachments']==12 and calibration['corner_count']==9
    paths=[f'scripts/{TAG}.py',f'scripts/{T.TAG}.py',f'results/{T.TAG}_rank22.json',
           'scripts/fresh_sqrt5_quad_filters_20261015.py','results/fresh_sqrt5_quad_filters_20261015_review.json']
    return dict(status='PASS_COMPLETE_SELECTED_MIXED_D_FINITE_TRANSITIONS_AND_TRIANGLES',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},targets=results,finite_labels=list(ALPHABET),
        endpoint_completeness='Every destination vertex is the specified simply-laced forest. The new root therefore has only label3 contacts and degree at most3. All subsets of the entire retained root set of sizes1,2,3 are tested; a disconnected new root is impossible because none of the target forests contains A1.',
        finite_label_completeness='The two-cycle forces4*cos(pi/m)^2 into Q(sqrt5). The real cyclotomic degree restriction leaves2,3,4,5,6,8,10,12;8 and12 require sqrt2 or sqrt3 and are excluded. The remaining exact labels are tested, with cycle-field and conjugate checks. A framing bridge need not itself lie in the ground field.',
        no_all_special_triangle='All three finite corners would force every label in the full triangle star to be2 or3. Its norm-two Gram is rational; its physical negative direction then survives the nontrivial embedding, contradicting conjugate positivity.',
        scope='Complete finite endpoint-attachment and two-special-triangle possibilities for the four literal target forests only. Accepted rows satisfy the exact necessary physical/conjugate arithmetic and all finite corners; their compact local face geometry is not independently reconstructed here. No all-parameter quadrilateral closure, vertex-frequency inequality, dimension bound or compact ambient realization is claimed.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
