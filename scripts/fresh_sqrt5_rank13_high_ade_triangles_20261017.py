#!/usr/bin/env python3
"""All two-special triangle corners for the joint high-baseline rank13 ADE set."""
import argparse,hashlib,json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_adjacent_discount_witness_20261016 as A

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_rank13_high_ade_triangles_20261017'
E='fresh_sqrt5_rank13_high_ade_endpoints_20261017';K=A.K;PHI=K(F(1,2),F(1,2))

def finite_third(g,port,other,contacts,label):
    new=len(g);h={i:g[i]-{other}for i in g if i!=other};h[new]=set(contacts)-{other}
    for i in h[new]:h[i]=h[i]|{new}
    h[new].add(port);h[port]=h[port]|{new}
    seen={port};queue=[port];limit=4 if label==5 else 2
    for i in queue:
        for j in h[i]:
            if j not in seen:
                seen.add(j);queue.append(j)
                if len(seen)>limit:return None
    degrees={i:len(h[i])for i in seen}
    if sum(degrees.values())!=2*(len(seen)-1)or max(degrees.values())>2:return None
    if min(degrees[port],degrees[new])!=1:return None
    if label==10:assert len(seen)==2;return 'I2(10)'
    return {2:'I2(5)',3:'H3',4:'H4'}[len(seen)]

def run():
    data=json.loads((ROOT/f'results/{E}.json').read_text())
    for p,h in data['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    targets=[];counts=Counter()
    for t in data['targets']:
        n=data['dimension'];g={i:set()for i in range(n)}
        for i,j in t['edges']:g[i].add(j);g[j].add(i)
        gram=[[K(2*(i==j)-int(j in g[i]))for j in range(n)]for i in range(n)]
        inverse=A.inverse(gram);assert all(not x.b for row in inverse for x in row)
        inv=[[x.a for x in row]for row in inverse]
        accepted=[];pairs=set();endpoint_records=[]
        for port_record in t['ports']:
            i=port_record['port'];a=inv[i][i];assert a>0
            for ai,attachment in enumerate(port_record['attachments']):
                contacts=attachment['contacts'];b=-sum((inv[i][j]for j in contacts),F(0))
                c=F(2)-sum((inv[j][k]for j in contacts for k in contacts),F(0))
                assert b<=0
                for label in(5,10):
                    counts['edge_label_tests']+=1
                    if label==10 and b:
                        counts['label10_nonfield_cycle']+=1;continue
                    residual=K(c)+2*b*PHI-a*(1+PHI)if label==5 else K(c)-a*(2+PHI)
                    physical=residual.sign();conjugate=residual.conj().sign()
                    if physical!=-1 or conjugate!=1:
                        counts['edge_signature_rejected']+=1;continue
                    endpoint_records.append(dict(port=i,attachment=ai,label=label,
                        coefficients=[str(a),str(b),str(c)],physical_residual=residual.enc(),
                        conjugate_residual=residual.conj().enc()))
                    for j in range(n):
                        if j==i:continue
                        counts['third_corner_tests']+=1
                        exceptional=finite_third(g,i,j,contacts,label)
                        if exceptional is None:counts['third_corner_nonspherical']+=1;continue
                        pair=tuple(sorted((i,j)));pairs.add(pair)
                        accepted.append(dict(port=i,other=j,attachment=ai,destination=attachment['destination'],
                            label=label,third_exceptional_component=exceptional))
                        counts['accepted_occurrences']+=1;counts[f'accepted_label{label}']+=1
        targets.append(dict(index=t['index'],forest=t['forest'],baseline=t['baseline'],
            corner_count=len(pairs),possible_two_special_triangle_corners=[list(p)for p in sorted(pairs)],
            accepted_occurrences=accepted,arithmetic_endpoint_records=endpoint_records,
            half_charge_corner_upper=str(F(len(pairs),2))))
        print('ADE13_TRIANGLES',t['index'],len(pairs),len(accepted),flush=True)
    paths=[f'scripts/{TAG}.py',f'results/{E}.json',f'scripts/{E}.py',
        'scripts/fresh_sqrt5_adjacent_discount_witness_20261016.py',
        'results/fresh_sqrt5_mixed_d4_d5_migration_20261017.json']
    return dict(status='PASS_COMPLETE_RANK13_HIGH_ADE_TWO_SPECIAL_TRIANGLE_CORNERS',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},counts=dict(counts),targets=targets,
        finite_label_reduction='The squared finite product lies in Q(sqrt5), giving labels2,3,4,5,6,10. Labels2/3 give a rational full edge Gram. For4/6, a fixed3-edge path violates the cycle field; without such a path b=0 and the edge Schur residual c-a*g^2 is rational, so its physical and conjugate signs coincide. Thus only5 and10 survive. Label10 also requires b=0.',
        arithmetic_formula='With old finite Gram inverse J and new contacts g at port i and -1 at S, a=J_ii, b=-sum_(j in S)J_ij, c=2-sum_(j,k in S)J_jk. The new Schur entry is c-2*b*g-a*g^2. Label5 gives g=-phi; label10 gives g^2=2+phi and b=0. Physical negativity and conjugate positivity are checked exactly.',
        finite_third_corner_completeness='After deleting the other released root, only the component containing the framing bond can differ from an induced subdiagram of the finite endpoint forest. Label5 requires a path of2,3,or4 roots with the5-edge at an end; label10 requires an isolated two-root component. These are exactly I2(5),H3,H4,I2(10).',
        all_special_triangle_exclusion='Normalize K-rational normals of totally positive norm a_i by adjoining sqrt(2/a_i). Each nonphysical embedding extends to a real embedding of this extension. If every triangle corner is ADE, the normalized full-star Gram is rational and fixed entrywise; its physical negative direction would then occur in the positive-definite conjugate form, a contradiction.',
        scope='Complete necessary local arithmetic and finite third-corner tests for triangles with two vertices in the selected215-type ADE set. No global ambient realization, all-special quadrilateral classification, or joint dimension bound is claimed. Corner-count upper bounds need not be simultaneously attainable.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    data=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['counts'],flush=True)
