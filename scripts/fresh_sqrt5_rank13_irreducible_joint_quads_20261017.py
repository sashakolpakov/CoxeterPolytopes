#!/usr/bin/env python3
"""Joint-set all-special quadrilaterals based at A13 or D13 vertices."""
import argparse,hashlib,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import sympy as s
import fresh_sqrt5_rank13_high_ade_endpoints_20261017 as E
import fresh_sqrt5_sealed_subcore_closure_20261015_review as J

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_rank13_irreducible_joint_quads_20261017'
INPUT=f'results/{E.TAG}.json';TRI='results/fresh_sqrt5_rank13_high_ade_triangles_20261017.json'
u,v=s.symbols('u v')

def encode(p):return[[list(m),str(c)]for m,c in sorted(s.Poly(p,u,v).terms())]

def factors(row):
    n=row['nodes'];g=2*s.eye(n)
    for(i,j),m in zip(combinations(range(n),2),row['labels']):
        q=u if(i,j)==(0,2)else v if(i,j)==(1,3)else 2-m
        g[i,j]=g[j,i]=q
    white=list(range(4,n));black=list(range(4));W=g.extract(white,white)
    D=g.extract(black,white);R=g.extract(black,black)-D*W.inv()*D.T
    determinant=s.Poly(R.det(method='domain-ge'),u,v);assert not determinant.is_zero
    coefficient,ff=s.factor_list(determinant.as_expr(),u,v);out=[];possible=False
    for f,mult in ff:
        p=s.Poly(f,u,v);a,b=p.degree(u),p.degree(v)
        if not a or not b:
            variable=u if a else v;q=s.Poly(f,variable)
            if q.degree()==1:
                root=-q.nth(0)/q.nth(1);assert root.is_Rational
                status='EXCLUDED_RATIONAL_EDGE';extra=dict(variable=str(variable),root=str(root))
            else:status='UNRESOLVED_UNIVARIATE';extra={};possible=True
        else:status='UNRESOLVED_COUPLED';extra={};possible=True
        out.append(dict(polynomial=encode(p),multiplicity=mult,status=status,**extra))
    assert s.expand(coefficient*s.prod(f**m for f,m in ff)-determinant.as_expr())==0
    # A fixed3-edge path decides whether each opposite product itself is a
    # cyclic product; without it only the square is automatically in K.
    fixed={i:set()for i in range(n)}
    for(i,j),m in zip(combinations(range(n),2),row['labels']):
        if m==3:fixed[i].add(j);fixed[j].add(i)
    def connected(i,j):
        seen={i};queue=[i]
        for x in queue:
            for y in fixed[x]:
                if y not in seen:seen.add(y);queue.append(y)
        return j in seen
    return dict(residual=[[encode(R[i,j])for j in black]for i in black],determinant=encode(determinant),
        coefficient=str(coefficient),factors=out,possible=possible,u_in_K=connected(0,2),v_in_K=connected(1,3))

def run():
    data=json.loads((ROOT/INPUT).read_text());tri=json.loads((ROOT/TRI).read_text())
    for doc in(data,tri):
        for p,h in doc['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    lookup={tuple(x['forest']):x['index']for x in data['targets']}
    records=[];buckets=defaultdict(list);targets=[]
    for t in data['targets']:
        if t['forest'] not in[['D13'],['A13']]:continue
        n=13;g={i:set()for i in range(n)}
        for i,j in t['edges']:g[i].add(j);g[j].add(i)
        maps=[];counts=Counter();possible=set();finite_pairs=set()
        for i,j in combinations(range(n),2):
            for ai,A in enumerate(t['ports'][i]['attachments']):
                for aj,B in enumerate(t['ports'][j]['attachments']):
                    for bond in(0,1):
                        h={x:set(y)for x,y in g.items()};h[n]=set(A['contacts']);h[n+1]=set(B['contacts'])
                        for x in h[n]:h[x].add(n)
                        for x in h[n+1]:h[x].add(n+1)
                        if bond:h[n].add(n+1);h[n+1].add(n)
                        remaining=set(h)-{i,j};names=E.classify({x:h[x]&remaining for x in remaining})
                        counts['tested']+=1;destination=lookup.get(names)
                        if destination is None:counts['fourth_corner_rejected']+=1;continue
                        finite_pairs.add((i,j));order=[i,j,n,n+1]+[x for x in g if x not in(i,j)]
                        labels=[0 if(a,b)in((0,2),(1,3))else 3 if order[b]in h[order[a]]else 2 for a,b in combinations(range(n+2),2)]
                        row=dict(nodes=n+2,white=list(range(4,n+2)),labels=labels);matrix=J.matrix(row);key=J.invariant(matrix)
                        found=None
                        for index in buckets[key]:
                            mp=J.isomorphism(matrix,J.matrix(records[index]['row']))
                            if mp is not None:J.check_map(matrix,J.matrix(records[index]['row']),mp);found=(index,mp);break
                        if found is None:
                            index=len(records);records.append(dict(index=index,row=row,symbolic=factors(row)))
                            buckets[key].append(index);found=(index,list(range(n+2)))
                        index,mp=found
                        if records[index]['symbolic']['possible']:possible.add((i,j));counts['unresolved_literal_patterns']+=1
                        else:counts['rational_branch_excluded_patterns']+=1
                        maps.append(dict(ports=[i,j],first_attachment=ai,second_attachment=aj,new_new_label=3 if bond else 2,
                            adjacent_destination_types=[A['destination'],B['destination']],opposite_destination_type=destination,
                            record=index,map=mp,swap_free_products=tuple(sorted((mp[0],mp[2])))==(1,3)))
            print('JOINT_QUAD_PAIR',t['forest'][0],i,j,'records',len(records),'possible',len(possible),flush=True)
        triangles={tuple(p)for p in tri['targets'][t['index']]['possible_two_special_triangle_corners']}
        alpha=sum((max(F(1,2)if p in triangles else F(0),F(1,4)if p in possible else F(0))for p in combinations(range(n),2)),F(0))
        targets.append(dict(index=t['index'],forest=t['forest'],baseline=t['baseline'],counts=dict(counts),maps=maps,
            finite_corner_pairs=[list(p)for p in sorted(finite_pairs)],possible_quad_corners=[list(p)for p in sorted(possible)],
            possible_quad_corner_count=len(possible),two_special_triangle_corners=[list(p)for p in sorted(triangles)],
            combined_half_triangle_quarter_quad_bound=str(alpha),threshold='12',within_threshold=alpha<=12))
    paths=[f'scripts/{TAG}.py',INPUT,TRI,f'scripts/{E.TAG}.py',f'scripts/{J.TAG}.py',
        'results/fresh_sqrt5_mixed_d4_d5_quad_factor_gate_20261017.json']
    return dict(status='PASS_COMPLETE_JOINT_HIGH_ADE_QUAD_PATTERNS_AT_RANK13_IRREDUCIBLES',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},records=records,targets=targets,
        scope='Complete all-selected finite quadrilateral patterns and rational-factor exclusions based at A13 and D13, while all other corners range over the215-type joint set. Coupled and nonlinear univariate branches are retained without a parameter cutoff. Bounds are local per-vertex charges only; the other213 selected types and every unselected vertex still require control. No dimension or global joint incidence bound is asserted.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();data=run()
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],[(t['forest'],t['possible_quad_corner_count'],t['combined_half_triangle_quarter_quad_bound'])for t in data['targets']],flush=True)
