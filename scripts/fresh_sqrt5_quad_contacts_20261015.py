#!/usr/bin/env python3
"""Compactness, H4-payment and complete extra-white tests for Qsqrt5 cores."""
import argparse,importlib.util,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import networkx as nx
import fresh_sqrt5_quad_filters_20261015 as A

ROOT=Path(__file__).resolve().parents[1];BASE='fresh_sqrt5_quad_exhaustion_20261015';TAG='fresh_sqrt5_quad_contacts_20261015'
sp=importlib.util.spec_from_file_location('sqrt5_unions',ROOT/f'scripts/{BASE}.py');E=importlib.util.module_from_spec(sp);sp.loader.exec_module(E)

def edges(row):return dict(zip(combinations(range(row['nodes']),2),row['labels']))
def h4_payments(row):
    g=edges(row);out=[]
    for a,b in product((0,2),(1,3)):
        nodes=row['white']+[a,b];H=nx.Graph();H.add_nodes_from(nodes)
        for i,j in combinations(sorted(nodes),2):
            if g[(i,j)]!=2:H.add_edge(i,j,label=g[(i,j)])
        component=nx.node_connected_component(H,a)
        if b not in component or len(component)!=4:continue
        J=H.subgraph(component)
        if J.degree(a)==J.degree(b)==1 and sorted(dict(J.degree()).values())==[1,1,2,2] and sorted(d['label'] for _,_,d in J.edges(data=True))==[3,3,5]:
            path=nx.shortest_path(J,a,b);assert J[path[0]][path[1]]['label']==5 or J[path[-2]][path[-1]]['label']==5
            out.append(dict(marked=[a,b],entire_H4=list(sorted(component)),path=path))
    return out

def affine_witness(row):
    """Search all induced affine diagrams with labels3/4/6; a found witness
    is independently verified as PSD/corank1 in the Qphi Gram representation.
    Failure to find one is NOT promoted to an all-subdiagram assertion.
    """
    n=row['nodes'];g=edges(row)
    for size in range(3,n+1):
        for ids in combinations(range(n),size):
            adj={i:{} for i in ids};valid=True;edge_count=0
            for i,j in combinations(ids,2):
                m=g[(i,j)]
                if m in (0,5,10):valid=False;break
                if m!=2:adj[i][j]=adj[j][i]=m;edge_count+=1
            if not valid or edge_count<size-1:continue
            seen={ids[0]};stack=[ids[0]]
            while stack:
                for j in adj[stack.pop()]:
                    if j not in seen:seen.add(j);stack.append(j)
            if len(seen)!=size:continue
            affine=False;pivots=[]
            if edge_count==size:
                affine=all(len(adj[i])==2 and all(m==3 for m in adj[i].values()) for i in ids)
            elif edge_count==size-1:
                diag={i:F(1) for i in ids};work={i:dict(a) for i,a in adj.items()}
                while diag:
                    i=next(i for i in diag if len(work[i])<=1);pivot=diag[i];pivots.append(str(pivot))
                    if len(diag)==1:affine=pivot==0;break
                    if pivot<=0:break
                    j,m=next(iter(work[i].items()));diag[j]-={3:F(1,4),4:F(1,2),6:F(3,4)}[m]/pivot
                    del work[j][i],work[i],diag[i]
            if affine:
                subrow=dict(nodes=size,labels=[g[(i,j)] for i,j in combinations(ids,2)])
                colors=A.coloring(subrow);assert colors['passes']
                M=A.scaled_matrix(subrow,colors['colors']);assert A.inertia(M)==[size-1,0,1]
                return dict(nodes=list(ids),labels=subrow['labels'],leaf_pivots=pivots,exact_inertia=[size-1,0,1])
    return None

def profiles(row):
    n=row['nodes'];W=row['white']+[n];g=edges(row);contacts=[2]*n;out=[];nodes=0
    def rec(k):
        nonlocal nodes
        nodes+=1
        if k==n:
            if any(m!=2 for m in contacts):out.append(list(contacts))
            return
        p=(k,n);g[p]=2;contacts[k]=2;rec(k+1)
        for label in (3,4,5,6,10):
            g[p]=label;contacts[k]=label
            if E.no_short(W,g) and E.corners_finite(W,g):rec(k+1)
        del g[p];contacts[k]=2
    rec(0);return out,nodes

def extension(row,contacts):
    n=row['nodes'];g=edges(row)
    for i,m in enumerate(contacts):g[(i,n)]=m
    return dict(nodes=n+1,white=row['white']+[n],labels=[g[p] for p in combinations(range(n+1),2)])

def run():
    base=json.loads((ROOT/f'results/{BASE}.json').read_text())
    filters=json.loads((ROOT/'results/fresh_sqrt5_quad_filters_20261015.json').read_text())
    records=[];counts=Counter();profile_counts=Counter()
    for rec in filters['records']:
        if rec['status'] not in ('PASS_PRINCIPAL_UNION_NECESSARY_TESTS','DOTTED_PARAMETER_UNFILTERED'):continue
        row=base['labelled_union_supports'][rec['union_index']];out=dict(class_index=rec['class_index'],union_index=rec['union_index'])
        paid=h4_payments(row);out['H4_payments']=paid
        dotted=rec['status']=='DOTTED_PARAMETER_UNFILTERED'
        if paid:out['status']='ALREADY_PAID_BY_ENTIRE_H4'
        else:
            aff=affine_witness(row);out['affine_witness']=aff
            if aff:out['status']='REJECT_AFFINE_PRINCIPAL_SUBDIAGRAM'
            else:
                domain,nodes=profiles(row);out['structural_profiles']=len(domain);out['structural_nodes']=nodes;trials=[];accepted=[]
                colors=A.coloring(row) if not dotted else None
                M=A.scaled_matrix(row,colors['colors']) if colors else None
                kernel=A.nullspace([[x.conj() for x in r] for r in M]) if M else []
                out['conjugate_kernel']=[list(map(A.K.enc,v)) for v in kernel]
                for contact in domain:
                    trial=dict(contacts=contact);ext=extension(row,contact);ep=h4_payments(ext)
                    if ep:trial.update(status='BECOMES_PAID_BY_ENTIRE_H4',H4_payments=ep)
                    else:
                        witness=affine_witness(ext)
                        if witness:trial.update(status='REJECT_AFFINE_PRINCIPAL_SUBDIAGRAM',affine_witness=witness)
                        elif dotted:trial['status']='DOTTED_PARAMETER_UNFILTERED';accepted.append(contact)
                        else:
                            ec=A.coloring(ext)
                            if not ec['passes']:trial.update(status='REJECT_CYCLE_FIELD',conflict=ec)
                            else:
                                assert ec['colors'][:-1]==colors['colors']
                                G=A.scaled_matrix(ext,ec['colors']);v=[r[-1].conj() for r in G[:-1]]
                                products=[A.dot(k,v) for k in kernel]
                                if any(products):trial.update(status='REJECT_CONJUGATE_KERNEL',products=[x.enc() for x in products])
                                else:
                                    check=A.conditions(ext,G);trial.update(check)
                                    if not check['physical']['passes']:trial['status']='REJECT_TWO_NEGATIVE_DIRECTIONS'
                                    elif not check['conjugate']['passes']:trial['status']='REJECT_CONJUGATE_NOT_PSD'
                                    else:trial['status']='PASS_EXTRA_WHITE_NECESSARY_TESTS';accepted.append(contact)
                    profile_counts[trial['status']]+=1;trials.append(trial)
                out['trials']=trials;out['accepted']=accepted
                out['dotted']=dotted;out['singular']=rec.get('singular')
                out['status']='UNPAID_DOTTED_CONTACT_DOMAIN' if dotted else 'UNPAID_FINITE_CONTACT_DOMAIN'
        counts[out['status']]+=1;records.append(out)
        if len(records)%25==0:print('processed',len(records),dict(counts),dict(profile_counts),flush=True)
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',
           'scripts/fresh_sqrt5_quad_filters_20261015.py','results/fresh_sqrt5_quad_filters_20261015.json']
    return dict(status='PASS_QSQRT5_H4_AFFINE_AND_ONE_WHITE_NECESSARY_TESTS',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},counts=dict(counts),profile_counts=dict(profile_counts),records=records,
        scope='Complete one-white domains for unpaid cores surviving certified affine witnesses. H4-paid cases need no further quad weight. Dotted parameters and arbitrary collections of external whites remain open. No global dimension bound.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(p.read_text())
    else:p.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['counts'],data['profile_counts'])
