#!/usr/bin/env python3
"""All known literal representatives of the retained triangle rows, with exact graph maps."""
import argparse,hashlib,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_discriminant_gate_20261017 as G
import fresh_sqrt5_rank13_all_selected_nonade_separation_20261017 as E
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_triangle_all_representatives_20261017'
FULL='fresh_sqrt5_rank13_joint_corner_all_nonade_repair_20261017'
FLAGS='fresh_sqrt5_rank13_triangle_simplex_flags_20261017'
K,Z,O=G.K,G.Z,G.O

def graph(edges):
    g={i:{} for i in range(14)}
    for i,j,m in edges:g[i][j]=g[j][i]=m
    return g

def refinement(g):
    colors={i:0 for i in g}
    for _ in range(14):
        codes={i:(colors[i],tuple(sorted((m,colors[j]) for j,m in row.items()))) for i,row in g.items()}
        ids={v:i for i,v in enumerate(sorted(set(codes.values())))};colors={i:ids[codes[i]] for i in g}
    key=(tuple(sorted(Counter(colors.values()).items())),tuple(sorted((min(colors[i],colors[j]),max(colors[i],colors[j]),m)
        for i,row in g.items() for j,m in row.items() if i<j)))
    return colors,key

def isomorphism(g,h,cg,ch):
    mapping={};used=set()
    def visit():
        if len(mapping)==14:return dict(mapping)
        choices=[]
        for i in g:
            if i in mapping:continue
            candidates=[j for j in h if j not in used and cg[i]==ch[j] and
                all(g[i].get(k,2)==h[j].get(v,2) for k,v in mapping.items())]
            if not candidates:return None
            choices.append((len(candidates),-len(g[i]),i,candidates))
        _,_,i,candidates=min(choices)
        for j in candidates:
            mapping[i]=j;used.add(j);out=visit()
            if out is not None:return out
            used.remove(j);del mapping[i]
        return None
    answer=visit()
    if answer is not None:
        assert len(set(answer.values()))==14
        assert all(g[i].get(j,2)==h[answer[i]].get(answer[j],2) for i,j in combinations(range(14),2))
    return answer

def run():
    gate=B.read(G.TAG);flags=B.read(FLAGS);active=B.read(G.INPUT);full=B.read(FULL)
    flagmap={r['face']:r for r in flags['records']};catalogues={};finite={B.C.code(g):name for name,g in B.C.models() if len(g)<=13}
    classes=[dict(c) for c in gate['classes']];component_cache={};stars=[];stargraphs=[];colors=[];buckets=defaultdict(list)
    literal_cache={};references=[];occurrences=[];counts=Counter()
    def read(tag):
        if tag not in catalogues:catalogues[tag]=B.read(tag)
        return catalogues[tag]
    def literal(p):
        if p['kind']=='old_triangle':
            w=read(p['source'])['triangle_patterns'][p['index']]['witness']
            t=read(G.ADE_ENDPOINT)['targets'][w['target']]
            o=read(G.ADE_TRIANGLE)['targets'][w['target']]['accepted_occurrences'][w['occurrence']]
            g=graph([[i,j,3] for i,j in t['edges']]);port=o['port']
            for i in t['ports'][port]['attachments'][o['attachment']]['contacts']:g[i][13]=g[13][i]=3
            if o['label']!=2:g[port][13]=g[13][port]=o['label']
            return g,[port,o['other'],13]
        t=read(p['source'])['targets'][p['target']]
        return E.literal(t,t['triangle_occurrences'][p['occurrence']])
    def spherical(g):
        out=[]
        for removed in range(14):
            end={i:{j:m for j,m in row.items() if j!=removed} for i,row in g.items() if i!=removed};names=[]
            for p in B.C.pieces(end):
                if sum(map(len,p.values()))!=2*len(p)-2:break
                code=B.C.code(p)
                if code not in finite:break
                names.append(finite[code])
            else:out.append(dict(omitted=removed,forest=sorted(names)))
        return out
    def determinant(g):
        scale=[O]*14;special=[(i,j,m) for i,row in g.items() for j,m in row.items() if i<j and m not in(3,5)]
        assert len(special)<=1
        if special:
            i,j,m=special[0];assert j==13 and m in(4,6,10)
            cut={a:{b:z for b,z in row.items() if {a,b}!={i,j}} for a,row in g.items()}
            side=next(p for p in B.C.pieces(cut) if 13 in p);assert i not in side
            for a in side:scale[a]={4:K(2),6:K(3),10:K(2)+G.PHI}[m]
        old={i:{j:m for j,m in row.items() if j!=13} for i,row in g.items() if i!=13}
        vector=[]
        for i in range(13):
            m=g[i].get(13,2)
            if m==2:v=Z
            elif m in(3,5):assert scale[i]==scale[13];v=-scale[i]*{3:O,5:G.PHI}[m]
            else:v=-scale[13]
            vector.append(v)
        det=O;schur=2*scale[13]
        for piece in B.C.pieces(old):
            ids=sorted(piece);key=(len(ids),tuple((ids.index(i),ids.index(j),m) for i in ids for j,m in sorted(piece[i].items()) if i<j))
            if key not in component_cache:
                M=[[K(2) if i==j else -{3:O,5:G.PHI}[piece[i][j]] if j in piece[i] else Z for j in ids] for i in ids]
                component_cache[key]=(G.det_positive(M),G.A.inverse(M))
            dd,inv=component_cache[key];s=scale[ids[0]];assert all(scale[i]==s for i in ids)
            for _ in ids:det*=s
            det*=dd;v=[vector[i] for i in ids];schur-=G.A.dot(v,G.A.mv(inv,v))/s
        assert schur.sign()<0 and schur.conj().sign()>0
        return det*schur
    def classify(g):
        sizes=sorted(len(p) for p in B.C.pieces(g));ss=spherical(g);assert len(ss)>=3
        d=dict(component_sizes=sizes,fixed_discriminant=all(n%2==0 for n in sizes),spherical_deletions=ss,
            simplex_dimension=len(ss)-1)
        if d['fixed_discriminant']:
            det=determinant(g);ci=None;root=None
            for i,c in enumerate(classes):
                root=G.ksqrt(det/K(*c['representative']))
                if root is not None:ci=i;break
            if ci is None:ci=len(classes);classes.append(dict(index=ci,representative=det.enc(),first_new_star=len(stars)));root=O
            d.update(canonical_determinant=det.enc(),discriminant_class=ci,square_ratio_root=root.enc())
        return d
    def register(g,known=None):
        exact=tuple((i,j,m) for i,row in g.items() for j,m in sorted(row.items()) if i<j)
        if exact in literal_cache:return literal_cache[exact]
        cg,key=refinement(g)
        for si in buckets[key]:
            mp=isomorphism(g,stargraphs[si],cg,colors[si])
            if mp is not None:
                answer=(si,[mp[i] for i in range(14)]);literal_cache[exact]=answer;return answer
        si=len(stars);record=dict(index=si,edges=[list(x) for x in exact])
        if known is None:
            record.update(classify(g));record['new_beyond_selected_representatives']=True
        else:
            r,f=known;record.update(component_sizes=r['component_sizes'],fixed_discriminant=r['fixed_discriminant'],
                spherical_deletions=[dict(omitted=x['omitted'],forest=x['forest']) for x in f['spherical_deletions']],
                simplex_dimension=f['simplex_dimension'],new_beyond_selected_representatives=False,reference_face=r['face'])
            for k in ('canonical_determinant','discriminant_class','square_ratio_root'):
                if k in r:record[k]=r[k]
        stars.append(record);stargraphs.append(g);colors.append(cg);buckets[key].append(si)
        answer=(si,list(range(14)));literal_cache[exact]=answer;return answer
    for r in gate['records']:
        si,mp=register(graph(r['edges']),(r,flagmap[r['face']]))
        s=stars[si];assert s['fixed_discriminant']==r['fixed_discriminant']
        if r['fixed_discriminant']:assert s['discriminant_class']==r['discriminant_class']
        assert sorted(mp[x['omitted']] for x in flagmap[r['face']]['spherical_deletions'])==sorted(x['omitted'] for x in s['spherical_deletions'])
        references.append(dict(face=r['face'],star=si,isomorphism=mp))
    original_stars=len(stars)
    for ri,r in enumerate(gate['records']):
        source_face=full['faces'][active['faces'][r['face']]['original_face']]
        for pi,p in enumerate(source_face['provenance']):
            g,black=literal(p);si,mp=register(g)
            assert all(g[i].get(j,2)==stargraphs[si][mp[i]].get(mp[j],2) for i,j in combinations(range(14),2))
            assert set(mp[i] for i in black)<={x['omitted'] for x in stars[si]['spherical_deletions']}
            occurrences.append(dict(face=r['face'],original_face=source_face['index'],provenance_index=pi,provenance=p,
                star=si,isomorphism=mp,black=black))
            counts[p['kind']]+=1
        if ri%250==0:print('ALL_REPRESENTATIVES',ri,len(occurrences),len(stars),len(classes),flush=True)
    assert len(occurrences)==46032
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{G.TAG}.py',f'scripts/{E.TAG}.py',
        f'results/{G.TAG}.json',f'results/{FLAGS}.json',f'results/{G.INPUT}.json',f'results/{FULL}.json']+[f'results/{t}.json' for t in sorted(catalogues)]
    return dict(status='PASS_EXACT_ALL_KNOWN_TRIANGLE_REPRESENTATIVES',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        counts=dict(literal_representatives=len(occurrences),selected_triangle_rows=len(gate['records']),
            old_unique_stars=original_stars,all_unique_stars=len(stars),new_stars=len(stars)-original_stars,
            old_discriminant_classes=len(gate['classes']),all_discriminant_classes=len(classes),provenance_kinds=dict(counts)),
        classes=classes,stars=stars,selected_references=references,occurrences=occurrences,
        scope='All46032 provenance entries of the3649 retained dual rows are reconstructed. Exact graph isomorphisms transfer normalized Gram, maximal-simplex and determinant data; every genuinely new star is classified separately. This does not expand compressed old-ADE provenance into all underlying catalogue occurrences, cover other triangle rows, classify all simplex stars, or cover quadrilaterals and absent vertex forests.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()))
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
