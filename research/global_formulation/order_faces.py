"""Exact scalar faces from moment-order cycles and positive row circuits.

An edge f->g is admitted only if g-f is a positive multiple of a retained
source scalar cone. A circulation gives a checked positive zero identity;
all its rows must vanish. This is equivalent conic presolve, not a new
geometric assumption or a selected list of diagrams.
"""
from collections import deque
from fractions import Fraction as F
import json
from pathlib import Path
import time

from . import conic,feasibility,sdpa_double


def find(original,problem,substitutions):
    scalar={}
    for b,block in enumerate(problem['blocks']):
        if len(block['entries'])!=1:continue
        q=block['entries'][0][0]
        if not q:continue
        size=max(abs(c) for c in q.values())
        scalar.setdefault(conic.canonical(conic.times(1/size,q)),(b,size))
    nodes=[];node_index={}
    def node(q):
        key=conic.canonical(q)
        if key not in node_index:node_index[key]=len(nodes);nodes.append(q)
        return node_index[key]
    edges=[];seen=set()
    for block in original['blocks']:
        if len(block['entries'])!=1:continue
        q=block['entries'][0][0]
        if len(q)!=2 or '' in q:continue
        positive=[v for v,c in q.items() if c>0];negative=[v for v,c in q.items() if c<0]
        if len(positive)!=1 or len(negative)!=1 or q[positive[0]]!=-q[negative[0]]:continue
        high=substitutions[positive[0]];low=substitutions[negative[0]]
        difference=conic.plus(high,conic.times(-1,low))
        if not difference:continue
        size=max(abs(c) for c in difference.values())
        match=scalar.get(conic.canonical(conic.times(1/size,difference)))
        if match is None:continue
        b,source_size=match;u=node(low);v=node(high);key=(u,v,b)
        if key in seen:continue
        seen.add(key);scale=size/source_size
        if conic.times(scale,problem['blocks'][b]['entries'][0][0])!=difference:
            raise ValueError('An order edge is not its asserted positive scalar cone')
        edges.append((u,v,b,scale))
    outgoing=[[] for _ in nodes];incoming=[[] for _ in nodes]
    for e,(u,v,_,_) in enumerate(edges):outgoing[u].append((v,e));incoming[v].append((u,e))
    visited=set();finish=[]
    for start in range(len(nodes)):
        if start in visited:continue
        visited.add(start);stack=[(start,0)]
        while stack:
            u,k=stack[-1]
            if k==len(outgoing[u]):finish.append(u);stack.pop();continue
            stack[-1]=(u,k+1);v,_=outgoing[u][k]
            if v not in visited:visited.add(v);stack.append((v,0))
    component={};groups=[]
    for start in reversed(finish):
        if start in component:continue
        c=len(groups);group=[];stack=[start];component[start]=c
        while stack:
            u=stack.pop();group.append(u)
            for v,_ in incoming[u]:
                if v not in component:component[v]=c;stack.append(v)
        groups.append(group)
    internal={}
    for e,(u,v,_,_) in enumerate(edges):
        if component[u]==component[v] and u!=v:internal.setdefault(component[u],[]).append(e)
    weights={};summaries=[]
    for c,edge_ids in internal.items():
        group=groups[c];root=group[0]
        def tree(adjacency):
            order=[root];parents={};found={root};queue=deque([root])
            while queue:
                u=queue.popleft()
                for v,e in adjacency[u]:
                    if component[v]==c and v not in found:
                        found.add(v);parents[v]=(u,e);order.append(v);queue.append(v)
            if len(found)!=len(group):raise ValueError('An asserted strong component has no spanning tree')
            return order,parents
        out_order,out_parent=tree(outgoing);in_order,in_parent=tree(incoming)
        flow={e:1 for e in edge_ids};out_count={v:0 for v in group};in_count=dict(out_count)
        for e in edge_ids:u,v,_,_=edges[e];out_count[u]+=1;in_count[v]+=1
        for order,parents,counts in ((out_order,out_parent,out_count),(in_order,in_parent,in_count)):
            for v in reversed(order[1:]):
                parent,e=parents[v];flow[e]+=counts[v];counts[parent]+=counts[v]
        balance={v:0 for v in group}
        for e,value in flow.items():
            u,v,b,scale=edges[e];balance[u]-=value;balance[v]+=value
            weights[b]=weights.get(b,F(0))+value*scale
        if any(balance.values()):raise ValueError('The positive circulation is not balanced')
        summaries.append(dict(nodes=len(group),edges=len(edge_ids)))
    identity=conic.plus(*(conic.times(w,problem['blocks'][b]['entries'][0][0]) for b,w in weights.items()))
    if identity or any(w<=0 for w in weights.values()):raise ValueError('The complete order-face identity is false')
    equations=[problem['blocks'][b]['entries'][0][0] for b in sorted(weights)]
    proof=dict(positive_zero_identity=[dict(block=b,weight=weights[b]) for b in sorted(weights)],
        every_scalar_source_identity_checked=True,strong_components=summaries,
        total_order_nodes=len(nodes),total_order_edges=len(edges),forced_scalar_rows=len(equations))
    return equations,proof


def three_row_circuits(problem):
    """Positive three-row circuits with a one/two-term summand.

    The two longer rows are searched when they have 3--64 terms.
    Deleting one coefficient from each of those rows leaves
opposite proportional forms. Their residual is looked up in the full
scalar row table. Every resulting rational identity is checked afresh.
"""
    rows={b:block['entries'][0][0] for b,block in enumerate(problem['blocks'])
        if len(block['entries'])==1 and block['entries'][0][0]}
    normalized={}
    for b,q in rows.items():
        size=max(abs(c) for c in q.values())
        normalized.setdefault(conic.canonical(conic.times(1/size,q)),(b,size))
    signatures={};seen=set();weights={};circuits=0
    def record(combination):
        nonlocal circuits
        key=tuple(sorted(b for b,_ in combination))
        if key in seen:return
        if any(w<=0 for _,w in combination) or conic.plus(*(
                conic.times(w,rows[b]) for b,w in combination)):
            raise ValueError('A proposed positive scalar circuit failed its exact identity')
        seen.add(key);circuits+=1
        for b,w in combination:weights[b]=weights.get(b,F(0))+w
    for b,q in rows.items():
        if not 3<=len(q)<=64:continue
        for deleted in q:
            remaining={v:c for v,c in q.items() if v!=deleted}
            size=max(abs(c) for c in remaining.values())
            key=conic.canonical(conic.times(1/size,remaining))
            opposite=conic.canonical(conic.times(-1/size,remaining))
            for old,old_size in signatures.get(opposite,[]):
                residual=conic.plus(conic.times(1/size,q),conic.times(1/old_size,rows[old]))
                if not residual:record([(b,1/size),(old,1/old_size)]);continue
                if len(residual)>2:raise ValueError('A deleted-coordinate signature left a large residual')
                rsize=max(abs(c) for c in residual.values())
                match=normalized.get(conic.canonical(conic.times(-1/rsize,residual)))
                if match is not None:
                    other,other_size=match
                    record([(b,1/size),(old,1/old_size),(other,rsize/other_size)])
            signatures.setdefault(key,[]).append((b,size))
    if conic.plus(*(conic.times(w,rows[b]) for b,w in weights.items())):
        raise ValueError('The complete scalar-circuit exposure is false')
    return [rows[b] for b in sorted(weights)],dict(
        positive_zero_identity=[dict(block=b,weight=weights[b]) for b in sorted(weights)],
        every_positive_circuit_checked=True,positive_circuits=circuits,forced_scalar_rows=len(weights))


def run(root,out,source,candidate_run):
    out.mkdir();started=time.monotonic()
    original_path=root/'results/global_shared_cusp_current/necessary_primal.exact.json'
    model_path=source/'equilibrated.exact.json';subs_path=source/'equilibrated_substitutions.exact.json'
    native_path=candidate_run/'native.result'
    paths=[original_path,model_path,subs_path,native_path,Path(__file__),root/'scripts/recover_sdpa_order_faces_20261125.py',
        root/'research/global_formulation/conic.py',root/'research/global_formulation/feasibility.py']
    before=sdpa_double.binding(paths,root);conic.write_json(out/'source_bindings.json',before)
    original=conic.read_model(original_path);problem=conic.read_model(model_path)
    substitutions={v:{w:F(c) for w,c in q.items()} for v,q in json.loads(subs_path.read_text()).items()}
    raw,phase=conic.native_vector(native_path)
    if len(raw)!=len(problem['variables'])+1:raise ValueError('Use the matching common-margin native point')
    values=dict(zip(problem['variables'],raw[:-1],strict=True))
    equations,proof=find(original,problem,substitutions)
    additional,circuits=three_row_circuits(problem)
    equations=list({conic.canonical(q):q for q in [*equations,*additional]}.values())
    proof['positive_scalar_circuits']=circuits
    proof['forced_scalar_rows']=len(equations)
    conic.write_json(out/'verified_order_faces.exact.json',proof)
    print(json.dumps(dict(stage='exact_order_faces',forced_scalar_rows=len(equations),
        strong_components=len(proof['strong_components']))),flush=True)
    refined,updates=conic.eliminate(problem,equations)
    combined={v:feasibility.substitute(q,updates) for v,q in substitutions.items()}
    report=dict(status=('EXACT_ORDER_FACES_EXPOSED_ORIGINAL_FEASIBILITY_UNDECIDED' if equations
        else 'NO_NEW_SCALAR_FACE_OR_ORIGINAL_POINT'),
        original_input_sha256=before[str(original_path.relative_to(root))]['sha256'],
        forced_scalar_rows=len(equations),source_variables=len(problem['variables']),
        refined_variables=len(refined['variables']),all_reductions_implied_by_source_cones=True,
        accepted_bounds_changed=False,attempts=[])
    def try_point(P,S):
        try:
            point,checked=conic.reconstruct(original,P,S,[values[v] for v in P['variables']])
            conic.write_json(out/'original_feasible_point.exact.json',point)
            report.update(status='EXACT_ORIGINAL_SDP_FEASIBLE',exact_check=checked)
            return True
        except ValueError as error:report['attempts'].append(str(error));return False
    if not try_point(refined,combined):
        def progress(stage,data):
            with (out/'progress.jsonl').open('a') as stream:
                stream.write(json.dumps(conic.plain(dict(stage=stage,**data)))+'\n')
        next_model,next_subs,_,trace=feasibility.prepare(refined,progress,linear_faces=False,seconds=180)
        combined={v:feasibility.substitute(q,next_subs) for v,q in combined.items()};refined=next_model
        conic.write_json(out/'additional_exact_cone_reductions.json',trace)
        try_point(refined,combined)
    conic.write_json(out/'equilibrated.exact.json',refined)
    conic.write_json(out/'equilibrated_substitutions.exact.json',combined)
    report['refined_variables']=len(refined['variables'])
    if sdpa_double.binding(paths,root)!=before:raise ValueError('An exact order-recovery input changed')
    report['elapsed_seconds']=time.monotonic()-started
    conic.write_json(out/'result.json',report)
    print(json.dumps(conic.plain({k:v for k,v in report.items() if k!='exact_check'})),flush=True)
