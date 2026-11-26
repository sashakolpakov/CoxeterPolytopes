"""Exact, repeatable cone presolve with replayable equivalence witnesses.

A pass is a checkpoint, never evidence of convergence.  No coefficient
is discarded by tolerance.  The recognized faces are scalar positive
circuits, zero PSD diagonals, and constant symbolic column kernels.
"""
from fractions import Fraction as F

from . import conic, feasibility


class Contradiction(ValueError):
    def __init__(self, proof):
        super().__init__('Exact conic contradiction')
        self.proof = proof


def substitute(q, substitutions):
    answer = {}
    for v,c in q.items():
        replacement = substitutions[v] if v else {'':F(1)}
        for w,a in replacement.items():
            answer[w] = answer.get(w,F(0))+c*a
            if not answer[w]: del answer[w]
    return answer


def linear_elimination(variables, equations):
    """Triangular sparse elimination; avoid quadratic scans of unit pivots."""
    order = {v:i for i,v in enumerate(variables)}
    basis, history = {}, []
    for k,q in sorted(enumerate(equations),key=lambda item:len(item[1])):
        row, operations = dict(q), []
        while True:
            present = [v for v in row if v in basis]
            if not present: break
            pivot = min(present,key=order.get)
            factor = row[pivot]
            row = conic.plus(row,conic.times(-factor,basis[pivot]))
            operations.append((pivot,factor))
        names = [v for v in row if v]
        if not names:
            if row:
                raise Contradiction(dict(kind='inconsistent_linear_equation',equation=k,
                                         operations=operations,remainder=row,history=history))
            continue
        pivot = min(names,key=order.get)
        scale = row[pivot]
        basis[pivot] = conic.times(1/scale,row)
        history.append(dict(equation=k,operations=operations,pivot=pivot,scale=scale,
                            normalized_row=basis[pivot]))
    substitutions = {}
    for v in reversed(variables):
        substitutions[v] = substitute(conic.times(-1,{w:c for w,c in basis[v].items() if w!=v}),
                                      substitutions) if v in basis else {v:F(1)}
    verify_elimination(variables,equations,substitutions,history)
    return substitutions,history


def verify_elimination(variables,equations,substitutions,history):
    basis = {}
    for item in history:
        row = dict(equations[item['equation']])
        for v,c in item['operations']:
            if v not in basis or row.get(v,F(0))!=c:
                raise ValueError('Invalid elimination operation')
            row = conic.plus(row,conic.times(-c,basis[v]))
        pivot,scale = item['pivot'],item['scale']
        if not scale or pivot in basis or row.get(pivot)!=scale:
            raise ValueError('Invalid elimination pivot')
        row = conic.times(1/scale,row)
        if row!=item['normalized_row']:
            raise ValueError('Changed normalized elimination row')
        basis[pivot] = row
    expected = {}
    for v in reversed(variables):
        expected[v] = substitute(conic.times(-1,{w:c for w,c in basis[v].items() if w!=v}),expected) \
            if v in basis else {v:F(1)}
    if expected!=substitutions or any(substitute(q,substitutions) for q in equations):
        raise ValueError('The elimination map is not exactly equivalent')


def test_row(problem,ref):
    b,i = ref
    return problem['blocks'][b]['entries'][i][i]


def discover(problem):
    """Cheap exact face tests, including every PSD diagonal as a nonnegative row."""
    rows,units,seen,equations,proofs = [],{},{},{},[]
    def equation(q):
        if q:
            first = q[min(q)]
            normalized = conic.times(1/first,q)
            equations.setdefault(conic.canonical(normalized),normalized)
    def circuit(combination):
        if any(w<=0 for _,w in combination) or conic.plus(*(
                conic.times(w,test_row(problem,ref)) for ref,w in combination)):
            raise ValueError('Invalid positive cone identity')
        proofs.append(dict(kind='positive_zero_sum',combination=combination))
        for ref,_ in combination: equation(test_row(problem,ref))
    for b,block in enumerate(problem['blocks']):
        M = block['entries']
        for i,row in enumerate(M):
            q = row[i]
            if not q:
                if any(row):
                    proofs.append(dict(kind='PSD_zero_diagonal',block=b,index=i))
                    for entry in row: equation(entry)
                continue
            if not any(qv for qv in q):
                if q.get('',F(0))<0:
                    raise Contradiction(dict(kind='negative_diagonal_constant',ref=(b,i),row=q))
                continue
            ref = (b,i)
            rows.append((ref,q))
            if len(q)==1 and '' not in q:
                v,c = next(iter(q.items()))
                units.setdefault(v,{})[1 if c>0 else -1] = (ref,c)
    for ref,q in rows:
        scale = abs(q[min(q)])
        normalized = conic.times(1/scale,q)
        key = conic.canonical(normalized)
        negative = conic.canonical(conic.times(-1,normalized))
        if negative in seen:
            old,old_scale = seen[negative]
            circuit([(old,1/old_scale),(ref,1/scale)])
        seen.setdefault(key,(ref,scale))
        if q.get('',F(0))>0 or len(q)==1: continue
        combination = [(ref,F(1))]
        for v,c in q.items():
            if not v: continue
            option = units.get(v,{}).get(-1 if c>0 else 1)
            if option is None: break
            other,a = option
            combination.append((other,-c/a))
        else:
            total = conic.plus(*(conic.times(w,test_row(problem,r)) for r,w in combination))
            if total.get('',F(0))<0:
                raise Contradiction(dict(kind='negative_positive_combination',combination=combination,
                                         constant=total['']))
            circuit(combination)
    return list(equations.values()),proofs


def column_compression(M):
    if len(M)==1:
        return ([0],[[F(1)]]) if M[0][0] else ([],[[]])
    kernels = feasibility.symbolic_column_kernels(M)
    relations = {max(i for i,_ in d):dict(d) for d in kernels}
    keep = [i for i in range(len(M)) if i not in relations]
    positions = {v:i for i,v in enumerate(keep)}
    E = []
    for i in range(len(M)):
        row = [F(0)]*len(keep)
        if i in positions: row[positions[i]]=F(1)
        else:
            d = relations[i]
            for j,c in d.items():
                if j==i: continue
                if j>=i: raise ValueError('Nontriangular symbolic kernel')
                row = [a-c*b/d[i] for a,b in zip(row,E[j],strict=True)]
        E.append(row)
    K = [[M[i][j] for j in keep] for i in keep]
    verify_compression(M,K,E)
    return keep,E


def verify_compression(M,K,E):
    if len(E)!=len(M) or any(len(row)!=len(K) for row in E):
        raise ValueError('Wrong congruence shape')
    for i,row in enumerate(M):
        for j in range(i,len(M)):
            expected = conic.plus(*(conic.times(a*b,K[k][l]) for k,a in enumerate(E[i]) if a
                                   for l,b in enumerate(E[j]) if b))
            if row[j]!=expected: raise ValueError('A cone reconstruction identity failed')


def block_key(M):
    first = next((q[min(q)] for row in M for q in row if q),None)
    if first is None: return None,F(1)
    scale = abs(first)
    normalized = [[conic.times(1/scale,q) for q in row] for row in M]
    return conic.digest_bytes(conic.canonical(normalized)),scale


def pass_once(problem):
    inferred,faces = discover(problem)
    equations = [*problem['equations'],*inferred]
    updates,history = linear_elimination(problem['variables'],equations)
    blocks,maps,seen = [],[],{}
    counts = dict(forced_equations=len(inferred),eliminated_variables=len(history),
                  removed_constant_cones=0,removed_duplicate_cones=0,removed_PSD_dimensions=0)
    for b,block in enumerate(problem['blocks']):
        old = block['entries']
        M = [[{} for _ in old] for _ in old]
        for i,row in enumerate(old):
            for j in range(i,len(old)):
                M[i][j]=M[j][i]=substitute(row[j],updates)
        if not any(v for row in M for q in row for v in q):
            matrix = [[q.get('',F(0)) for q in row] for row in M]
            try: pivots = conic.psd_pivots(matrix)
            except ValueError as error:
                raise Contradiction(dict(kind='negative_constant_after_elimination',block=b,matrix=matrix,
                    faces=faces,equations=equations,elimination=history,substitutions=updates)) from error
            maps.append(dict(source=b,kind='constant_PSD',pivots=pivots))
            counts['removed_constant_cones']+=1
            continue
        keep,E = column_compression(M)
        K = [[M[i][j] for j in keep] for i in keep]
        counts['removed_PSD_dimensions']+=len(M)-len(K)
        key,scale = block_key(K)
        if key in seen:
            target,old_scale = seen[key]
            ratio = scale/old_scale
            target_M = blocks[target]['entries']
            if K!=[[conic.times(ratio,q) for q in row] for row in target_M]:
                raise ValueError('False positive scalar-ray identity')
            counts['removed_duplicate_cones']+=1
        else:
            target,ratio = len(blocks),F(1)
            seen[key]=(target,scale)
            blocks.append(dict(name=block['name'],entries=K))
        maps.append(dict(source=b,kind='principal_congruence',target=target,scale=ratio,
                         keep=keep,embedding=E))
    used = {v for block in blocks for row in block['entries'] for q in row for v in q if v}
    free = [v for v in problem['variables'] if updates[v]=={v:F(1)}]
    null = [v for v in free if v not in used]
    if null:
        removed = set(null)
        updates = {v:{w:c for w,c in q.items() if w not in removed} for v,q in updates.items()}
    output = conic.model([v for v in free if v in used],blocks)
    proof = dict(schema='exact-cone-reduction-pass-v1',faces=faces,inferred_equations=inferred,
                 elimination=history,substitutions=updates,null_directions_fixed_zero=null,
                 block_maps=maps,counts=counts)
    verify_pass(problem,output,proof)
    return output,proof


def verify_pass(source,target,proof):
    # Independently replay source-cone identities before trusting any substitution.
    allowed = set()
    for item in proof['faces']:
        if item['kind']=='positive_zero_sum':
            combination = item['combination']
            if any(w<=0 for _,w in combination) or conic.plus(*(
                    conic.times(w,test_row(source,r)) for r,w in combination)):
                raise ValueError('Invalid retained positive-face proof')
            allowed.update(conic.canonical(conic.times(1/q[min(q)],q))
                           for r,_ in combination if (q:=test_row(source,r)))
        elif item['kind']=='PSD_zero_diagonal':
            row = source['blocks'][item['block']]['entries'][item['index']]
            if row[item['index']]: raise ValueError('A claimed zero diagonal is nonzero')
            allowed.update(conic.canonical(conic.times(1/q[min(q)],q)) for q in row if q)
        else: raise ValueError('Unknown cone-face proof')
    if any(conic.canonical(q) not in allowed for q in proof['inferred_equations']):
        raise ValueError('Unproved equation inserted in presolve')
    equations = [*source['equations'],*proof['inferred_equations']]
    raw_updates,replayed = linear_elimination(source['variables'],equations)
    if conic.canonical(replayed)!=conic.canonical(proof['elimination']):
        raise ValueError('Changed retained elimination history')
    null = set(proof['null_directions_fixed_zero'])
    if any(raw_updates[v]!={v:F(1)} for v in null):
        raise ValueError('A removed gauge is not a free coordinate')
    expected = {v:{w:c for w,c in q.items() if w not in null} for v,q in raw_updates.items()}
    if expected!=proof['substitutions']:
        raise ValueError('Changed reconstruction map')
    if [m['source'] for m in proof['block_maps']]!=list(range(len(source['blocks']))):
        raise ValueError('Source constraints lost from reduction proof')
    referenced = set()
    for item,block in zip(proof['block_maps'],source['blocks'],strict=True):
        if null and any(v in null for row in block['entries'] for q in row
                        for v in substitute(q,raw_updates)):
            raise ValueError('A removed gauge changes a source cone')
        M = [[substitute(q,expected) for q in row] for row in block['entries']]
        if item['kind']=='constant_PSD':
            if any(v for row in M for q in row for v in q): raise ValueError('Nonconstant deleted cone')
            conic.psd_pivots([[q.get('',F(0)) for q in row] for row in M])
        elif item['kind']=='principal_congruence':
            K = target['blocks'][item['target']]['entries']
            if item['scale']<=0: raise ValueError('Nonpositive duplicate scale')
            scaled = [[conic.times(item['scale'],q) for q in row] for row in K]
            if [[M[i][j] for j in item['keep']] for i in item['keep']]!=scaled:
                raise ValueError('Changed principal submatrix')
            verify_compression(M,scaled,item['embedding'])
            referenced.add(item['target'])
        else: raise ValueError('Unknown block deletion rule')
    if referenced!=set(range(len(target['blocks']))): raise ValueError('Unproved target cone')
    if any(v not in target['variables'] for q in expected.values() for v in q if v):
        raise ValueError('A reconstruction variable is missing')
    for v in target['variables']:
        if expected[v]!={v:F(1)}: raise ValueError('A target coordinate is not free')
    return dict(all_source_constraints_accounted_for=True,all_congruences_verified_exactly=True,
                all_face_equations_have_positive_source_proofs=True)


def audit(problem):
    inferred,_ = discover(problem)
    counts = dict(equalities=len(problem['equations']),pending_face_equations=len(inferred),
                  constant_cones=0,zero_PSD_dimensions=0,duplicate_cones=0,unused_variables=0)
    seen,used = set(),set()
    for block in problem['blocks']:
        M = block['entries']
        used.update(v for row in M for q in row for v in q if v)
        if not any(v for row in M for q in row for v in q):
            conic.psd_pivots([[q.get('',F(0)) for q in row] for row in M])
            counts['constant_cones']+=1
            continue
        keep,_ = column_compression(M)
        counts['zero_PSD_dimensions']+=len(M)-len(keep)
        key,_ = block_key(M)
        if key in seen: counts['duplicate_cones']+=1
        seen.add(key)
    counts['unused_variables']=sum(v not in used for v in problem['variables'])
    return dict(schema='conic-fixed-point-audit-v1',passed=not any(counts.values()),counts=counts,
                scope='exact scalar positive circuits, zero diagonals, constant symbolic kernels, constants and duplicates',
                no_claim_of_minimal_face_or_numerical_condition_number=True)
