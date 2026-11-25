"""Complete family count/rank marginals under spherical/cusp deletions.

Long A/B/D and affine A/B/C/D components use all-rank polynomials. Small
connected inputs use their entire one/two-node deletion histograms; no
vertex forests or face-star catalogue is chosen or priced.
"""
from fractions import Fraction as F
from itertools import combinations

from . import data
from .polynomial import as_poly, var

FAMILIES = ("A","B","D","E6","E7","E8","F4","H3","H4","I2")


def labelled_graph(kind,rank,affine=False):
    adjacency = data.graph(kind,rank,affine)
    G = [{j:3 for j in row} for row in adjacency]
    m = len(G)
    def label(i,j,value): G[i][j] = G[j][i] = value
    if kind == "B": label(m-4 if affine else m-2,m-3 if affine else m-1,4)
    elif kind == "C": label(0,1,4); label(m-2,m-1,4)
    elif kind == "F": label(2 if affine else 1,3 if affine else 2,4)
    elif kind in ("H","I"): label(0,1,5)
    elif kind == "G": label(0,1,6)
    return G


def family(G,nodes):
    n = len(nodes)
    if n == 1: return "A"
    edges = [(i,j,G[i][j]) for i in nodes for j in G[i] if i < j and j in nodes]
    if n == 2:
        label = edges[0][2]
        return "A" if label == 3 else "B" if label == 4 else "I2"
    degree = {i:sum(j in nodes for j in G[i]) for i in nodes}
    labels = [v for _,_,v in edges if v != 3]
    if not labels:
        branches = [i for i in nodes if degree[i] == 3]
        if not branches and max(degree.values()) <= 2: return "A"
        if len(branches) == 1:
            center = branches[0]; arms = []
            for neighbor in G[center]:
                if neighbor not in nodes: continue
                previous,current,length = center,neighbor,1
                while degree[current] == 2:
                    following = next(j for j in G[current] if j in nodes and j != previous)
                    previous,current,length = current,following,length+1
                arms.append(length)
            arms.sort()
            if arms[:2] == [1,1]: return "D"
            if arms in ([1,2,2],[1,2,3],[1,2,4]): return "E"+str(n)
    if len(labels) == 1 and max(degree.values()) <= 2:
        i,j,label = next(edge for edge in edges if edge[2] != 3)
        endpoint = degree[i] == 1 or degree[j] == 1
        if label == 4:
            if endpoint: return "B"
            if n == 4: return "F4"
        if label == 5 and endpoint and n in (3,4): return "H"+str(n)
    raise ValueError("Unclassified connected spherical deletion component")


def deletion_histograms(kind,rank,affine=False):
    G = labelled_graph(kind,rank,affine)
    return graph_histograms(G)


def graph_histograms(G):
    result = [{(T,j):F(0) for T in FAMILIES for j in (0,1)} for _ in range(2)]
    for count in (1,2):
        for deleted in combinations(range(len(G)),count):
            remaining = set(range(len(G)))-set(deleted)
            while remaining:
                root = remaining.pop(); nodes = {root}; stack = [root]
                while stack:
                    for neighbor in G[stack.pop()]:
                        if neighbor in remaining:
                            remaining.remove(neighbor); nodes.add(neighbor); stack.append(neighbor)
                T = family(G,nodes)
                result[count-1][T,0] += 1
                result[count-1][T,1] += len(nodes)
    return result


def classical(kind,affine=False):
    r = var("r"); m = r+1
    a,b = ({(T,j):as_poly(0) for T in FAMILIES for j in (0,1)} for _ in range(2))
    if not affine:
        if kind == "A":
            a["A",0],a["A",1] = 2*(r-1),r*(r-1)
            b["A",0],b["A",1] = 3*(r-1)*(r-2)/2,r*(r-1)*(r-2)/2
        elif kind == "B":
            a["A",0],a["B",0] = r,r-2
            a["A",1],a["B",1] = r*(r-1)/2+1,r*(r-1)/2-1
            b["A",0],b["B",0] = r*(r-2),(r-2)*(r-3)/2
            b["A",1],b["B",1] = (r-2)*(r*r-r+3)/3,(r-2)*(r-3)*(r+2)/6
        else:
            a["A",0],a["D",0] = r+2,r-4
            a["A",1],a["D",1] = r*(r-1)/2+6,r*(r-1)/2-6
            b["A",0],b["D",0] = r*r-7,(r-4)*(r-5)/2
            b["D",1] = (r-4)*(r-5)*(r+6)/6
            b["A",1] = r*(r-1)*(r-2)/2-b["D",1]
    else:
        total_count,total_rank = 3*(m-1)*(m-2)/2,m*(m-1)*(m-2)/2
        if kind == "A":
            a["A",0],a["A",1] = m,m*(m-1)
            b["A",0],b["A",1] = m*(m-2),total_rank
        elif kind == "B":
            a["A",0],a["B",0],a["D",0] = as_poly(4),m-2,m-4
            a["A",1],a["B",1],a["D",1] = as_poly(6),m*(m-1)/2,m*(m-1)/2-6
            b["B",0],b["D",0] = (m-2)*(m-3)/2,(m-4)*(m-5)/2
            b["B",1],b["D",1] = (m-2)*(m-3)*(m+2)/6,(m-4)*(m-5)*(m+6)/6
            for j,total in ((0,total_count),(1,total_rank)): b["A",j] = total-b["B",j]-b["D",j]
        elif kind == "C":
            a["A",0],a["B",0] = as_poly(2),2*m-4
            a["A",1],a["B",1] = as_poly(2),m*(m-1)-2
            b["A",0],b["B",0] = (m-2)*(m+3)/2,(m-2)*(m-3)
            b["A",1] = (m-2)*(2+m*(m-1)/6)
            b["B",1] = total_rank-b["A",1]
        else:
            a["A",0],a["D",0] = as_poly(6),2*(m-4)
            a["A",1],a["D",1] = as_poly(10),m*(m-1)-10
            b["D",0],b["D",1] = (m-4)*(m-5), (m-4)*(m-5)*(m+6)/3
            b["A",0],b["A",1] = total_count-b["D",0],total_rank-b["D",1]
    return a,b


def name(kind,rank): return kind if kind in ("A","B","D") else "I2" if kind == "I" else kind+str(rank)


def fixed_normalize(q,counter,power,divisor=1):
    t = var("t")
    return sum((c/divisor*t**(power-m.count("N"))*counter for m,c in q.items()),as_poly(0))


def finite_tail_normalize(q,kind):
    t = var("t"); result = as_poly(0)
    for m,c in q.items():
        ndegree,rdegree = m.count("N"),m.count("r")
        field = "count" if not rdegree else ("rank","square","cube")[rdegree-1]
        result += c*t**(4-ndegree if not rdegree else 5-ndegree-rdegree)*var(f"tail_{kind}_{field}")
    return result


def ideal_tail_normalize(q,kind):
    # Exact division by r+1, leaving only the existing reciprocal moment.
    coefficients = {}
    for m,c in q.items():
        degree = m.count("r")
        coefficients[degree] = coefficients.get(degree,as_poly(0))+c*var("N")**m.count("N")
    quotient = {}
    for degree in range(max(coefficients,default=0),0,-1):
        value = coefficients.get(degree,as_poly(0))
        quotient[degree-1] = value
        coefficients[degree-1] = coefficients.get(degree-1,as_poly(0))-value
    t = var("t"); result = as_poly(0)
    for degree,qpart in quotient.items():
        field = "count" if not degree else ("rank","square","cube")[degree-1]
        for m,c in qpart.items():
            nd = m.count("N")
            result += c*t**(2-nd if not degree else 3-nd-degree)*var(f"tail_{kind}_{field}")
    for m,c in coefficients.get(0,as_poly(0)).items():
        result += c*t**(2-m.count("N"))*var(f"tail_{kind}_reciprocal")
    return result


def vertex_columns():
    N,r = var("N"),var("r")
    finite,ideal = ({(T,j):as_poly(0) for T in FAMILIES for j in (0,1)} for _ in range(2))
    for kind,rank in data.finite_types():
        a,b = deletion_histograms(kind,rank)
        u = var(f"u_{kind}{rank}")
        for T,j in finite:
            retained = (1 if not j else rank) if name(kind,rank) == T else 0
            raw = (N-rank)*(N-rank-1)*retained/2+(N-rank)*a[T,j]+b[T,j]
            finite[T,j] += fixed_normalize(raw,u,4)
    for kind in data.FINITE_TAIL:
        a,b = classical(kind)
        for T,j in finite:
            retained = (as_poly(1) if not j else r) if T == kind else as_poly(0)
            raw = (N-r)*(N-r-1)*retained/2+(N-r)*a[T,j]+b[T,j]
            finite[T,j] += finite_tail_normalize(raw,kind)
    for kind,rank in data.affine_types():
        a,b = deletion_histograms(kind,rank,True)
        u = var(f"u_{kind}{rank}")
        for T,j in ideal:
            raw = (N-1-rank)*a[T,j]/2+b[T,j]
            ideal[T,j] += fixed_normalize(raw,u,2,rank+1)
    for kind in ("A","B","C","D"):
        a,b = classical(kind,True)
        for T,j in ideal:
            ideal[T,j] += ideal_tail_normalize((N-1-r)*a[T,j]/2+b[T,j],kind)
    return finite,ideal


def self_checks():
    def value(q,rank): return sum((c*F(rank)**m.count("r") for m,c in q.items()),F(0))
    checked = 0
    for affine,kinds in ((False,("A","B","D")),(True,("A","B","C","D"))):
        for kind in kinds:
            ranks = (6,9,25,31) if affine else (4,7,21,29)
            for rank in ranks:
                literal = deletion_histograms(kind,rank,affine)
                formulas = classical(kind,affine)
                for q,table in zip(formulas,literal,strict=True):
                    for key in table:
                        if value(q[key],rank) != table[key]:
                            raise ValueError(f"Uniform deletion formula disagrees with the labelled graph: {kind}{rank} {affine} {key}")
                checked += 1
    # Every lower/exceptional connected input must classify completely.
    for kind,rank in data.finite_types(): deletion_histograms(kind,rank)
    for kind,rank in data.affine_types(): deletion_histograms(kind,rank,True)
    columns = vertex_columns()
    def evaluate(q,point):
        total = F(0)
        for monomial,c in q.items():
            for v in monomial: c *= point.get(v,F(0))
            total += c
        return total
    for types in ((("E",6),("A",1),("A",1)),(("A",21),("D",36))):
        n = sum(rank for _,rank in types); G = []; point = dict(t=F(1,n))
        for kind,rank in types:
            H = labelled_graph(kind,rank); offset = len(G)
            G += [{j+offset:label for j,label in row.items()} for row in H]
            if (kind,rank) in data.finite_types():
                key = f"u_{kind}{rank}"; point[key] = point.get(key,F(0))+F(1,n)
            else:
                for degree,field in enumerate(("count","rank","square","cube")):
                    key = f"tail_{kind}_{field}"
                    point[key] = point.get(key,F(0))+(F(1,n) if not degree else F(rank,n)**degree)
        literal = graph_histograms(G)[1]
        for key,q in columns[0].items():
            if evaluate(q,point)*n**5 != literal[key]:
                raise ValueError("A whole-vertex spherical family normalization failed: "+str(key))
    for types in ((("F",4),("A",1)),(("A",25),("C",26))):
        n = 1+sum(rank for _,rank in types); V = 1; point = dict(t=F(1,n)); literal = {(T,j):F(0) for T in FAMILIES for j in (0,1)}
        for kind,rank in types: V *= rank+1
        for kind,rank in types:
            a,b = deletion_histograms(kind,rank,True)
            for key in literal:
                literal[key] += V*(F(n-1-rank,2*(rank+1))*a[key]+F(1,rank+1)*b[key])
            if (kind,rank) in data.affine_types():
                key = f"u_{kind}{rank}"; point[key] = point.get(key,F(0))+F(1,n)
            else:
                for degree,field in enumerate(("count","rank","square","cube")):
                    key = f"tail_{kind}_{field}"
                    point[key] = point.get(key,F(0))+(F(1,n) if not degree else F(rank,n)**degree)
                point[f"tail_{kind}_reciprocal"] = point.get(f"tail_{kind}_reciprocal",F(0))+F(1,n*(rank+1))
        for key,q in columns[1].items():
            if evaluate(q,point)*V*n**3 != literal[key]:
                raise ValueError("A product-cusp family normalization failed: "+str(key))
    return dict(all_connected_deletion_inputs_classified=True,
        classical_deletion_polynomial_regressions=checked,
        whole_vertex_family_normalization_fixtures=4,
        formulas_proved_by_retained_edges_and_fork_cores_not_polynomial_fitting=True)
