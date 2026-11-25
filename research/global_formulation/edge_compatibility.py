"""One common spherical edge population for endpoints and incident faces.

All equations are incidence pushforwards of the SAME edge diagrams. Rank
tails are symbolic. No independent upper envelope of a deletion statistic
is substituted for its value on that population.
"""
from . import conic, data, supporting_counts as counts, vertices, white_profiles as white
from .moments import dimension_domain
from .polynomial import as_poly, var


def spherical_tail(q,kind,power):
    """t^(power+1) times the sum of q(n,r) over the tail components."""
    t = var("t"); result = as_poly(0)
    for m,c in q.items():
        nd,rd = m.count("N"),m.count("r")
        if any(v not in ("N","r") for v in m) or rd>3: raise ValueError("Unsupported component statistic")
        field = "count" if not rd else ("rank","square","cube")[rd-1]
        exponent = power-nd if not rd else power+1-nd-rd
        if exponent<0: raise ValueError("A tail normalization would introduce a pole")
        result += c*t**exponent*var(f"tail_{kind}_{field}")
    return result


def endpoint_columns():
    """Normalized sums over all edge ends at a finite vertex or full cusp."""
    N,r = var("N"),var("r")
    finite,ideal = ({(T,j):as_poly(0) for T in white.FAMILIES for j in (0,1)} for _ in range(2))
    for kind,rank in data.finite_types():
        deleted,_ = counts.deletion_histograms(kind,rank)
        u = var(f"u_{kind}{rank}")
        for T,j in finite:
            retained = (1 if not j else rank) if counts.name(kind,rank)==T else 0
            finite[T,j] += counts.fixed_normalize((N-rank)*retained+deleted[T,j],u,4)
    for kind in data.FINITE_TAIL:
        deleted,_ = counts.classical(kind)
        for T,j in finite:
            retained = (as_poly(1) if not j else r) if T==kind else as_poly(0)
            finite[T,j] += spherical_tail((N-r)*retained+deleted[T,j],kind,4)
    for kind,rank in data.affine_types():
        deleted,_ = counts.deletion_histograms(kind,rank,True)
        u = var(f"u_{kind}{rank}")
        for key,q in deleted.items(): ideal[key] += counts.fixed_normalize(as_poly(q),u,2,rank+1)
    for kind in ("A","B","C","D"):
        deleted,_ = counts.classical(kind,True)
        for key,q in deleted.items(): ideal[key] += counts.ideal_tail_normalize(q,kind)
    return finite,ideal


def extend_edge(domain):
    """The whole diagram supporting an edge is spherical of rank n-1."""
    t = var("t")
    u = vertices.counters(domain,data.finite_types())
    tails = {T:vertices.tail(domain,T,minimum) for T,(minimum,_) in data.FINITE_TAIL.items()}
    vertices.rank_data(domain,u,tails,1-t)
    family = {(T,j):as_poly(0) for T in white.FAMILIES for j in (0,1)}
    deletion = dict(family); N,r = var("N"),var("r")
    for (kind,rank),counter in u.items():
        T = counts.name(kind,rank)
        family[T,0] += counter; family[T,1] += rank*counter
        one,_ = counts.deletion_histograms(kind,rank)
        for U,j in deletion:
            retained = (1 if not j else rank) if T==U else 0
            raw = (N-1-rank)*retained+one[U,j]
            deletion[U,j] += counts.fixed_normalize(raw,counter,1)
    for kind,tail in tails.items():
        family[kind,0] += tail["count"]; family[kind,1] += tail["rank"]
        one,_ = counts.classical(kind)
        for T,j in deletion:
            retained = (as_poly(1) if not j else r) if T==kind else as_poly(0)
            deletion[T,j] += spherical_tail((N-1-r)*retained+one[T,j],kind,1)
    domain.eq("same_complete_edge_component_count",sum((family[T,0] for T in white.FAMILIES),as_poly(0))-var("edge_components"))
    domain.edge_family_columns,domain.edge_deleted_columns = family,deletion


def tests(domain,side=None):
    """Declared finite common-state tests; no arbitrary incidence products."""
    if side is None: names=domain.star_common
    else:
        names=("star_white_components",f"star_scale{side}",f"star_degree{side}")
        names += tuple(n for n in domain.star_common if n.startswith(("eta","right","c")))
    return [as_poly(1),*[var(name) for name in names]]


def couple_face(program,face,profile,sides):
    t = var("t")
    def eq(name,*forms): program.eq(name,conic.plus(*forms))
    def negative(q): return conic.times(-1,q)
    state = face.domain.edge_host_data
    for k,test in enumerate(tests(face.domain)):
        for (T,j),whole in state["white"].items():
            local = var("white_type_"+T)*(t if not j else var("white_rank"))
            eq(f"edge_white_population{sides}_{T}_{j}_test{k}",
                profile.linear(test*local),negative(face.linear(test*whole)))
    for i,host in enumerate(state["hosts"]):
        for k,test in enumerate(tests(face.domain,i)):
            for (T,j),deleted in host["deletion"].items():
                local = var(f"attach{i}")*var("white_type_"+T)*(t if not j else var("white_rank"))
                eq(f"edge_host_deleted{sides}_{i}_{T}_{j}_test{k}",
                    profile.linear(test*local),negative(face.linear(test*deleted)))


def couple(program,finite,ideal,all_faces,profiles,edge):
    t = var("t"); endpoint = endpoint_columns()
    def eq(name,*forms): program.eq(name,conic.plus(*forms))
    def negative(q): return conic.times(-1,q)
    for sides,face in all_faces.items(): couple_face(program,face,profiles[sides],sides)
    for key,column in edge.domain.edge_family_columns.items():
        T,j = key
        # Both endpoints carry exactly this same spherical edge diagram.
        eq(f"same_edge_family_endpoints_{T}_{j}",conic.times(2,edge.linear(column)),
            negative(finite.linear(endpoint[0][key])),negative(ideal.linear(endpoint[1][key])))
        # Every edge has n-1 incident two-faces, including at ideal ends.
        side_terms=[]; deleted_terms=[]
        for sides,face in all_faces.items():
            state = face.domain.edge_host_data
            indices=range(sides) if sides in (3,4) else (0,)
            side=sum((state["edge"][i][key] for i in indices),as_poly(0))
            side_terms.append(face.linear(t*side))
            deleted_terms.append(face.linear(t*len(tuple(indices))*state["white"][key]))
        eq(f"same_edge_family_incident_faces_{T}_{j}",*side_terms,negative(edge.linear((1-t)*column)))
        eq(f"same_edge_family_white_deletions_{T}_{j}",*deleted_terms,
            negative(edge.linear(edge.domain.edge_deleted_columns[key])))


def self_checks():
    """Exact product-cusp/forest normalization and deletion commutation."""
    from fractions import Fraction as F
    from itertools import product
    def evaluate(q,point):
        result=F(0)
        for m,c in q.items():
            for v in m: c*=point.get(v,F(0))
            result+=c
        return result
    def population(components,n):
        point=dict(t=F(1,n))
        for kind,rank in components:
            if (kind,rank) in data.finite_types():
                name=f"u_{kind}{rank}";point[name]=point.get(name,F(0))+F(1,n)
            else:
                for k,suffix in enumerate(("count","rank","square","cube")):
                    name=f"tail_{kind}_{suffix}"
                    point[name]=point.get(name,F(0))+(F(1,n) if not k else F(rank,n)**k)
        return point
    def block_graph(components,affine=False):
        graph=[]
        for kind,rank in components:
            offset=len(graph)
            graph += [{j+offset:label for j,label in row.items()}
                for row in counts.labelled_graph(kind,rank,affine)]
        return graph
    columns=endpoint_columns(); forests=0
    for components in ((("A",1),("E",6),("D",4)),(("A",21),("D",36)),(("B",24),("B",3),("F",4))):
        n=sum(r for _,r in components);point=population(components,n)
        literal=counts.graph_histograms(block_graph(components))[0]
        for key,q in columns[0].items():
            if evaluate(q,point)*n**5!=literal[key]: raise ValueError("Finite endpoint normalization failed")
        forests+=1
    cusps=0
    for components in ((("A",1),("F",4)),(("A",25),("C",26)),(("B",5),("D",6))):
        n=1+sum(r for _,r in components);V=1;point=dict(t=F(1,n));literal={key:F(0) for key in columns[1]}
        for _,r in components: V*=r+1
        for kind,r in components:
            one,_=counts.deletion_histograms(kind,r,True)
            for key,q in one.items(): literal[key]+=V*F(q,r+1)
            if (kind,r) in data.affine_types():
                name=f"u_{kind}{r}";point[name]=point.get(name,F(0))+F(1,n)
            else:
                for k,suffix in enumerate(("count","rank","square","cube")):
                    name=f"tail_{kind}_{suffix}"
                    point[name]=point.get(name,F(0))+(F(1,n) if not k else F(r,n)**k)
                name=f"tail_{kind}_reciprocal";point[name]=point.get(name,F(0))+F(1,n*(r+1))
        for key,q in columns[1].items():
            if evaluate(q,point)*V*n**3!=literal[key]: raise ValueError("Cusp endpoint normalization failed")
        cusps+=1
    from .moments import dimension_domain
    edges=0
    for components in ((("A",1),("D",4)),(("A",21),("D",36)),(("E",7),("B",3))):
        n=1+sum(r for _,r in components);point=population(components,n)
        d=dimension_domain("edge_fixture",n);d.variable("edge_components");extend_edge(d)
        literal=counts.graph_histograms(block_graph(components))[0]
        for key,q in d.edge_deleted_columns.items():
            if evaluate(q,point)*n*n!=literal[key]: raise ValueError("Edge deletion normalization failed")
        edges+=1
    return dict(finite_endpoint_fixtures=forests,complete_product_cusp_endpoint_fixtures=cusps,
        whole_edge_deletion_fixtures=edges,all_families_have_one_common_edge_population=True)
