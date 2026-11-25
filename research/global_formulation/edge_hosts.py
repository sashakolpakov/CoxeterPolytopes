"""Complete rooted spherical components for a wall adjoining a two-face.

The distinguished wall belongs to one connected component of the spherical
edge diagram. Deleting it gives precisely the white components it attaches
to. Long A/B/D ranks and root positions stay symbolic. The exceptional
tables contain every node of each complete connected classification input;
they are not selected vertex forests or face stars.
"""
from fractions import Fraction as F

from . import data, supporting_counts as counts, white_profiles as white
from .polynomial import as_poly, var


def deleted_at(graph, node):
    remaining = set(range(len(graph)))-{node}
    result = {(T,j):F(0) for T in white.FAMILIES for j in (0,1)}
    while remaining:
        first = remaining.pop(); component = {first}; stack = [first]
        while stack:
            for neighbor in graph[stack.pop()]:
                if neighbor in remaining:
                    remaining.remove(neighbor); component.add(neighbor); stack.append(neighbor)
        T = counts.family(graph,component)
        result[T,0] += 1; result[T,1] += len(component)
    return result


def exceptional_roots():
    """All node outcomes, grouped only when deletion and inverse data agree."""
    result = {}
    for T in ("E6","E7","E8","F4","H3","H4"):
        rank = white.FIXED[T]
        graph = counts.labelled_graph(T[0],rank)
        lower,upper = white.fixed_inverse_bounds(T)
        lengths = [2,2,1,1] if T == "F4" else [2]*rank
        grouped = {}
        for node in range(rank):
            deletion = deleted_at(graph,node)
            lo,hi = lengths[node]*lower[node][node],lengths[node]*upper[node][node]
            key = (tuple(deletion.items()),lo,hi)
            grouped.setdefault(key,[]).append(node)
        for index,((items,lo,hi),nodes) in enumerate(grouped.items()):
            result[T+"_orbit"+str(index)] = dict(family=T,rank=rank,
                nodes=nodes,deletion=dict(items),inverse_lower=lo,inverse_upper=hi)
    return result


def rooted(domain, side, exceptional):
    prefix = f"edge_host{side}_"
    t,a = var("t"),var(f"star_scale{side}")
    def variable(name,boolean=False): return domain.variable(prefix+name,boolean=boolean)
    def eq(name,q): domain.eq(prefix+name,q)
    def ge(name,q): domain.ge(prefix+name,q)
    kinds = ("A","B_short","B_near","B_tail","D_leaf","D_fork","D_near","D_tail","I2",*exceptional)
    flags = {name:variable(name,True) for name in kinds}
    eq("one_root_case",sum(flags.values(),as_poly(0))-1)
    # Sum of Boolean variables equal to one is pointwise one-hot. Their
    # moment basis and equation closure retain necessary consequences.
    family = {T:as_poly(0) for T in white.FAMILIES}
    for name,flag in flags.items():
        T = exceptional[name]["family"] if name in exceptional else name.split("_")[0]
        family[T] += flag
    R,x,v = variable("rank"),variable("position"),variable("projection_square")
    left,right = variable("left",True),variable("right",True)
    ge("rank_minimum",R-t); ge("edge_rank_capacity",1-t-R)
    eq("same_projection_square",v-a*a)
    ge("position_minimum",x-t); ge("position_maximum",R-x)
    # x=t for cases where a chain position is unused, including a D leaf.
    inactive = 1-family["A"]-family["B"]-(family["D"]-flags["D_leaf"])
    eq("unused_chain_position",inactive*(x-t))
    ge("left_nonempty",x-t-t*left); ge("left_empty",left-(x-t))
    eq("right_only_A",right*(1-family["A"]))
    ge("A_right_nonempty",family["A"]*(R-x-t*right))
    ge("A_right_empty",family["A"]*(right-R+x))
    eq("A_root_inverse",family["A"]*(2*v*x*(R+t-x)-t*(R+t)))
    ge("B_rank",family["B"]*(R-2*t))
    eq("B_short_position",flags["B_short"]*(x-R))
    eq("B_near_position",flags["B_near"]*(x-R+t))
    ge("B_tail_position",flags["B_tail"]*(R-x-2*t))
    eq("B_short_inverse",flags["B_short"]*(v*R-t))
    eq("B_long_inverse",(family["B"]-flags["B_short"])*(2*v*x-t))
    ge("D_rank",family["D"]*(R-4*t))
    eq("D_fork_position",flags["D_fork"]*(x-R+2*t))
    eq("D_near_position",flags["D_near"]*(x-R+3*t))
    ge("D_tail_position",flags["D_tail"]*(R-x-4*t))
    eq("D_leaf_inverse",flags["D_leaf"]*(v*R-2*t))
    eq("D_chain_inverse",(family["D"]-flags["D_leaf"])*(2*v*x-t))
    eq("I2_rank",family["I2"]*(R-2*t))
    # At a root of I2(m), a^2=sin^2(pi/m); m>=5 in the
    # canonical family. The continuous envelope includes every label.
    ge("I2_projection_upper",family["I2"]*(F(3,8)-v))
    deletion = {(T,j):as_poly(0) for T in white.FAMILIES for j in (0,1)}
    deletion["A",0] = t*(family["A"]*(left+right)+family["B"]*left+flags["B_near"]
        +family["D"]*left+flags["D_near"]+2*flags["D_fork"]+flags["D_leaf"]+family["I2"])
    deletion["A",1] = family["A"]*(R-t)+family["B"]*(x-t)+t*flags["B_near"]
    deletion["A",1] += family["D"]*(x-t)+3*t*flags["D_near"]+2*t*flags["D_fork"]
    deletion["A",1] += (R-t)*flags["D_leaf"]+t*family["I2"]
    deletion["B",0],deletion["B",1] = t*flags["B_tail"],(R-x)*flags["B_tail"]
    deletion["D",0],deletion["D",1] = t*flags["D_tail"],(R-x)*flags["D_tail"]
    for name,entry in exceptional.items():
        flag = flags[name]
        eq(name+"_rank",flag*(R-entry["rank"]*t))
        lo,hi = entry["inverse_lower"],entry["inverse_upper"]
        if lo == hi: eq(name+"_root_inverse",flag*(lo*v-1))
        else:
            ge(name+"_root_inverse_lower",flag*(1-lo*v))
            ge(name+"_root_inverse_upper",flag*(hi*v-1))
        for key,count in entry["deletion"].items(): deletion[key] += t*count*flag
    eq("deleted_rank",sum((deletion[T,1] for T in family),as_poly(0))-R+t)
    eq("same_side_degree",sum((deletion[T,0] for T in family),as_poly(0))-t*var(f"star_degree{side}"))
    # A spherical Gram inverse dominates the diagonal contribution of
    # every attachment, whose squared cosine is at least 1/4.
    gap = 1-a*a-var(f"star_degree{side}")/4
    ge("degree_projection_budget",gap)
    degree = var(f"star_degree{side}")
    eq("unattached_projection",(degree-1)*(degree-2)*(degree-3)*(1-a*a))
    domain.matrix(prefix+"degree_projection_localizer",[[gap,gap*a],[gap*a,gap*a*a]])
    domain.basis(prefix+"classification",[*(prefix+name for name in kinds),prefix+"rank",
        prefix+"position",prefix+"projection_square",prefix+"left",prefix+"right",f"star_scale{side}",f"star_degree{side}"],1)
    return dict(family=family,rank=R,deletion=deletion,flags=flags)


def extend_face(domain,sides):
    """Classify each side's whole connected host in the same face state."""
    t = var("t"); white_columns = {}
    for T in white.FAMILIES:
        for j in (0,1):
            name = f"white_population_{T}_{j}"
            white_columns[T,j] = domain.variable(name)
        c,r = white_columns[T,0],white_columns[T,1]
        minimum = {"A":1,"B":2,"D":4}.get(T,white.FIXED.get(T))
        domain.integer_counts[f"white_population_{T}_0"] = minimum
        domain.upper_bounds[f"white_population_{T}_0"] = F(1,minimum)
        domain.ge(f"white_population_{T}_minimum",r-minimum*c)
        domain.ge(f"white_population_{T}_rank_bound",(1-2*t)*c-t*r)
        if T in white.FIXED: domain.eq(f"white_population_{T}_fixed_rank",r-minimum*c)
    domain.rank_total = 1-2*t
    domain.eq("same_white_population_count",sum((white_columns[T,0] for T in white.FAMILIES),as_poly(0))-var("star_white_components"))
    domain.eq("same_white_population_rank",sum((white_columns[T,1] for T in white.FAMILIES),as_poly(0))-1+2*t)
    domain.basis("white_family_population",[f"white_population_{T}_{j}" for T,j in white_columns],1)
    hosts = [rooted(domain,i,exceptional_roots()) for i in range(sides)]
    edge_columns = []
    for i,host in enumerate(hosts):
        columns = {}
        for (T,j),q in white_columns.items():
            spectator = q-host["deletion"][T,j]
            domain.ge(f"edge_host{i}_spectator_{T}_{j}",spectator)
            columns[T,j] = spectator+host["family"][T]*(t if not j else host["rank"])
            domain.ge(f"edge_host{i}_edge_population_{T}_{j}",columns[T,j])
        edge_columns.append(columns)
    # The component measure copies only common population observables, not
    # independent copies of the marked host choices of this same face.
    domain.star_common += tuple(f"white_population_{T}_{j}" for T,j in white_columns)
    domain.edge_host_data = dict(white=white_columns,hosts=hosts,edge=edge_columns)


def self_checks():
    """All roots of complete inputs plus symbolic classical inverse laws."""
    from functools import reduce
    n = 101; t = F(1,n); checked = 0
    # Deletion identities are checked at every marked root, including all
    # fork/short-root boundaries; ranks test formulas, not exhaustiveness.
    exceptional = exceptional_roots()
    for kind,ranks in (("A",(1,2,3,8,31)),("B",(2,3,4,9,32)),("D",(4,5,6,10,33))):
        for rank in ranks:
            G = data.rational_root_Gram(kind,rank); inverse = white.inverse(G)
            graph = counts.labelled_graph(kind,rank)
            for node in range(rank):
                p = node+1
                if kind == "A": case="A"; x=p*t
                elif kind == "B": case="B_short" if p==rank else "B_near" if p==rank-1 else "B_tail"; x=p*t
                else:
                    case="D_leaf" if p>=rank-1 else "D_fork" if p==rank-2 else "D_near" if p==rank-3 else "D_tail"
                    x=t if case=="D_leaf" else p*t
                d = white.dimension_domain("root_fixture",n)
                d.variable("star_scale0");d.variable("star_degree0")
                host = rooted(d,0,exceptional)
                v = 1/(G[node][node]*inverse[node][node])
                # Costs use a^2 only; omit equations/blocks involving a,
                # which need not be rational on these fixtures.
                point={name:F(0) for name in d.variables}
                point.update(t=t,edge_host0_rank=rank*t,edge_host0_position=x,
                    edge_host0_projection_square=v,edge_host0_left=F(x>t),
                    edge_host0_right=F(kind=="A" and p<rank),star_degree0=F(len(graph[node])))
                point["edge_host0_"+case]=F(1)
                def evaluate(q):
                    return sum((c*reduce(lambda a,b:a*point[b],m,F(1)) for m,c in q.items()),F(0))
                for name,q in d.equations:
                    if "star_scale0" not in q.support and evaluate(q): raise ValueError("Root identity failed: "+str((kind,rank,node,name)))
                for name,q in d.inequalities:
                    if "star_scale0" not in q.support and evaluate(q)<0: raise ValueError("Root inequality failed: "+str((kind,rank,node,name)))
                literal=deleted_at(graph,node)
                if any(evaluate(q)!=t*literal[key] for key,q in host["deletion"].items()): raise ValueError("Root deletion failed")
                checked += 1
    exceptional_checked=0
    for name,entry in exceptional.items():
        T,rank=entry["family"],entry["rank"]
        for node in entry["nodes"]:
            graph=counts.labelled_graph(T[0],rank)
            if entry["deletion"]!=deleted_at(graph,node): raise ValueError("An exceptional root was incorrectly grouped")
            if not T.startswith("H"):
                G=data.rational_root_Gram(T[0],rank)
                actual=G[node][node]*white.inverse(G)[node][node]
                if actual!=entry["inverse_lower"] or actual!=entry["inverse_upper"]: raise ValueError("Wrong normalized root inverse")
            exceptional_checked+=1
    if exceptional_checked!=32: raise ValueError("An exceptional node was omitted")
    return dict(complete_exceptional_root_orbits=len(exceptional),
        every_exceptional_node_covered=True,exceptional_nodes_checked=exceptional_checked,
        classical_marked_root_regressions=checked,selected_diagram_lists_used=False,
        regression_fixtures_do_not_replace_classical_inverse_and_deletion_proofs=True)
