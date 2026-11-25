"""Common supporting forests, component responses and closed-face marginals.

Every connected spherical family is present. A/B/D ranks and marked chain
positions are symbolic; no forest or marked-star list is generated. The
finite moment formulation remains an outer relaxation, including rational
envelopes for H3/H4 and necessary gaps for the unbounded Coxeter labels.
"""
from fractions import Fraction as F
from itertools import combinations
from itertools import product

from . import conic, data
from .moments import dimension_domain
from .polynomial import as_poly, var

FAMILIES = ("A", "B", "D", "E6", "E7", "E8", "F4", "H3", "H4", "I2")
FIXED = {"E6": 6, "E7": 7, "E8": 8, "F4": 4, "H3": 3, "H4": 4, "I2": 2}


def inverse(matrix):
    n = len(matrix)
    A = [list(row)+[F(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = A[j][j]
        if not pivot: raise ValueError("Zero pivot in positive component Gram matrix")
        A[j] = [x/pivot for x in A[j]]
        for i in range(n):
            if i != j:
                c = A[i][j]
                A[i] = [a-c*b for a, b in zip(A[i], A[j], strict=True)]
    return [row[n:] for row in A]


class Golden:
    """Exact a+b*phi, with phi^2=phi+1; only fixed H inverses use it."""
    def __init__(self, a=0, b=0): self.a, self.b = F(a), F(b)
    def __bool__(self): return bool(self.a or self.b)
    def __add__(self, other):
        o = other if isinstance(other, Golden) else Golden(other)
        return Golden(self.a+o.a, self.b+o.b)
    __radd__ = __add__
    def __neg__(self): return Golden(-self.a, -self.b)
    def __sub__(self, other): return self+-other if isinstance(other, Golden) else self+(-F(other))
    def __rsub__(self, other): return -self+other
    def __mul__(self, other):
        o = other if isinstance(other, Golden) else Golden(other)
        return Golden(self.a*o.a+self.b*o.b, self.a*o.b+self.b*o.a+self.b*o.b)
    __rmul__ = __mul__
    def __truediv__(self, other):
        o = other if isinstance(other, Golden) else Golden(other)
        norm = o.a*o.a+o.a*o.b-o.b*o.b
        if not norm: raise ZeroDivisionError("Zero golden-field divisor")
        return self*Golden((o.a+o.b)/norm, -o.b/norm)
    def __rtruediv__(self, other): return Golden(other)/self


def fixed_inverse_bounds(family):
    rank = FIXED[family]
    if family.startswith("H"):
        phi = Golden(0, 1)
        G = [[Golden(2*int(i == j)-int(abs(i-j) == 1)) for j in range(rank)] for i in range(rank)]
        G[0][1] = G[1][0] = -phi
        K = inverse(G)
        lo, hi = F(1618033988749894, 10**15), F(1618033988749895, 10**15)
        if not (lo > 1 and lo*lo-lo-1 < 0 < hi*hi-hi-1):
            raise ValueError("Invalid exact golden-ratio bracket")
        lower = [[q.a+q.b*(lo if q.b >= 0 else hi) for q in row] for row in K]
        upper = [[q.a+q.b*(hi if q.b >= 0 else lo) for q in row] for row in K]
    else:
        K = inverse(data.rational_root_Gram(family[0], rank))
        lower, upper = K, K
    return lower, upper


def corner_pairs(sides):
    return [(i, (i+1) % sides) for i in range(sides)] if sides in (3, 4) else [(0, 1)]


def corner_name(stem, i, sides): return stem+str(i) if sides in (3, 4) else stem
def response(i, j): return var(f"star_response_{min(i,j)}_{max(i,j)}")
def original(i, j): return var(f"star_original_{min(i,j)}_{max(i,j)}")


def corner_patterns(sides):
    for states in product(range(3),repeat=len(corner_pairs(sides))):
        p = as_poly(1)
        for i,state in enumerate(states):
            eta,right = (var(corner_name(n,i,sides)) for n in ("eta","right"))
            p *= (eta,right,1-eta-right)[state]
        yield "".join(map(str,states)),p


def component_fraction(domain):
    return sum((var(name) for name in domain.integer_counts), as_poly(0))


def vertex_statistics(finite_domain, ideal_domain):
    """Exact forest/link deletion sums in the existing whole-vertex measures."""
    t = var("t")
    k = component_fraction(finite_domain)
    finite = dict(white=t*t*(1-2*t)*((1-3*t)*k/2+t),
        finite_components=t*t*(1-t)*k/2,
        edge_components=t**3*((1-2*t)*k+t))
    h, g = as_poly(0), as_poly(0)
    for family, rank in data.affine_types():
        u = var(f"u_{family}{rank}")
        cycle = family == "A" and rank >= 2
        edges = rank+1 if cycle else rank
        h += F(edges*(rank-1), rank+1)*u
        g += (1 if cycle else F(2*rank, rank+1))*u
    for family in ("A", "B", "C", "D"):
        U, R, V = (var(f"tail_{family}_{suffix}") for suffix in ("count", "rank", "reciprocal"))
        h += R-U if family == "A" else R-2*U+2*V
        g += U if family == "A" else 2*U-2*V
    ideal = dict(white=t*((1-t)*(1-2*t)-(1-3*t)*h)/2,
        finite_components=as_poly(0), edge_components=t*t*g)
    for domain, statistics in ((finite_domain, finite), (ideal_domain, ideal)):
        for name, q in statistics.items(): domain.ge("supporting_forest_"+name, q)
    return finite, ideal


def extend_face(domain, sides, joint=False):
    t = var("t")
    k = domain.variable("star_white_components")
    domain.ge("white_component_minimum", k-t)
    domain.ge("white_component_rank_capacity", 1-2*t-k)
    a = [domain.variable(f"star_scale{i}") for i in range(sides)]
    degree = [domain.variable(f"star_degree{i}") for i in range(sides)]
    for i in range(sides):
        domain.upper_bounds[f"star_degree{i}"] = F(3)
        z = domain.variable(f"star_scale_inverse{i}", bounded=False)
        domain.eq(f"strict_edge_projection{i}", a[i]*z-1)
        domain.eq(f"integer_side_degree{i}", degree[i]*(degree[i]-1)*(degree[i]-2)*(degree[i]-3))
        domain.ge(f"edge_component_minimum{i}", k+t*(1-degree[i])-t)
        domain.ge(f"edge_component_maximum{i}", 1-t-k-t*(1-degree[i]))
        domain.basis(f"projection_pair{i}", (f"star_scale{i}", f"star_scale_inverse{i}"), 1)
    for i in range(sides):
        for j in range(i, sides): domain.variable(f"star_response_{i}_{j}")
        domain.eq(f"white_response_diagonal{i}", response(i,i)+a[i]*a[i]-1)
    domain.matrix("total_white_response", [[response(i,j) for j in range(sides)] for i in range(sides)])
    adjacent = {tuple(sorted(pair)) for pair in corner_pairs(sides)}
    for i, j in combinations(range(sides), 2):
        domain.variable(f"star_original_{i}_{j}", bounded=(i,j) in adjacent)
    common = ["star_white_components", *(f"star_scale{i}" for i in range(sides)),
        *(f"star_degree{i}" for i in range(sides))]
    for corner, (i, j) in enumerate(corner_pairs(sides)):
        eta = var(corner_name("eta", corner, sides))
        right = var(corner_name("right", corner, sides))
        c = var(corner_name("c", corner, sides))
        e = domain.variable(f"star_edge{corner}", boolean=True)
        w = domain.variable(f"star_common{corner}")
        domain.upper_bounds[f"star_common{corner}"] = F(2)
        domain.ge(f"original_edge_lower{corner}", original(i,j)-e/2)
        domain.ge(f"original_edge_upper{corner}", e-original(i,j))
        domain.eq(f"Schur_adjacent{corner}", original(i,j)+response(i,j)-c*a[i]*a[j])
        domain.eq(f"integer_common{corner}", w*(w-1)*(w-2))
        domain.ge(f"corner_cycle_capacity{corner}", 1+eta-e-w)
        domain.eq(f"finite_corner_connectivity{corner}", (1-eta)*(e+w-1+right))
        domain.ge(f"corner_connectivity{corner}", e+w-1+right)
        count = k+t*(2-degree[i]-degree[j]-e)
        domain.ge(f"finite_corner_component_minimum{corner}", (1-eta)*(count-t))
        domain.ge(f"finite_corner_component_maximum{corner}", (1-eta)*(1-count))
        common += [corner_name(n, corner, sides) for n in ("c", "eta", "right")]
        common += [f"star_edge{corner}", f"star_common{corner}"]
    if sides == 4:
        for opposite, (i,j) in enumerate(((0,2),(1,3))):
            domain.eq(f"Schur_opposite{opposite}",
                original(i,j)+response(i,j)-var(f"opposite{opposite}")*a[i]*a[j])
            if joint:
                z = domain.variable(f"opposite_strict_inverse{opposite}",bounded=False)
                gap = var(f"opposite{opposite}")-1
                domain.eq(f"strict_opposite_incidence{opposite}",gap*z-1)
                domain.matrix(f"opposite_gap_localizer{opposite}",[[gap,gap*z],[gap*z,gap*z*z]])
                domain.basis(f"opposite_strict_witness{opposite}",(f"opposite_strict_inverse{opposite}",),2)
    domain.basis("shared_supporting_state", common, 1)
    domain.basis("white_response_and_scales", [f"star_response_{i}_{j}" for i in range(sides)
        for j in range(i,sides)]+[f"star_scale{i}" for i in range(sides)], 1)
    domain.star_common = tuple(common)
    if joint:
        spatial = next((M for name,M in domain.matrices if name == "common_spatial_Gram"),None)
        for label,p in corner_patterns(sides):
            domain.ge("closed_pattern_mass_"+label,p)
            domain.matrix("closed_pattern_white_"+label,
                [[p*response(i,j) for j in range(sides)] for i in range(sides)])
            if spatial is not None:
                domain.matrix("closed_pattern_spatial_"+label,[[p*q for q in row] for row in spatial])


def component_domain(sides, threshold, common_domain, joint):
    d, t = dimension_domain(f"supporting_component{sides}", threshold), var("t")
    family = {T: d.variable("white_type_"+T, boolean=True) for T in FAMILIES}
    d.eq("one_connected_family", sum(family.values(), as_poly(0))-1)
    for T, U in combinations(FAMILIES, 2): d.eq("disjoint_families_"+T+"_"+U, family[T]*family[U])
    R = d.variable("white_rank")
    for T, minimum in (("A",1),("B",2),("D",4)):
        d.ge(T+"_minimum_rank", family[T]*(R-minimum*t))
    for T, rank in FIXED.items(): d.eq(T+"_rank", family[T]*(R-rank*t))
    d.ge("component_rank_capacity", 1-2*t-R)
    fixed = sum((family[T] for T in FIXED), as_poly(0))
    classical = 1-fixed
    g, b, x, l1, l2, v1, v2, z, edge_square = [], [], [], [], [], [], [], [], []
    for i in range(sides):
        gi = d.variable(f"attach{i}", boolean=True)
        bi = d.variable(f"amplitude{i}"); d.upper_bounds[f"amplitude{i}"] = F(2)
        xi = d.variable(f"position{i}")
        li, lj = (d.variable(f"leaf{k}_{i}", boolean=True) for k in (1,2))
        vi, vj = (d.variable(f"leaf_amplitude{k}_{i}") for k in (1,2))
        d.upper_bounds[f"leaf_amplitude1_{i}"] = d.upper_bounds[f"leaf_amplitude2_{i}"] = F(2)
        chain = d.variable(f"chain{i}", boolean=True)
        ei = d.variable(f"edge_square{i}")
        zi = [d.variable(f"node_amplitude{i}_{p}") for p in range(8)]
        for p in range(8): d.upper_bounds[f"node_amplitude{i}_{p}"] = F(2)
        d.eq(f"amplitude_support{i}", bi*(1-gi))
        d.eq(f"classical_mark{i}", chain+li+lj-classical*gi)
        d.eq(f"first_leaf_family{i}", li*(1-family["B"]-family["D"]))
        d.eq(f"second_leaf_family{i}", lj*(1-family["D"]))
        d.eq(f"first_leaf_amplitude{i}", vi-bi*li)
        d.eq(f"second_leaf_amplitude{i}", vj-bi*lj)
        d.ge(f"chain_position_lower{i}", xi-t*chain)
        d.ge(f"chain_position_upper{i}", R*chain-xi)
        d.ge(f"B_chain_end{i}", family["B"]*(R-t-xi))
        d.ge(f"D_chain_end{i}", family["D"]*(R-2*t-xi))
        d.eq(f"fixed_node_amplitude{i}", sum(zi, as_poly(0))-fixed*bi)
        for p,q in combinations(range(8),2): d.eq(f"one_mark{i}_{p}_{q}", zi[p]*zi[q])
        for T, rank in FIXED.items():
            for p in range(rank,8): d.eq(f"{T}_node_range{i}_{p}", family[T]*zi[p])
        for T in ("A","D"):
            d.eq(f"{T}_edge_square{i}", family[T]*(2*ei-bi*bi))
        d.eq(f"B_edge_square{i}", family["B"]*(2*ei-(1+li)*bi*bi))
        for T, rank in FIXED.items():
            lengths = [2,2,1,1] if T == "F4" else [2]*rank
            d.eq(f"{T}_edge_square{i}", family[T]*(ei-sum((zi[p]*zi[p]/lengths[p]
                for p in range(rank)), as_poly(0))))
        d.ge(f"attachment_label_lower{i}", ei-gi/4)
        d.ge(f"attachment_label_upper{i}", gi-ei)
        d.ge(f"attachment_label_first_gap{i}", (ei-gi/4)*(ei-gi/2))
        d.ge(f"attachment_label_second_gap{i}", (ei-gi/2)*(ei-5*gi/8))
        g.append(gi); b.append(bi); x.append(xi); l1.append(li); l2.append(lj)
        v1.append(vi); v2.append(vj); z.append(zi); edge_square.append(ei)
        d.basis(f"fixed_mark{i}", [f"node_amplitude{i}_{p}" for p in range(8)]+
            [f"amplitude{i}", *("white_type_"+T for T in FIXED)], 1)
    minimum = {}
    for i,j in combinations(range(sides),2):
        m = d.variable(f"minimum{i}_{j}")
        order = d.variable(f"order{i}_{j}", boolean=True)
        d.eq(f"marked_minimum{i}_{j}", m-order*x[i]-(1-order)*x[j])
        d.ge(f"marked_order_left{i}_{j}", order*(x[j]-x[i]))
        d.ge(f"marked_order_right{i}_{j}", (1-order)*(x[i]-x[j]))
        minimum[i,j] = m
        d.basis(f"marked_pair{i}_{j}", (f"position{i}",f"position{j}",f"minimum{i}_{j}",f"order{i}_{j}"),1)
    Q = [[None]*sides for _ in range(sides)]
    cI = d.variable("dihedral_parameter")
    d.ge("I2_parameter_lower", cI-F(1,2))
    d.ge("I2_canonical_label_lower",family["I2"]*(cI*cI-F(5,8)))
    fixed_bounds = {T: fixed_inverse_bounds(T) for T in FIXED if T != "I2"}
    for i in range(sides):
        for j in range(i,sides):
            q = d.variable(f"component_response_{i}_{j}"); Q[i][j] = Q[j][i] = q
            d.eq(f"response_support_left{i}_{j}", q*(1-g[i]))
            d.eq(f"response_support_right{i}_{j}", q*(1-g[j]))
            m = x[i] if i == j else minimum[i,j]
            d.eq(f"A_inverse_response{i}_{j}", family["A"]*(t*(R+t)*q-
                b[i]*b[j]*m*(R+t-x[i]-x[j]+m)))
            ui, uj = b[i]-v1[i], b[j]-v1[j]
            d.eq(f"B_inverse_response{i}_{j}", family["B"]*(t*q-b[i]*b[j]*m-
                ui*v1[j]*x[i]-v1[i]*uj*x[j]-v1[i]*v1[j]*R))
            ui, uj = b[i]-v1[i]-v2[i], b[j]-v1[j]-v2[j]
            vi, vj = v1[i]+v2[i], v1[j]+v2[j]
            d.eq(f"D_inverse_response{i}_{j}", family["D"]*(4*t*q-4*b[i]*b[j]*m-
                2*ui*vj*x[i]-2*vi*uj*x[j]-R*(v1[i]*v1[j]+v2[i]*v2[j])-
                (R-2*t)*(v1[i]*v2[j]+v2[i]*v1[j])))
            d.eq(f"I2_inverse_response{i}_{j}", family["I2"]*(2*(1-cI*cI)*q-
                z[i][0]*z[j][0]-z[i][1]*z[j][1]-cI*(z[i][0]*z[j][1]+z[i][1]*z[j][0])))
            for T,(lower,upper) in fixed_bounds.items():
                rank = FIXED[T]
                lower_q = sum((lower[p][r]*z[i][p]*z[j][r] for p in range(rank) for r in range(rank)), as_poly(0))
                upper_q = sum((upper[p][r]*z[i][p]*z[j][r] for p in range(rank) for r in range(rank)), as_poly(0))
                if lower == upper:
                    d.eq(f"{T}_inverse_response{i}_{j}", family[T]*(q-lower_q))
                else:
                    d.ge(f"{T}_inverse_response_lower{i}_{j}", family[T]*(q-lower_q))
                    d.ge(f"{T}_inverse_response_upper{i}_{j}", family[T]*(upper_q-q))
        d.ge(f"inverse_diagonal_cost{i}", Q[i][i]-edge_square[i])
    d.matrix("component_inverse_response", Q)
    local = [n for n in d.variables if n.startswith(("white_","attach","amplitude","position","leaf","chain","edge_square"))]
    d.basis("component_inverse_state", local+[f"component_response_{i}_{j}" for i in range(sides) for j in range(i,sides)],1)
    d.star_common = common_domain.star_common if joint else ()
    for name in d.star_common:
        d.variable(name, boolean=name in common_domain.boolean)
        d.upper_bounds[name] = common_domain.upper_bounds[name]
    if joint:
        for i in range(len(corner_pairs(sides))):
            eta,right,c = (var(corner_name(n,i,sides)) for n in ("eta","right","c"))
            d.ge(f"shared_corner_kind{i}",1-eta-right)
            d.eq(f"shared_corner_disjoint{i}",eta*right)
            d.eq(f"shared_ideal_angle{i}",eta*(1-c))
            d.eq(f"shared_right_angle{i}",right*c)
        for label,p in corner_patterns(sides):
            d.ge("closed_pattern_mass_"+label,p)
            d.matrix("closed_pattern_component_response_"+label,[[p*q for q in row] for row in Q])
        ports = ["white_rank", *("white_type_"+T for T in FAMILIES), *(f"attach{i}" for i in range(sides)),
            *(f"component_response_{i}_{j}" for i in range(sides) for j in range(i,sides))]
        d.basis("joint_face_attachment_state", [*d.star_common,*ports],1)
    return d


def component_couplings(program,face,profile,sides,joint):
    t = var("t")
    tests = [as_poly(1)]
    if joint:
        tests += [var(n) for n in face.domain.star_common]
        tests += [p for _,p in corner_patterns(sides)]
    tests = list({conic.canonical(p):p for p in tests}.values())
    for k,test in enumerate(tests):
        def match(name,component,whole):
            program.eq(f"support{sides}_{name}_joint{k}",conic.plus(profile.linear(test*component),
                conic.times(-1,face.linear(test*whole))))
        match("component_count",t,var("star_white_components"))
        match("component_rank",var("white_rank"),1-2*t)
        for i in range(sides):
            match(f"side_degree{i}",var(f"attach{i}"),var(f"star_degree{i}"))
            for j in range(i,sides):
                match(f"inverse_response{i}_{j}",var(f"component_response_{i}_{j}"),response(i,j))
        for corner,(i,j) in enumerate(corner_pairs(sides)):
            match(f"common_component{corner}",var(f"attach{i}")*var(f"attach{j}"),var(f"star_common{corner}"))


def couple(program, finite, ideal, faces, residual, edge_count, add_measure, statistics, level,family_statistics=None):
    t = var("t")
    all_faces = {**faces, 2: residual}
    profiles = {}
    for sides,face in all_faces.items():
        profile = add_measure(f"supporting_components{sides}",
            component_domain(sides,face.domain.t_scale,face.domain,joint=level >= 4))
        profiles[sides] = profile
        component_couplings(program,face,profile,sides,joint=level >= 4)
    if family_statistics is not None:
        for (T,j),column in family_statistics[0].items():
            terms = [finite.linear(column),ideal.linear(family_statistics[1][T,j])]
            ideal_terms = [ideal.linear(family_statistics[1][T,j])]
            for sides,profile in profiles.items():
                statistic = var("white_type_"+T)*(t if not j else var("white_rank"))
                terms.append(conic.times(-len(corner_pairs(sides)),profile.linear(statistic)))
                if level >= 4:
                    ideal_corners = sum((var(corner_name("eta",c,sides))
                        for c in range(len(corner_pairs(sides)))),as_poly(0))
                    ideal_terms.append(conic.times(-1,profile.linear(ideal_corners*statistic)))
            program.eq(f"global_white_family_{T}_statistic{j}",conic.plus(*terms))
            if level >= 4:
                # The same family population must be a deletion marginal
                # of complete classified cusp factors at ideal corners.
                program.eq(f"ideal_white_family_{T}_statistic{j}",conic.plus(*ideal_terms))
    for name in ("white", "finite_components"):
        terms = [finite.linear(statistics[0][name]), ideal.linear(statistics[1][name])]
        for sides, face in all_faces.items():
            k = var("star_white_components")
            if name == "white":
                value = len(corner_pairs(sides))*k
            else:
                value = sum(((1-var(corner_name("eta",c,sides)))*(k+t*(2-var(f"star_degree{i}")-
                    var(f"star_degree{j}")-var(f"star_edge{c}")))
                    for c,(i,j) in enumerate(corner_pairs(sides))), as_poly(0))
            terms.append(conic.times(-1,face.linear(value)))
        program.eq("global_supporting_"+name, conic.plus(*terms))
    if level >= 4:
        edge_domain = dimension_domain("shared_edge_support", finite.domain.t_scale)
        k = edge_domain.variable("edge_components")
        edge_domain.ge("edge_component_minimum", k-t)
        edge_domain.ge("edge_component_maximum", 1-t-k)
        edge_domain.basis("edge_support",("t","edge_components"),2)
        edge = add_measure("shared_edge_support", edge_domain)
        program.eq("shared_edge_count", conic.plus(edge.linear(1),conic.times(-1,edge_count.linear(1))))
        program.eq("shared_edge_endpoint_components", conic.plus(conic.times(2,edge.linear(k)),
            conic.times(-1,finite.linear(statistics[0]["edge_components"])),
            conic.times(-1,ideal.linear(statistics[1]["edge_components"]))))
        terms = [conic.times(-1,edge.linear((1-t)*k))]
        for sides,face in all_faces.items():
            # Residual corners select their outgoing side, once per corner.
            indices = range(sides) if sides in (3,4) else (0,)
            value = sum((var("star_white_components")+t*(1-var(f"star_degree{i}")) for i in indices),as_poly(0))
            terms.append(face.linear(t*value))
        program.eq("closed_face_shared_edge_components",conic.plus(*terms))
    return profiles


def local_star_fixture():
    """An honest rational moment lift of a shared nonzero-response local star.

    Its supporting forest is A3 plus n-5 spectator A1 components. It is
    a local Gram fixture, not a finite-volume polytope existence claim.
    """
    from . import faces, strengthening
    from .moments import Program
    n = 57
    d = faces.face_domain(3,n)
    extend_face(d,3,joint=True)
    d.fixed.update(t=F(1,n),e=F(0)); d.equation_closure = True
    strengthening.enrich(d,2)
    p = Program(); face = p.measure("local_face",d)
    cd = component_domain(3,n,d,True)
    cd.fixed.update(t=F(1,n),e=F(0)); cd.equation_closure = True
    profile = p.measure("local_components",cd)
    component_couplings(p,face,profile,3,True)
    problem = p.compile(dict(local_fixture_only=True))
    common = dict(star_white_components=F(n-4,n))
    for i in range(3):
        common.update({f"c{i}":F(1),f"eta{i}":F(1),f"right{i}":F(0),
            f"star_scale{i}":F(1,2) if i < 2 else F(1),
            f"star_degree{i}":F(int(i < 2)),f"star_edge{i}":F(int(i > 0)),
            f"star_common{i}":F(int(i == 0))})
    face_values = dict(common)
    for i,h in enumerate((F(1,2),F(1,2),F(3,4))):
        face_values.update({f"h{i}":h,f"inverse_h{i}":1/h,
            f"star_scale_inverse{i}":1/common[f"star_scale{i}"]})
    Q = [[F(3,4),F(1,4),F(0)],[F(1,4),F(3,4),F(0)],[F(0),F(0),F(0)]]
    for i in range(3):
        for j in range(i,3):
            face_values[f"star_response_{i}_{j}"] = Q[i][j]
            if i != j: face_values[f"star_original_{i}_{j}"] = F(0) if j == 1 else F(1,2)
    active = dict(common,white_type_A=F(1),white_rank=F(3,n),dihedral_parameter=F(1,2))
    for i,position in enumerate((1,3,0)):
        active.update({f"attach{i}":F(int(i < 2)),f"amplitude{i}":F(int(i < 2)),
            f"chain{i}":F(int(i < 2)),f"position{i}":F(position,n),
            f"edge_square{i}":F(1,2) if i < 2 else F(0)})
    active.update(minimum0_1=F(1,n),order0_1=F(1))
    for i in range(3):
        for j in range(i,3): active[f"component_response_{i}_{j}"] = Q[i][j]
    spectator = dict(common,white_type_A=F(1),white_rank=F(1,n),dihedral_parameter=F(1,2))
    def value(monomial,values):
        result = F(1)
        for v in monomial: result *= values.get(v,F(0))
        return result
    point = {name:value(m,face_values) for m,name in face.moments.items()}
    point.update({name:value(m,active)+(n-5)*value(m,spectator) for m,name in profile.moments.items()})
    checked = conic.check_point(problem,point)
    return dict(moment_variables=len(problem["variables"]),equalities=checked["equalities_checked"],
        cone_blocks=len(checked["blocks"]),all_three_layers_have_an_exact_local_Gram_lift=True,
        global_polytope_existence_claimed=False)


def self_checks():
    """Independent Gram inversions and direct small flag counts, not search cases."""
    def evaluate(q, values):
        total = F(0)
        for monomial,c in q.items():
            for v in monomial: c *= values.get(v,F(0))
            total += c
        return total
    from . import faces
    from .moments import Domain
    base = faces.corner_domain(57)
    extend_face(base,2)
    domain = component_domain(2,57,base,False)
    checked = 0
    for kind,rank in [(k,r) for k in ("A","B","D") for r in (4,7,13)]+[("E",6),("E",7),("E",8),("F",4),("I",2)]:
        K = inverse([[F(2),-F(8,5)],[-F(8,5),F(2)]]) if kind == "I" else inverse(data.rational_root_Gram(kind,rank))
        family = kind if kind in ("A","B","D") else kind+str(rank)
        for i,j in combinations(range(rank),2):
            values = dict(t=F(1,57),white_rank=F(rank,57),**{"white_type_"+family:F(1)})
            if kind == "I": values["dihedral_parameter"] = F(4,5)
            nodes = (i,j)
            for port,node in enumerate(nodes):
                values[f"attach{port}"] = values[f"amplitude{port}"] = F(1)
                if kind in ("E","F","I"):
                    values[f"node_amplitude{port}_{node}"] = F(1)
                elif kind == "B" and node == rank-1:
                    values[f"leaf1_{port}"] = values[f"leaf_amplitude1_{port}"] = F(1)
                elif kind == "D" and node >= rank-2:
                    leaf = node-rank+3
                    values[f"leaf{leaf}_{port}"] = values[f"leaf_amplitude{leaf}_{port}"] = F(1)
                else:
                    values[f"position{port}"] = F(node+1,57)
            values["minimum0_1"] = min(values.get("position0",F(0)),values.get("position1",F(0)))
            for a in range(2):
                for b in range(a,2): values[f"component_response_{a}_{b}"] = K[nodes[a]][nodes[b]]
            for label,q in domain.equations:
                if label.startswith(family+"_inverse_response") and evaluate(q,values):
                    raise ValueError("A parameterized component response disagrees with direct Gram inversion: "+label)
            checked += 1
    for rank in (3,4):
        lower,upper = fixed_inverse_bounds("H"+str(rank))
        if any(a <= 0 or a > b for lo,hi in zip(lower,upper,strict=True) for a,b in zip(lo,hi,strict=True)):
            raise ValueError("An H inverse enclosure failed its positive interval check")
    def components(graph,deleted=()):
        remaining = set(range(len(graph)))-set(deleted)
        result = 0
        while remaining:
            result += 1; frontier = [remaining.pop()]
            while frontier:
                for v in graph[frontier.pop()] & remaining:
                    remaining.remove(v); frontier.append(v)
        return result
    def forest(types):
        result = []
        for kind,rank in types:
            offset = len(result)
            result += [{v+offset for v in row} for row in data.graph(kind,rank)]
        return result
    cases = ((("A",4),),(('D',4),),(('A',1),('B',3),('D',4)))
    for types in cases:
        G = forest(types); n = len(G)
        fd,iden = Domain("fixture_finite"),Domain("fixture_ideal")
        fd.integer_counts = {f"u_{k}{r}":r for k,r in types}
        statistics,_ = vertex_statistics(fd,iden)
        values = dict(t=F(1,n))
        for k,r in types: values[f"u_{k}{r}"] = values.get(f"u_{k}{r}",F(0))+F(1,n)
        white = sum(components(G,pair) for pair in combinations(range(n),2))
        edge = sum(components(G,(i,)) for i in range(n))
        if evaluate(statistics["white"],values)*n**5 != white:
            raise ValueError("Finite whole-forest two-deletion marginal failed")
        if evaluate(statistics["edge_components"],values)*n**5 != edge:
            raise ValueError("Finite whole-forest one-deletion marginal failed")
    ideal_cases = ((("A",1),("A",2)),(("B",3),("C",2)),(("A",2),("D",4)))
    for types in ideal_cases:
        graphs = [data.graph(k,r,True) for k,r in types]
        n = 1+sum(r for _,r in types); V = 1
        for G in graphs: V *= len(G)
        fd,iden = Domain("fixture_finite"),Domain("fixture_ideal")
        _,statistics = vertex_statistics(fd,iden)
        values = dict(t=F(1,n))
        for k,r in types: values[f"u_{k}{r}"] = values.get(f"u_{k}{r}",F(0))+F(1,n)
        white = 0
        for edge_factor,G in enumerate(graphs):
            choices = [list(combinations(range(len(H)),2)) if i == edge_factor else
                [(j,) for j in range(len(H))] for i,H in enumerate(graphs)]
            for deletions in product(*choices):
                white += sum(components(H,deletion) for H,deletion in zip(graphs,deletions,strict=True))
        edge = sum(sum(components(H,(node,)) for H,node in zip(graphs,nodes,strict=True))
            for nodes in product(*(range(len(H)) for H in graphs)))
        if evaluate(statistics["white"],values)*V*n**3 != white:
            raise ValueError("Complete product-cusp two-face deletion marginal failed")
        if evaluate(statistics["edge_components"],values)*V*n**3 != edge:
            raise ValueError("Complete product-cusp edge deletion marginal failed")
    return dict(component_inverse_pairs_checked=checked,
        finite_forest_flag_count_fixtures=len(cases), product_cusp_flag_count_fixtures=len(ideal_cases),
        classical_formulas_are_all_rank_identities_not_fitted_tables=True,
        no_component_rank_or_label_cutoff=True)
