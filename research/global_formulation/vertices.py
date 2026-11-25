"""Whole-vertex polynomial statistics on complete component-count domains."""
from fractions import Fraction as F

from . import data
from .moments import dimension_domain
from .polynomial import as_poly, var

FEATURES = tuple(f"phi{j}" for j in range(26))+("chi58", "ideal_corner", "right_corner", "angle_square")


def counters(d, types):
    result = {}
    for kind, rank in types:
        name = f"u_{kind}{rank}"
        result[kind, rank] = d.variable(name)
        d.upper_bounds[name] = F(1, rank)
        d.ge("integer_population_"+name, result[kind, rank]*(result[kind, rank]-var("t")))
        d.basis("count_"+name, ("t", "e", name), 1)
        d.basis("scaled_"+name, ("t", name), 2)
    return result


def tail(d, kind, minimum, reciprocal=False):
    names = [f"tail_{kind}_{suffix}" for suffix in ("count", "rank", "square", "cube")]
    u, r, s, q = (d.variable(n) for n in names)
    d.upper_bounds[names[0]] = F(1, minimum)
    t = var("t")
    d.ge(kind+"_tail_minimum", r-minimum*u)
    d.ge(kind+"_integer_tail_population", u*(u-t))
    d.ge(kind+"_tail_maximum", u-t*r)
    d.ge(kind+"_square_minimum", s-minimum*t*r)
    d.ge(kind+"_square_maximum", r-s)
    d.ge(kind+"_cube_minimum", q-minimum*t*s)
    d.ge(kind+"_cube_maximum", s-q)
    d.ge(kind+"_count_rank_Cauchy", u*s-t*r*r)
    d.matrix(kind+"_rank_moments", [[r, s], [s, q]])
    d.matrix(kind+"_count_rank_moments", [[u, t*r], [t*r, t*s]])
    d.basis(kind+"_tail", ("t", "e", *names), 1)
    answer = dict(count=u, rank=r, square=s, cube=q)
    if reciprocal:
        v = d.variable(f"tail_{kind}_reciprocal")
        d.upper_bounds[f"tail_{kind}_reciprocal"] = F(1, minimum*(minimum+1))
        d.ge(kind+"_reciprocal_upper", u-(minimum+1)*v)
        d.ge(kind+"_reciprocal_lower", (1+t)*v-t*u)
        d.matrix(kind+"_reciprocal_Cauchy", [[r+u, u], [u, v]])
        answer["reciprocal"] = v
    return answer


def rank_data(d, counts, tails, total):
    t = var("t")
    rank = sum((r*u for (_, r), u in counts.items()), as_poly(0))+sum((v["rank"] for v in tails.values()), as_poly(0))
    sigma2 = t*sum((r*r*u for (_, r), u in counts.items()), as_poly(0))+sum((v["square"] for v in tails.values()), as_poly(0))
    sigma3 = t*t*sum((r*r*r*u for (_, r), u in counts.items()), as_poly(0))+sum((v["cube"] for v in tails.values()), as_poly(0))
    d.eq("whole_vertex_rank", rank-total)
    d.ge("rank_square_upper", total*total-sigma2)
    d.ge("rank_square_lower", sigma2-t*total)
    d.ge("rank_cube_lower", sigma3-t*sigma2)
    d.ge("rank_cube_upper", total*sigma2-sigma3)
    d.ge("rank_third_elementary_symmetric", total**3-3*total*sigma2+2*sigma3)
    return sigma2, sigma3


def add_feature_cone(d, columns):
    # Linear inverse of the manuscript's 26-atom matrix, applied to the
    # aggregate of ALL corners of this same vertex.
    beta = [as_poly(0)]*26
    for j in range(1, 26): beta[j] = columns[f"phi{j}"]
    beta[1] -= sum((b*columns[f"phi{j}"] for j, b in enumerate(data.SAME_BASELINES, 4)), as_poly(0))
    beta[22] -= beta[23]
    beta[0] = columns["phi0"]-sum(beta[1:], as_poly(0))
    for j, value in enumerate(beta): d.ge(f"complete_vertex_atom{j}", value)
    for name in FEATURES[26:]:
        d.ge(name+"_nonnegative", columns[name])
        d.ge(name+"_corner_capacity", columns["phi0"]-columns[name])


def finite_domain(threshold, parity, tables):
    d = dimension_domain("finite_vertex", threshold)
    t = var("t")
    u = counters(d, data.finite_types())
    tails = {k: tail(d, k, minimum) for k, (minimum, _) in data.FINITE_TAIL.items()}
    sigma2, sigma3 = rank_data(d, u, tails, as_poly(1))
    primitive = {}
    def indicator(name, value, k):
        g = d.variable("guard_"+name, boolean=True)
        d.ge(name+"_active", g*(value-k*t))
        d.ge(name+"_inactive", (1-g)*((k-1)*t-value))
        d.basis(name+"_guard", tuple(sorted(value.support | {"t", "guard_"+name})), 1)
        primitive[name] = g
        return g
    a = {}
    for rank in (1, 2, 9, 10, 11): a[rank] = indicator("a"+str(rank), u["A", rank], 1)
    small = indicator("small_a", u["A", 1]+u["A", 2]+u["A", 3], 1)
    two_a1 = indicator("two_a1", u["A", 1], 2)
    rank15 = u["A", 15]+u["D", 15]
    two15, three15 = indicator("two15", rank15, 2), indicator("three15", rank15, 3)
    guards = dict(has_a1=a[1], a1_a10=a[1]*a[10], a1_a9=a[1]*a[9],
        a9_a2=a[9]*a[2], a11_or_a9_a2=a[11]+a[9]*a[2]-a[11]*a[9]*a[2], two15=two15)
    columns = {name: as_poly(0) for name in FEATURES}
    columns["phi0"] = t*t*(1-t)/2
    ordinary = sum((tables["finite"][c]["ordinary"]*v for c, v in u.items()), as_poly(0))
    W = sum((tables["finite"][c]["W"]*v for c, v in u.items()), as_poly(0))
    J = sum((tables["finite"][c]["J"]*v for c, v in u.items()), as_poly(0))
    for kind, values in tails.items():
        U, R = values["count"], values["rank"]
        ordinary += 5*R+data.FINITE_TAIL[kind][1]*U
        W += F(29, 3)*R+(F(-176, 3)+int(kind == "D"))*U
        J += F(205, 3)*R+(F(-1696, 3)+3*int(kind == "D"))*U
    columns["phi1"] = t**3*ordinary
    columns["phi3"] = t*t*sum((value*u[left]*u[right] for (left, right), value in tables["ordinary_cross"].items()), as_poly(0))
    for j, component in enumerate(data.SAME_TYPES, 4):
        g = a[1] if j <= 13 else small if j <= 15 else three15 if j <= 17 else two15
        columns[f"phi{j}"] = t**3*u[component]*g
    for j, left, right, orbit, guard in data.CROSS:
        columns[f"phi{j}"] = orbit*t*t*u[left]*u[right]*guards[guard]
    columns["phi25"] = 2*t*t*(u["A", 2]**2-t*u["A", 2])*a[10]*two_a1
    for j, minimum in data.feature_minimum_ranks().items():
        d.ge(f"feature{j}_whole_rank_support", (1-minimum*t)*columns[f"phi{j}"])
    columns["right_corner"] = t*t*(1-sigma2)/2
    angle = d.variable("angle_square_total")
    lower = t*sum((tables["finite"][c]["angle_square_lower"]*v for c, v in u.items()), as_poly(0))
    upper = t*sum((tables["finite"][c]["angle_square_upper"]*v for c, v in u.items()), as_poly(0))
    for kind, values in tails.items():
        U, R, S = values["count"], values["rank"], values["square"]
        if kind == "A":
            square_order = S+2*t*R+t*U
            lower += (data.HARMONIC2_LOWER-F(3, 2))*square_order-t*(R+U)/2-t*U/data.FINITE_TAIL[kind][0]
            upper += (data.HARMONIC2_UPPER-F(3, 2))*square_order-t*(R+U)/2+t*U
        elif kind == "B":
            lower += (S-t*R)/4
            upper += (S-t*R)/4
        else:
            lower += (S-t*R)/4-(F(1, 2)+F(2, data.FINITE_TAIL[kind][0]))*t*U
            upper += (S-t*R)/4-t*U/2
    d.ge("complete_inverse_Gram_angle_lower", angle-lower)
    d.ge("complete_inverse_Gram_angle_upper", upper-angle)
    d.ge("angle_square_rank_capacity", (sigma2-t)/2-angle)
    columns["angle_square"] = t*t*angle
    allowance = t*t*(1-6*t+(5*t*t if parity == "odd" else 0))/8
    ports = dict(P=t*(1-t)*(1-2*t)/6, Q=as_poly(0),
        A=t*(1-3*sigma2+2*sigma3)/6, S=t*t*(W-t*J)/2,
        edge_endpoints=t**3, allowance=allowance)
    add_feature_cone(d, columns)
    for name, p in ports.items(): d.ge("whole_vertex_"+name, p)
    return d, columns, ports


def ideal_domain(threshold, tables):
    d = dimension_domain("ideal_vertex", threshold)
    t, e = var("t"), var("e")
    u = counters(d, data.affine_types())
    tails = {k: tail(d, k, data.AFFINE_MINIMUM, reciprocal=True) for k in ("A", "B", "C", "D")}
    sigma2, _ = rank_data(d, u, tails, 1-t)
    b, chi = d.variable("absolute_height"), d.variable("chi_fraction")
    d.ge("absolute_height_unimodal", F(1, 2)-b)
    d.ge("absolute_height_variance", sigma2+2*t*(1-t)-3*b*b)
    d.eq("chi_only_at58", chi*(1-e))
    q = sum((tables["affine"][c]["ordinary_per_order"]*v for c, v in u.items()), as_poly(0))
    H = sum((tables["affine"][c]["deletion_per_order"]*v for c, v in u.items()), as_poly(0))
    for kind, values in tails.items():
        U, R, V = values["count"], values["rank"], values["reciprocal"]
        q += 5*U+{"A": 0, "B": 46, "C": 41, "D": 51}[kind]*V
        beta, reciprocal = {"A": (-176, 0), "B": (-349, 1687),
            "C": (-352, 1696), "D": (-346, 1678)}[kind]
        H += (29*R+beta*U+reciprocal*V)/3
    columns = {name: as_poly(0) for name in FEATURES}
    columns["phi0"] = t*(1-t)/2
    columns["phi2"] = t*q
    columns["ideal_corner"] = columns["phi0"]
    columns["angle_square"] = columns["phi0"]
    columns["chi58"] = columns["phi0"]*chi
    ports = dict(P=(sigma2-t*(1-t))/6, Q=((1-t)**2-sigma2)/8,
        A=as_poly(0), S=t*H/2, edge_endpoints=t*t, allowance=t*(3*(1-t)-5*b)/4)
    add_feature_cone(d, columns)
    for name, p in ports.items(): d.ge("whole_vertex_"+name, p)
    return d, columns, ports
