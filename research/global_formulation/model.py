"""A global tail outer relaxation; all feature and incidence ports are bound."""
from fractions import Fraction as F

from . import conic, data, faces, incidence, vertices
from .moments import Program, dimension_domain
from .polynomial import Poly, as_poly, monomials, var

GROUPS = {
    "neutral": (0,),
    "finite_connected": (1, *range(4, 19)),
    "ideal": (2,),
    "finite_cross": (3, *range(19, 26)),
}
REPRESENTATIVE = dict(neutral=0, finite_connected=1, ideal=2, finite_cross=3)


def combine(*terms):
    return conic.plus(*terms)


def scaled(coefficient, form):
    return conic.times(coefficient, form)


def build(threshold, parity, tables, upper_dimension=None):
    if threshold < 40 or parity not in ("odd", "even"):
        raise ValueError("This initial global formulation is specified for N>=40")
    if upper_dimension not in (None, 57):
        raise ValueError("The retained upper-dimension anchor available here is57")
    endpoint = threshold if threshold % 2 == (1 if parity == "odd" else 0) else threshold+1
    p, t, e = Program(), var("t"), var("e")
    def add_measure(name, domain, mass_scale=1):
        # Restrict domains BEFORE registering any moments or marginal rows.
        # Late substitution would leave the earlier columns disconnected
        # from the fixed dimension/exceptional-indicator quotient.
        if parity == "odd":
            domain.fixed["e"] = F(0)
        else:
            domain.ge("even_integer_gap_after58", (1-e)*(1-60*t))
        if upper_dimension is not None:
            domain.ge("retained_dimension58_tail_anchor", t-F(1, upper_dimension))
            if endpoint == upper_dimension:
                domain.fixed["t"] = F(1, upper_dimension)
                domain.fixed["e"] = F(int(upper_dimension == 58))
        return p.measure(name, domain, mass_scale)

    multipliers = (as_poly(1), t, e)
    finite_d, finite_features, finite_ports = vertices.finite_domain(endpoint, parity, tables)
    ideal_d, ideal_features, ideal_ports = vertices.ideal_domain(endpoint, tables)
    finite = add_measure("finite_vertices", finite_d, endpoint**2)
    ideal = add_measure("ideal_vertices", ideal_d, endpoint)
    face = {s: add_measure(f"faces{s}", faces.face_domain(s, endpoint)) for s in (3, 4)}
    remainder = add_measure("residual_corners", faces.corner_domain(endpoint, angle=False))
    five = add_measure("faces_at_least5", dimension_domain("large_faces", endpoint))
    negative_excess = add_measure("large_face_excess", dimension_domain("large_face_excess", endpoint))
    atoms = data.atoms()
    channels = {}
    for s in (3, 4, 0):
        common = face[s] if s else remainder
        for i in range(s if s else 1):
            tag = f"face{s}_corner{i}" if s else "residual_corner"
            atom_measures = {a: add_measure(tag+f"_atom{a}", dimension_domain("atom", endpoint)) for a in range(26)}
            channels[s, i] = atom_measures
            groups = {}
            for group, members in GROUPS.items():
                domain = faces.corner_domain(endpoint, REPRESENTATIVE[group], angle=bool(s))
                groups[group] = add_measure(tag+"_"+group, domain)
                for j, multiplier in enumerate(monomials(("t", "e"), 2)):
                    p.eq(tag+"_atom_group_"+group+f"_{j}", combine(
                        groups[group].linear(multiplier),
                        *(scaled(-1, atom_measures[a].linear(multiplier)) for a in members)))
                if s and group == "finite_connected":
                    for j, multiplier in enumerate(multipliers):
                        p.ge(tag+f"_weighted_inverse_Gram_{j}", combine(
                            groups[group].linear(var("c")*multiplier),
                            *(scaled(-F(1, data.SAME_DISTANCE[a-4]+1) if a >= 4 else -F(1, 17),
                                atom_measures[a].linear(multiplier)) for a in members)))
            names = ("t", "e", "c", "eta", "right", "chi") if s else ("t", "e", "eta", "right", "chi")
            rename = {n: n+str(i) for n in names if n not in ("t", "e")} if s else {}
            for j, multiplier in enumerate(monomials(names, 2)):
                p.eq(tag+f"_shared_corner_marginal{j}", combine(
                    common.linear(multiplier.rename(rename)),
                    *(scaled(-1, groups[g].linear(multiplier)) for g in GROUPS)))
    # Every feature, including the full old 26-vector, uses exactly these
    # same occupation measures on the vertex and face sides.
    for name in vertices.FEATURES:
        for k, multiplier in enumerate(multipliers):
            terms = [finite.linear(finite_features[name]*multiplier), ideal.linear(ideal_features[name]*multiplier)]
            if name.startswith("phi"):
                j = int(name[3:])
                for atom_measures in channels.values():
                    terms.extend(scaled(-atoms[a][j], atom_measures[a].linear(multiplier))
                        for a in range(26) if atoms[a][j])
            else:
                key = {"chi58": "chi", "ideal_corner": "eta", "right_corner": "right"}[name]
                terms.append(scaled(-1, remainder.linear(var(key)*multiplier)))
                for s, measure in face.items():
                    terms.append(scaled(-1, measure.linear(sum((var(key+str(i)) for i in range(s)), as_poly(0))*multiplier)))
            p.eq(f"feature_marginal_{name}_potential{k}", combine(*terms))
    for k, multiplier in enumerate(multipliers):
        p.eq(f"large_face_corner_incidence{k}", combine(remainder.linear(multiplier),
            scaled(-5, five.linear(multiplier)), scaled(-1, negative_excess.linear(multiplier))))
    scales = {"f1": F(1, endpoint), "f2": 1, "f3": endpoint, "B": endpoint,
        "E3": endpoint, "K3": endpoint}
    counts = {name: add_measure("count_"+name, dimension_domain("count", endpoint), scales.get(name, endpoint**2))
        for name in incidence.VARIABLES if name not in incidence.PORTS}
    for k, multiplier in enumerate(multipliers):
        p.eq(f"all_two_face_counts{k}", combine(counts["f2"].linear(multiplier),
            scaled(-1, face[3].linear(multiplier)), scaled(-1, face[4].linear(multiplier)), scaled(-1, five.linear(multiplier))))
        p.eq(f"all_edge_endpoints{k}", combine(scaled(2, counts["f1"].linear(multiplier)),
            scaled(-1, finite.linear(finite_ports["edge_endpoints"]*multiplier)),
            scaled(-1, ideal.linear(ideal_ports["edge_endpoints"]*multiplier))))
    for row in incidence.polynomial_rows(parity):
        degree = max(len(c)-1 for c in row["coefficients"].values())
        terms = []
        for name, coefficients in row["coefficients"].items():
            coefficient = sum((c*t**(degree-i) for i, c in enumerate(coefficients)), as_poly(0))
            if name in incidence.PORTS:
                terms += [finite.linear(coefficient*finite_ports[name]), ideal.linear(coefficient*ideal_ports[name])]
            else:
                terms.append(counts[name].linear(coefficient))
        (p.eq if row["sense"] == "=" else p.ge)("global_"+row["name"], combine(*terms))
    p.eq("positive_incidence_excess", combine(conic.form(-1),
        scaled(2, face[3].linear(1)), face[4].linear(1), scaled(-1, negative_excess.linear(1)),
        scaled(-1, finite.linear(finite_ports["allowance"])), scaled(-1, ideal.linear(ideal_ports["allowance"]))))
    model = p.compile(dict(
        problem_role="necessary_data_primal_for_the_whole_dimension_tail",
        threshold=threshold, parity=parity, parity_endpoint=endpoint,
        upper_dimension=upper_dimension, existing_dimension58_tail_used=upper_dimension == 57,
        intrinsic_feature_coordinates=26, total_feature_coordinates=len(vertices.FEATURES),
        feature_potential_basis=("1", "1/n", "[n=58]"),
        all_vertex_feature_and_Prokhorov_ports_use_same_measures=True,
        face_domain="one universal Lorentz Gram schema for all triangles and quadrilaterals",
        spherical_and_affine_rank_parameters_unbounded=True,
        selected_vertex_or_face_diagram_lists_used=False,
        finite_SDP_exact_realizability_claimed=False, finite_sparse_outer_moments=True,
        global_bound_changed=False, independent_source_review_complete=False,
        actual_polytope_lift="notes/lp/global_formulation.md"))
    return model
