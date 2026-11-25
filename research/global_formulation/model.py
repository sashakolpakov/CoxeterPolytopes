"""A global tail outer relaxation; all feature and incidence ports are bound."""
from fractions import Fraction as F

from . import conic, data, faces, incidence, vertices, strengthening, white_profiles, supporting_counts
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


def build(threshold, parity, tables, upper_dimension=None, level=1):
    if threshold < 40 or parity not in ("odd", "even"):
        raise ValueError("This initial global formulation is specified for N>=40")
    if upper_dimension not in (None, 57):
        raise ValueError("The retained upper-dimension anchor available here is57")
    if level not in (1,2,3,4,5): raise ValueError("Use a retained global hierarchy level")
    endpoint = threshold if threshold % 2 == (1 if parity == "odd" else 0) else threshold+1
    p, t, e = Program(), var("t"), var("e")
    def add_measure(name, domain, mass_scale=1):
        # Restrict domains BEFORE registering any moments or marginal rows.
        # Late substitution would leave the earlier columns disconnected
        # from the fixed dimension/exceptional-indicator quotient.
        if parity == "odd":
            domain.fixed["e"] = F(0)
        else:
            if endpoint >= 58:
                domain.ge("even_integer_gap_after58", (1-e)*(1-60*t))
            else:
                domain.ge("even_integer_gap_around58", (1-e)*(1-56*t)*(1-60*t))
        if upper_dimension is not None:
            domain.ge("retained_dimension58_tail_anchor", t-F(1, upper_dimension))
            if endpoint == upper_dimension:
                domain.fixed["t"] = F(1, upper_dimension)
                domain.fixed["e"] = F(int(upper_dimension == 58))
        strengthening.enrich(domain, min(level,2))
        domain.equation_closure = level >= 5
        return p.measure(name, domain, mass_scale)

    multipliers = (as_poly(1), t, e)
    finite_d, finite_features, finite_ports = vertices.finite_domain(endpoint, parity, tables)
    ideal_d, ideal_features, ideal_ports = vertices.ideal_domain(endpoint, tables, parity)
    supporting_statistics = white_profiles.vertex_statistics(finite_d,ideal_d) if level >= 3 else None
    family_statistics = supporting_counts.vertex_columns() if level >= 3 else None
    if level >= 3:
        for domain,columns in zip((finite_d,ideal_d),family_statistics,strict=True):
            for (T,j),q in columns.items(): domain.ge(f"supporting_family_{T}_statistic{j}",q)
    finite = add_measure("finite_vertices", finite_d, endpoint**2)
    ideal = add_measure("ideal_vertices", ideal_d, endpoint)
    face = {}
    for s in (3,4):
        domain = faces.face_domain(s,endpoint)
        if level >= 3: white_profiles.extend_face(domain,s,joint=level >= 4)
        face[s] = add_measure(f"faces{s}",domain)
    domain = faces.corner_domain(endpoint,angle=True)
    if level >= 3: white_profiles.extend_face(domain,2,joint=level >= 4)
    remainder = add_measure("residual_corners",domain)
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
                domain = faces.corner_domain(endpoint, REPRESENTATIVE[group], angle=True)
                groups[group] = add_measure(tag+"_"+group, domain)
                for j, multiplier in enumerate(monomials(("t", "e"), 2)):
                    p.eq(tag+"_atom_group_"+group+f"_{j}", combine(
                        groups[group].linear(multiplier),
                        *(scaled(-1, atom_measures[a].linear(multiplier)) for a in members)))
                if group == "finite_connected":
                    for j, multiplier in enumerate(multipliers):
                        p.ge(tag+f"_weighted_inverse_Gram_{j}", combine(
                            groups[group].linear(var("c")*multiplier),
                            *(scaled(-F(1, data.SAME_DISTANCE[a-4]+1) if a >= 4 else -F(1, 17),
                                atom_measures[a].linear(multiplier)) for a in members)))
            names = ("t", "e", "c", "eta", "right", "chi")
            rename = {n: n+str(i) for n in names if n not in ("t", "e")} if s else {}
            for j, multiplier in enumerate(monomials(names, 2)):
                p.eq(tag+f"_shared_corner_marginal{j}", combine(
                    common.linear(multiplier.rename(rename)),
                    *(scaled(-1, groups[g].linear(multiplier)) for g in GROUPS)))
    # The retained minimal-indefinite-core triangle proof selects two
    # distinct pure ordinary-unit corners. Its spectator-rank argument is
    # uniform; see the explicit source analysis in the specification.
    for k, multiplier in enumerate(multipliers):
        p.ge(f"uniform_two_unit_triangle_cover{k}", combine(
            scaled(-2, face[3].linear(multiplier)),
            *(channels[3, i][a].linear(multiplier) for i in range(3) for a in (1, 2))))
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
                key = {"chi58": "chi", "ideal_corner": "eta", "right_corner": "right", "angle_square": "c"}[name]
                power = 2 if name == "angle_square" else 1
                terms.append(scaled(-1, remainder.linear(var(key)**power*multiplier)))
                for s, measure in face.items():
                    terms.append(scaled(-1, measure.linear(sum((var(key+str(i))**power for i in range(s)), as_poly(0))*multiplier)))
            p.eq(f"feature_marginal_{name}_potential{k}", combine(*terms))
    for k, multiplier in enumerate(multipliers):
        p.eq(f"large_face_corner_incidence{k}", combine(remainder.linear(multiplier),
            scaled(-5, five.linear(multiplier)), scaled(-1, negative_excess.linear(multiplier))))
    scales = {"f1": F(1, endpoint), "f2": 1, "f3": endpoint, "B": endpoint,
        "E3": endpoint, "K3": endpoint}
    counts = {name: add_measure("count_"+name, dimension_domain("count", endpoint), scales.get(name, endpoint**2))
        for name in incidence.VARIABLES if name not in incidence.PORTS}
    if level >= 3:
        white_profiles.couple(p,finite,ideal,face,remainder,counts["f1"],add_measure,supporting_statistics,level,
            family_statistics=family_statistics)
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
        old_uniform_two_unit_triangle_cover_included=True,
        all_old_quadrilateral_transfer_premises_imported=False,
        face_domain="one universal Lorentz Gram schema for all triangles and quadrilaterals",
        spherical_and_affine_rank_parameters_unbounded=True,
        selected_vertex_or_face_diagram_lists_used=False,
        finite_SDP_exact_realizability_claimed=False, finite_sparse_outer_moments=True,
        global_bound_changed=False, independent_source_review_complete=False,
        actual_polytope_lift="notes/lp/global_formulation.md"))
    if level > 1:
        model["metadata"]["uniform_strengthening_level"] = level
        model["metadata"]["uniform_strengthening"] = (
            "all mixed spatial rank minors; compound/tensor and weighted Gram PSD; "
            "integer factorial, gap and rank-capacity localizers for all population counters")
    if level >= 3:
        model["metadata"]["shared_supporting_profiles"] = dict(
            complete_families=white_profiles.FAMILIES,
            unbounded_A_B_D_ranks_and_I2_labels=True,
            H_inverse_entries_use_proved_rational_envelopes=True,
            all_sides_share_each_component_response=True,
            white_deletion_and_finite_vertex_population_marginals=True,
            all_white_family_count_and_rank_marginals=True,
            separate_ideal_and_finite_family_marginals=level >= 4,
            joint_closed_face_moment_matching=level >= 4,
            common_edge_component_state=level >= 4,
            full_discrete_realizability_claimed=False)
    return model
