"""One Lorentz Gram domain for every triangle or quadrilateral, at all ranks."""
from fractions import Fraction as F
from itertools import combinations

from . import data
from .moments import dimension_domain
from .polynomial import as_poly, determinant, var


def corner_flags(d, suffix="", angle=True):
    e = var("e")
    eta = d.variable("eta"+suffix, boolean=True)
    right = d.variable("right"+suffix, boolean=True)
    chi = d.variable("chi"+suffix, boolean=True)
    d.ge("vertex_kind"+suffix, 1-eta-right)
    d.eq("ideal_right_disjoint"+suffix, eta*right)
    d.eq("chi_ideal"+suffix, chi*(1-eta))
    d.eq("chi_dimension58"+suffix, chi*(1-e))
    if angle:
        c = d.variable("c"+suffix)
        d.eq("ideal_angle"+suffix, eta*(1-c))
        d.eq("right_angle"+suffix, right*c)
        d.ge("connected_finite_angle"+suffix, c-var("t")*(1-eta-right))
        return c, eta, right, chi
    return eta, right, chi


def face_domain(sides, threshold):
    if sides not in (3, 4):
        raise ValueError("Use the uniform small-face schema; residual corners have their own measure")
    d = dimension_domain(f"face{sides}", threshold)
    c, h, eta, right, chi = [], [], [], [], []
    for i in range(sides):
        ci, ei, ri, bi = corner_flags(d, str(i))
        c.append(ci); eta.append(ei); right.append(ri); chi.append(bi)
        hi = d.variable(f"h{i}", bounded=False)
        zi = d.variable(f"inverse_h{i}", bounded=False)
        fi = d.variable(f"finite_gap_inverse{i}", bounded=False)
        h.append(hi)
        d.eq(f"strict_common_interior{i}", hi*zi-1)
        d.eq(f"finite_corner_gap{i}", (1-ei)*((1-ci)*fi-1))
        for name in (f"h{i}", f"inverse_h{i}", f"finite_gap_inverse{i}"):
            d.basis("unbounded_"+name, (name,), 2)
    H = [[as_poly(int(i == j)) for j in range(sides)] for i in range(sides)]
    for i in range(sides):
        j = (i+1) % sides
        H[i][j] = H[j][i] = -c[i]
    if sides == 4:
        for name, (i, j) in zip(("opposite0", "opposite1"), ((0, 2), (1, 3)), strict=True):
            value = d.variable(name, bounded=False)
            H[i][j] = H[j][i] = -value
        d.eq("shared_four_wall_rank", determinant(H))
    spatial = [[H[i][j]+h[i]*h[j] for j in range(sides)] for i in range(sides)]
    d.matrix("common_spatial_Gram", spatial)
    for indices in combinations(range(sides), 3):
        selected = [[spatial[i][j] for j in indices] for i in indices]
        d.eq("spatial_rank_two_"+"".join(map(str, indices)), determinant(selected))
        selected_H = [[H[i][j] for j in indices] for i in indices]
        d.ge("Lorentz_three_minor_"+"".join(map(str, indices)), -determinant(selected_H))
    d.ge("right_angle_count", (1 if sides == 3 else 3)-sum(right, as_poly(0)))
    d.basis("shared_variables", tuple(sorted(d.variables)), 1)
    gram_variables = [f"c{i}" for i in range(sides)]+[f"h{i}" for i in range(sides)]
    if sides == 4: gram_variables += ["opposite0", "opposite1"]
    d.basis("shared_Gram_degree_four", gram_variables, 2)
    return d


def corner_domain(threshold, atom=None, angle=True):
    d = dimension_domain("corner", threshold)
    values = corner_flags(d, angle=angle)
    if angle:
        c, eta, right, _ = values
    else:
        eta, right, _ = values
        c = None
    if atom is not None and atom:
        if atom == 2:
            d.eq("atom_ideal", eta-1); d.eq("atom_not_right", right)
        else:
            d.eq("atom_finite", eta)
            cross = atom == 3 or atom >= 19
            d.eq("atom_component_relation", right-int(cross))
            if angle and not cross:
                minimum = F(1, data.SAME_DISTANCE[atom-4]+1) if 4 <= atom <= 18 else F(1, 17)
                d.ge("atom_inverse_Gram_angle", c-minimum)
    d.basis("corner", tuple(sorted(d.variables)), 1)
    return d
