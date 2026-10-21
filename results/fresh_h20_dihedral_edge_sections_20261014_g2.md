# G2 sections: finite continuous reduction and three exact exclusions

**Exact result:** three specified continuous one-contact `A2 -> G2` source
edge families fail for every transverse parameter. Their complete necessary
face-closure census contains 131 distinct section normals, each rejected at
an actual finite neighboring vertex: 105 have an obtuse ridge and 26 have
an angle strictly between `pi/3` and `pi/2`.

The broader scout traversed all 5,048 anchored literal source families and
found no positive proposal. Its 481,016 ridge rejections are **numerical**,
not an exact theorem for every family. All 66 numerically unresolved or
near-boundary recipes were checked exactly and are double roots at old
endpoints, outside the transverse interior family. The label-6 construction
as a whole remains open until the other numerical rejections have been
certified or replaced by a uniform argument.

## Why a continuous family reduces to finitely many candidates

For old norm-two roots `a,b` with `(a,b)=-1`, normalize
`(e,a)=1`, `(e,b)=0`. Label 6 gives `e²=3/4`. Decompose

\[
 e=w_a+z,\qquad w_a^2=\frac23,\qquad z^2=\frac1{12},
\]

where `z` is in the Lorentzian two-plane orthogonal to the old 20-wall edge.

The actual `G2` face has dimension 18. It need not be compact: `G2` has
affine extensions. Every one of its bounding facets at a finite point,
however, is orthogonal to both `G2` normals, since `G2` has no proper
connected finite spherical extension. The initial 18 other walls cannot
bound a finite-volume dimension-18 face. Therefore at least one additional
old wall `l` must bound it. With `g=(l,a)`, this requires

\[
 (l,b)=0,\qquad (z,l)=\frac{g}{12}.
\]

The latter linear condition and `z²=1/12` have at most two solutions in the
edge plane. The old wall's plane projection is nonzero: otherwise it would
contain the entire source edge, contrary to the strict source inequalities
off the known 20 incident walls. This is a finite necessary closure
reduction, not a grid of edge positions.

## Complete exact families

The specified recipes are anchored literal supports from the closed source
corner walk; their IDs refer to the saved scout catalogue.

| Family | Source corner / edge index | A2 anchor | Endpoint types | Distinct transverse normals |
|---|---|---|---|---:|
| 0 | 3 / 1 | 22 in `{22,40}` | ideal–finite | 39 |
| 34 | 3 / 57 | 11 in `{11,39}` | ideal–ideal | 76 |
| 660 | 15 / 1 | 3 in `{3,40}` | finite–finite | 16 |

The [exact checker](../scripts/fresh_h20_dihedral_edge_sections_20261014_g2_selected.py)
reconstructs **every** eligible old wall and every real root of its closure
quadratic. It tests the crossing signs exactly, includes both branches,
and identifies duplicate normals with exact coefficient keys. Floating-point
data only proposes which existing wall supplies a defect witness; no root,
sign, or completeness decision is trusted from that data.

For each candidate it rebuilds a future timelike transverse point, the
section normal, and the dual direction for the proposed released wall.
The first contact is positively scaled to avoid field division. It checks
all 210 source inequalities and the exact 20-wall active set. The retained
19 projected normals are positive definite; the new wall is independent
because its pairing with the released dual direction is nonzero. Thus the
20 active walls form an actual finite vertex. A positive two-by-two Gram
determinant and the exact zero set establish a genuine defective ridge.
Potentially redundant remote walls are not used as a substitute for this
actuality check.

The [saved exact certificate](fresh_h20_dihedral_edge_sections_20261014_g2_selected.json)
records all 131 normals, contact points, signs and defects. Both generation
and its saved-result replay passed.

## Broader numerical census and exact endpoint cleanup

The [numerical producer](../scripts/fresh_h20_dihedral_edge_sections_20261014_g2_scout.py)
examined all 5,048 anchored A2 families in the 6,492 saved literal source
edge supports. Its [summary](fresh_h20_dihedral_edge_sections_20261014_g2_summary.json)
and [compressed complete archive](fresh_h20_dihedral_edge_sections_20261014_g2_scout.json.gz)
record 481,082 numerically distinct proposals:

| Numerical status | Count |
|---|---:|
| Proposed actual-ridge defect | 481,016 |
| No defect certified numerically | 61 |
| Boundary or conditioning ambiguity | 5 |

The [separate exact endpoint checker](../scripts/fresh_h20_dihedral_edge_sections_20261014_g2_remaining.py)
reconstructs all 66 remaining quadratics and signs from the discrete recipes.
Every one has discriminant exactly zero and a normal orthogonal to one old
endpoint. These are not interior transverse sections. Its
[certificate](fresh_h20_dihedral_edge_sections_20261014_g2_remaining.json)
and saved replay pass. This resolves those 66 recipes only; it does not turn
the other 481,016 numerical failures into exact exclusions. Floating-point
deduplication and near-boundary filters also remain numerical in the broad
scout.

```sh
python scripts/fresh_h20_dihedral_edge_sections_20261014_g2_selected.py --verify
python scripts/fresh_h20_dihedral_edge_sections_20261014_g2_remaining.py --verify
```

No job remains running. No finite-volume Coxeter example was obtained.

## A constructive next reduction for label 4

The corresponding `B2` face admits finite `B3` boundary extensions, so the
all-orthogonal `G2` reduction cannot simply be reused. It still gives a
finite route. Put `e²=1`, `g_a=(l,a)`, `g_b=(l,b)`, `d=l²`, `x=(e,l)`.
An additional bounding wall must be one of:

1. Orthogonal: `g_b=0`, `x=g_a`.
2. Label 3 at the modified end: `g_b=0`,
   `5x²-8g_a x+4g_a²-d=0`, with `g_a-x<0` and positive projected norm.
3. Label 3 at the unchanged end: `x=g_a`,
   `d-g_a²=2g_b²`, with `g_b<0`.

Each contact value imposes one linear condition on the Lorentzian edge
plane, leaving at most two normals. In particular the irrational `B3`-join
values must be included. This is a proposed next finite census, not a search
or completeness certificate already executed here.
