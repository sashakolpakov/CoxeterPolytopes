# A bounded next target in the nonarithmetic H9 chamber

The next concrete experiment is a **two-cap extension of the existing 33-wall
octic bridge chamber**, supported on

`T = [0,1,2,5,6,10,30]`.

The support is positive definite of type `H3+A3+A1`, and the surviving strict
spacelike escape has a strictly spacelike projection to its Lorentzian
three-dimensional complement. Two actual finite vertices already lie on
this two-face. These are useful prerequisites for an escape-cutting cap
cluster; they do not prove that such a cluster exists. The proposed search
below has **not** been run.

The [checker](../scripts/fresh_nonarithmetic_escape_target_20261016.py) and
[exact certificate](fresh_nonarithmetic_escape_target_20261016.json) perform
only the bounded diagnostics described here. Generation and a separate saved
replay pass. They reuse the frozen independent octic arithmetic and sign
engine; this is not another independent arithmetic review of the old chamber.
No frozen predecessor was changed.

## Why this route differs from the excluded lattices

The [33-wall bridge](fresh_h9_joint_corner_20261015_bridge_review.md) retains
the original changed H8 lift, all its later geometric walls, and a normal
in an octic field. It has 92 finite vertices, 253 finite edges and 322 open
finite-corner directions. Its exact escape remains strict. A complete
exterior-ray census of this 33-wall chamber was not claimed in its certificate.

There is also a literal obstruction to an arithmetic completion. Retained
walls11 and12 have Cartan two-cycle

`4 (n11,n12)^2 / (n11^2 n12^2) = 4 sqrt(5)`.

The checker reconstructs this equality. Its nonidentity conjugate is negative,
contradicting the necessary conjugate-positive-semidefinite condition for a
quasi-arithmetic finite-covolume completion. This is the same inherited
[two-wall obstruction](fresh_h9_h8_disjoint_cap_20261014.md). Thus the target
is a nonarithmetic geometric completion, with no integral lattice or
primitive-reflector restriction imposed on the added normals. The present
infinite-volume reflection chamber is not being called a nonarithmetic lattice.

The retained H20 chamber gives a less useful immediate target: its compact
S61 polygon excludes every additional perpendicular acute wall from cutting
the saved escape, and the exact
[847-profile nonright bank](fresh_h20_nonright_s61_20261014.md) still admits a
common strict spacelike escape. That bank does not exclude all nonright
walls. The H9 target below changes the support and uses its actual finite
corners, without repeating the H20 bank or a new integral-host search.

## An exact obstruction for the first tempting support

Let `S=[0,1,2,5,6,10,29]`. For the literal saved escape E, project orthogonally
to `S` by

`P_S = E - sum_i (G_S^-1 [(E,n_j)]_j)_i n_i`.

The seven-root Gram is positive definite. The checker reconstructs the actual
future finite rays with indices20,21,90, normalized by `(V_i,H)=-1`, and
proves the exact identity

`P_S = lambda20 V20 + lambda21 V21 + lambda90 V90`,

with **all three coefficients strictly positive**. It stores the vectors and
coefficients. Their full nine-root incident Grams are positive definite, and
each vector has exactly its recorded zero set and strict negative pairings
with all other 33 walls.

This excludes every coherent acute additional wall perpendicular to S,
without an angle cutoff or a field assumption. To see the strictness, let
v be any existing finite corner, h an interior point of the extended chamber,
and n_i its old incident roots. Write

`h = a v + sum_i b_i n_i`, with `a>0` and
`b=G_v^-1[(h,n_i)] <= 0`.

The inverse of the positive-definite Z-matrix G_v is entrywise nonnegative.
For a proposed outward normal r acute to every old root,
`(r,h) >= a(r,v)`; hence `(r,v)<0`. This applies even if the old chosen
center H does not remain inside the extension. If r is perpendicular to S,
then `(r,E)=(r,P_S)<0` by the displayed positive combination. Such a wall
cannot cut E. This obstruction is stronger than merely observing a timelike
projection.

## The neighboring support and a failed polar shortcut

Replace wall29 by wall30 to obtain T. Its complete finite-vertex incidences
in the inherited census are

- vertex21: `T + {29,31}`;
- vertex91: `T + {29,32}`.

They are joined along the rank-eight support `T+{29}`. The checker verifies
both full corners and the strict spacelike norm of `P_T(E)`. The projections
of the three side normals to `T^perp` have exact squares

| side | projected square |
|---|---:|
| 29 | 3/2 |
| 31 | 1/2 |
| 32 | 3/2 |

Consequently a norm-two new wall perpendicular to T that forms a **finite**
corner with wall31 must meet31 at a right angle. At label3 the two-root
residual determinant is `2*(1/2)-1=0`; every larger finite label makes it
negative. This is an all-label necessary restriction on that endpoint.

The unique polar line perpendicular to `T+{31,32}` was tested exactly. Choose
the orientation with `(r,H)=-1`. It is spacelike, but `(r,n29)>0` and
`(r,E)<0`. The opposite orientation is obtuse to wall4. Thus neither polar
orientation can supply the desired cap. No claim about other normals follows.

## The finite two-cap experiment

Retain all33 old walls. Seek norm-two normals u and v perpendicular to T.
Use the following finite recipe family:

1. Choose `j in {3,7}` and `m in {2,3,4,5}`. Impose
   `(u,n31)=0`, `(u,n_j)=-2 cos(pi/m)`, and `u^2=2`.
2. Choose `n,ell in {2,3,4,5}`. Impose
   `(v,n32)=-2 cos(pi/n)`, `(u,v)=-2 cos(pi/ell)`, and `v^2=2`.

In the three-dimensional complement, each step intersects two linear
conditions with one quadratic. There are128 recipes and at most512 isolated
ordered branch pairs. Any rank-deficient or continuous intersection must be
recorded separately; it is not covered by that isolated-branch count. The
two selected anchors3 and7 are outside T and belong to the active support of
the literal saved escape.

A successful proposal must satisfy all67 old/new and mutual Coxeter gates,
have a common strict future interior point, make both walls actual facets,
and have at least one positive pairing with E. Require the three supports

`T+{31,u}`, `T+{u,v}`, `T+{v,32}`

to be positive-definite full rank-nine Grams with feasible future vertices.
This asks for a connected finite-corner chain across the two existing open
directions, not just another wall count. A simple first implementation may
require every unspecified new pairing to be right or strictly dotted; that
is an explicit additional restriction on the experiment, not a complete
all-label test of the entire family.

The search is intentionally a finite proposal family. Failure would not
exclude arbitrary endpoint caps, larger clusters, other labels, changed
retained walls, or a finite-volume completion. Success would still require
a new escape/global-volume analysis. No proposed cap has been accepted here.

## Read-only checks requested for the current checkpoint

No logical scope inconsistency was found in these final source/proof reads:

- The [F25 independent review](fresh_sqrt2sqrt5_lanner_catalogue_20261016_review.md)
  uses **exact** ground field Q(sqrt2,sqrt5). Its higher-order Lanner exclusion
  invokes the same-ground-field condition, not containment alone. The
  triangle determinant argument applies to cyclic products even when
  individual normalized entries lie outside F. The negative-plane
  quadrilateral obstruction is a principal-submatrix argument, so extra
  white nodes cannot repair it.
- The [single-field composition](fresh_compact21_single_field_20261016.md)
  removes only F25 from the older two-field reduction. It retains Q(sqrt5)
  and explicitly sets `dimension21_excluded=false`. It uses the field list
  valid from dimension14, not the shorter restriction beginning at22.
- The [Witt setup review](fresh_h9_witt_twist_20261016_review.md) correctly
  limits its conclusion to local index, norm classes and the zero sector.
  The primitive norm argument is `beta G^-1` integral, hence `q | 2 beta`;
  the local comparison excludes similarity to E8-plus-binary forms. Neither
  conclusion asserts reflectivity. The subsequent
  [T20 symmetry review](fresh_h9_witt_twist_symmetry_20261016_review.md)
  therefore does not contradict this setup.

The first three source/report/certificate triples are hash-bound by this
diagnostic. Their arithmetic replays are inherited, not rerun here. The
published F25 inputs remain those identified in
[Nikulin sections2 and3.1](https://arxiv.org/pdf/0708.3991). These checks add no
new dimension, novelty or current-best claim.

```sh
python3 scripts/fresh_nonarithmetic_escape_target_20261016.py
python3 scripts/fresh_nonarithmetic_escape_target_20261016.py --verify
```
