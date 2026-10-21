# Independent review of the literal golden R5 chamber period

20 October 2026. Outside-author paper review: **PASS**, conditional on the
explicit inherited full-arrangement inputs. No scientific code, import,
generation, replay, Lean or TeX was executed for this review. This does not
record a successful fixed computational check.

Bound source:

* `results/golden_dyadic_r5_chamber_period_candidate_20261020.md`, SHA256
  `e2aeff93b42f87673e55399804df28a9eee3713ac2e7cf064c6f885e5bfce671`.

I read the complete source and independently checked the geometry and
algebra below. I also reread its actual-wall inputs
`golden_dyadic_unit_vertex_roofs_20261020.md` (SHA256
`e0c6b78e7c9c9ece0b77ef75cdfdcdc54d25c5177b1b474a3bcbead9515653e6`)
and `golden_dyadic_a3_b2_b3_actual_edge_20261020.md` (SHA256
`400ba9e7fb1369816a9cdfa10d51a22a75a5babec0a3b5aa7f128e2fbb9327dc`).
The latter already has my separate full edge review. The frozen M6/H3
packet and complete primitive norm pool {1,2,4} remain imported inputs.

## Full faces, interval and simultaneous mirrors

Both H4 frames have the stated Gram and unit determinant phi^-4.
The displayed inverse gives the two residue contractions, P1^2=-2phi,
g^2=1 and P1.g=0 exactly. Their explicit formulas are integral.
The finite Coxeter classification seals each actual H4 face: no new
actual wall meeting it can have a nonorthogonal finite extension.
The positive acute projection of an actual point to the positive
rank-eight system W is feasible for every actual wall. Projecting V0
onto its binary orthogonal space gives precisely S in equation (2).

The six rows exhaust the finite weighted completions of B3+A1 under
the imported norm/pairing rules. In particular a second 4-bond prevents
a bridge, and the norm1 A1 cannot acquire a same-length simple bond
with an integral lattice pairing. For any earlier actual wall, its
completion with W is positive definite and thus belongs to these rows.

I derived the collision bound b >= (phi^4/2)|a| and both embedding
bounds (4). The four negative A=2a values in the widest strip are
complete; the narrower q1 and q2 strips leave only -phi^-1. The
resulting q1/q2 values of b^2 are not in O. For q4, the actual
reflection divisibilities imply a in O and b in 2O, leaving a=-phi^-1;
both remaining b^2 values have valuation one at inert 2 and cannot be
squares of elements of 2O. No chosen height cutoff is used.

At a=0 the two surviving positive primitive candidates are g and Q.
Their simultaneous intersection is retained. The identity
g=Q+beta6 makes g redundant in the full chamber, since Q and beta6
are independent oriented integral mirror inequalities. Q has product
-1 with the norm1 beta6 and is the unique new simple wall. Thus the
new full H4+B3+B2 germ follows, rather than being assumed from a
binary endpoint or an unweighted diagram.

## Integral period and actual chamber matching

The old corner is actual by the same positive acute projection
principle and has the full stated rank-nine positive Gram. The ten
frame images preserve every block Gram and map the future P_old to
future P1, with beta5 mapped to Q and beta6 fixed. This includes the
inward B2 signs, not just the H4 subspace.

I independently derived Dp0 from the dyadic glue identity, then DV0
from p0=phi^3 V0+phi^2(e1+e2+e3+e4), and all remaining columns from
F=phi^-1 N-phi^-2 DV0. These reproduce every column of the displayed
matrix and its trace 10+11phi. The explicit original-coordinate basis
has determinant, up to sign, 2phi^2 relative to the time plus ordinary
D4 coordinate module; D4 has index 2 there. Thus it is indeed a basis
of the literal R5 lattice, not a sublattice of unspecified index.
Integral entries and exact Gram preservation give det D=+1 or -1;
the integral adjugate supplies the integral inverse.

D therefore preserves the full integral root arrangement. Matching
both complete actual inward vertex germs and the time orientation
implies D(C)=C. The trace exceeds 10, so D has infinite order. This
contradicts compactness of C, whose discrete stabilizer acts faithfully
on its finite spanning vertex set.

The infinite-volume conclusion also has a direct index formulation:
the reflection group W acts freely on arrangement chambers. If two
powers of D defined the same coset modulo W, their nonzero difference
power would belong to W and stabilize C, hence would be the identity.
Thus W has infinite index in the integral orthogonal group. The
arithmetic lattice input named in the source then implies infinite
covolume. This formulation avoids any separate inference about a
particular fundamental domain from compactness of a quotient.

No flaw or source correction was found. The conclusion concerns only
the literal R5 full reflection arrangement. It is not a similarity- or
commensurability-class exclusion, an all-dimensional family theorem,
a new example, or a uniform dimension bound. Fixed exact execution
remains a separate root-controlled verification step.
