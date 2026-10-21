# Hostile audit: equivariant minimalization and trace descent

**Verdict (2026-08-10): PASS after proof-completeness and scope repairs.**
For a full finite-covolume trace representation on a smooth projective
surface, a finite numerical operator orbit of one genuine exceptional
lattice is equivalent to orthogonality of that lattice to the trace block,
and hence to integral isometric descent of the trace block through the chosen
contraction.  Stability of the exceptional lattice additionally makes the
original integral operators descend.  No complement-rank bound, finite
enumeration of contractions, or raw-polytope-to-carrier theorem is obtained.

Audited theorem:
[equivariant_minimalization_trace_descent_theorem.md](equivariant_minimalization_trace_descent_theorem.md).

## Pinned artifacts

| artifact | SHA-256 |
|---|---|
| `results/equivariant_minimalization_trace_descent_theorem.md` | `ff8a89ae37a24ad2a10515f5acd931b15a280ed8f7fafb77e1db0c3d6b61ae5d` |
| `scripts/certify_equivariant_minimalization_trace_descent.py` | `4ee02fb4683e5d98f2cf8d458f37860a0cd4ec843b1a6b2b784732a3db09872e` |
| `data/enumerations/equivariant_minimalization_trace_descent.json` | `490c3f73c9ea311a3f3da5e45c2b4b27250ed588b446e2744a443256876ded37` |
| `tests/test_equivariant_minimalization_trace_descent.py` | `77a52428a9a307cc8f1c938be7ddc834e8402acbd5ec68b1f5292d9f5659d452` |

## 1. Exceptional lattice and integral descent

A birational morphism between smooth surfaces factors into blowups at
points.  In the total-transform basis of the successive exceptional classes,

```text
Num(Y) = pi^*Num(X) orthogonal_sum <-1>^r
```

integrally.  Thus `E_pi` is primitive, unimodular, and negative definite,
and `E_pi^perp=pi^*Num(X)` as integral lattices.  The prime exceptional
curves themselves need not be an orthogonal basis when centres are
infinitely near.  The revised proof explicitly uses this integral equality;
the rational orthogonal decomposition alone would not justify an integral
map `T -> Num(X)`.  The surface factorization input agrees with the
[Stacks Project, Tag 0C5Q](https://stacks.math.columbia.edu/tag/0C5Q).

Once `T` is orthogonal to `E_pi`, every integral trace vector therefore has
a unique integral preimage on `X`.  If each `u_i` preserves `E_pi`, it also
preserves its integral orthogonal summand and descends to an integral pairing
isometry.  Splitting

```text
K_Y = pi^*K_X + K_(Y/X)
```

between those stable summands proves that the descended operator fixes
`K_X`.  This part is valid for arbitrary sequences of point blowups, not only
blowups at mutually distinct points.

## 2. The finite-covolume representation step

The integral orthogonal group of a negative-definite lattice is finite: its
elements permute the finite set of lattice vectors below any fixed norm
bound, and a sufficiently large such set spans.  Hence the kernel
`Gamma_0` of the action on `E_pi` has finite index.

There was a hidden image step in the draft, now repaired.  Restriction to the
trace block maps `Gamma` onto the Coxeter group `W(P)`, since the displayed
generators restrict to all its simple reflections.  Therefore the image of
`Gamma_0` is finite-index in `W(P)` and remains a finite-covolume hyperbolic
lattice.  At the distinguished real place it fixes no nonzero vector: a
timelike, null, or spacelike fixed vector would force respectively a fixed
interior point, fixed ideal point, or invariant proper totally geodesic
subspace, contradicting its full-sphere limit set.

An element of `Res_(k/Q)(V)` whose distinguished component vanishes is zero:
write it in a `k`-basis and use injectivity of the distinguished field
embedding.  Consequently `T_Q^(Gamma_0)=0`.  The nondegenerate invariant
trace form identifies invariant covectors with invariant vectors, so every
equivariant map from `T_Q` to a trivial finite-dimensional target vanishes.
Applied to the exceptional projection, this proves `T perp E_pi` without an
irreducibility assumption on `T`.

Finite covolume and the full trace representation are essential here.  For
an infinite-covolume subgroup preserving a proper hyperbolic subspace, or
for a partial root subsystem with a fixed trace summand, the invariant-vector
step can fail and a nonzero exceptional projection is possible.

## 3. Finite-orbit equivalence

The strengthened Theorem 2.2 is sound.

- If the orbit `mathcal E={u(E_pi)}` is finite, take the external orthogonal
  direct sum of its negative-definite integral lattices.  `Gamma` permutes
  the finitely many summands, and after killing that permutation action each
  stabilizer acts through a finite group `O(E)`.  The total image is finite.
  The simultaneous projection `x -> (p_E(x))_E` is equivariant, so the same
  finite-index no-covector argument kills every component.

- If `T perp E_pi`, every translate of `E_pi` lies in
  `C=T_Q^perp intersect Num(Y)`.  By Hodge index, `Num(Y)_Q` has one positive
  direction and the nondegenerate trace block already contains it.  Thus
  `C` is negative definite, `O(C)` is finite, and the orbit is finite.

- Orthogonality is equivalent to containment in `pi^*Num(X)` by the integral
  exceptional splitting, giving the third condition.

The other orbit members are merely numerical sublattices.  They need not be
effective, contractible, or exceptional lattices of other MMP morphisms.
The theorem requires one initial genuine contraction and does not promote
its translates to geometry.

When the orbit is finite but not a singleton, the original integral `u_i`
need not descend.  What descends is the integral trace lattice; the report
then constructs rational orthogonal-sum extensions on `Num(X)_Q`.  This is
exactly enough for the earlier rank theorem, which requires rational pairing
isometries fixing `K_X`, but it is not an integral or geometric action of the
original group on the minimal model.

## 4. Exact matrix gate and its boundary

For a certified geometric contraction frame, with full trace-inclusion
matrix `J`, exact intersection matrix `G_Y`, and a full total-transform
exceptional-frame matrix `C_pi`,

```text
D_pi = J^t G_Y C_pi
```

is precisely the matrix of trace-exceptional pairings.  Therefore

```text
D_pi = 0
  <=> T perp E_pi
  <=> isometric descent through pi
  <=> the numerical Gamma-orbit of E_pi is finite.
```

This is a genuinely useful finite exact gate **per supplied contraction**.
It does not show that the set of geometric contraction frames is finite or
effectively enumerable, and it cannot be run on an arbitrary negative frame
or an incomplete list of prime exceptional curves.  Calling it a finite
MMP search theorem would be an overclaim; the report now excludes that
reading.

The calibration rows have also been checked at their actual strength.  The
six minimal K3 or canonical-K3-cover rows use `pi=id`, so their rank-zero
matrices are vacuous zero controls, not evidence that a nontrivial
contraction descends.  The `Gamma16` blowups are the nontrivial positive
controls.  For the five rational full-NS rows the replay chooses the standard
blowup basis and forms `J`, `G_Y`, `C_pi`, and `D_pi` exactly; the resulting
defect ranks are `6,10,11,14,16`.  The Type-VI rank-fourteen plane frames and
the Type-MII rank-sixteen `P1xP1` frame are independently pinned by the
Coble blowdown archive.  Thus those ranks are no longer inferred from Picard
rank alone.

The printed rational matrices are normal-form representatives, not a
reconstruction of `C_pi` in each upstream source's original coordinates.
This loses no information used here: rank and vanishing are invariant under
integral changes of the surface, trace, and exceptional-frame bases.  All
five rational rows remain defect-only controls because their archived data
do not supply the theorem's canonical-preserving reflection extensions.
The Type-VI and Type-MII rows are in characteristic three and therefore do
not enter the complex Picard ceiling `20`; under the analogous positive-
characteristic hypotheses the unconditional ceiling is `22`.

## 5. Actual automorphisms and the rational diagnostic

For `kappa(Y)>=0`, the smooth minimal model is unique in its birational
class.  If `g` is an automorphism of `Y`, the birational map
`pi g pi^(-1)` of the minimal model is biregular; hence `pi g=g_X pi` and
`g` preserves the contracted curves and `E_pi`.  The stable version of the
theorem applies, and the already-audited minimal-surface rank theorem gives

```text
d(n+1) <= 20
```

and Kodaira dimension zero.  This does not apply to arbitrary Cremona
actions on rational or ruled surfaces.

The rational diagnostic is exact.  A minimal rational surface has numerical
rank one (`P2`) or two (a minimal Hirzebruch surface), whereas the full trace
block has rank `d(n+1)>=3` for `n>=2`.  Thus no exceptional lattice of a
contraction to such a minimal surface can have finite numerical orbit, or
equivalently be orthogonal to the full trace block, under all the theorem's
operator hypotheses.  The conclusion is numerical: it does not say that
there are infinitely many effective MMP contractions.

## 6. Gamma16 blowup family and correspondence formula

For the rank-seventeen lattice

```text
S = U + E8(-1) + E7(-1),
```

the archived nineteen-wall Gram has rank 17, inertia `(1,16,2)`, and an
integral determinant-two wall basis.  The upstream maximality archive is
used as a pinned dependency; this certifier checks its hash and conclusion
but does not independently reprove reflection maximality.

Blowing up `r` distinct points gives

```text
NS(Y_r) = pi_r^*S orthogonal_sum <-1>^r,
K_(Y_r) = E_1 + ... + E_r.
```

For a root `a`, the integral cycle `Z=Delta_X+a times a` induces
`x -> x+(x.a)a`, the square-minus-two reflection on `NS(X)`.  Under the
covariant correspondence convention, the graph-composition term in (3.5)
acts as `pi_r^* s_a pi_(r,*)`: it is `s_a` on the pullback block and zero on
the exceptional block.  A product cycle `E_j times E_j` sends `x` to
`(x.E_j)E_j`; because `E_j^2=-1`, its negative is the identity on `E_j`.
Thus the displayed sum is exactly `s_a direct_sum 1_r`.  The revised
certifier checks this correspondence matrix for all nineteen walls and each
sample blowup count, rather than checking only the desired block matrix.

The construction is uniform for every `r>=0`, while samples
`0,1,2,5,11` are regression checks.  The complement rank **inside
`Num(Y_r)`** is exactly `r`, and `K_(Y_r)^2=-r`; the raw Gamma16 chamber and
its maximality remain unchanged.  This proves that raw reflection-maximal
Coxeter data do not bound complement rank across the stated redundant
correspondence-carrier category.

It does not make these correspondences automorphisms, assert that the full
nef cone of `Y_r` has nineteen walls, or produce new reflection overgroups.
Indeed, the pullback nef section is the old chamber while every exceptional
curve contributes extra geometry outside that section.  Those exclusions
are necessary to the sharpness claim.

## Replay

```text
python3 -m scripts.certify_equivariant_minimalization_trace_descent --check-archive
python3 -m pytest -q tests/test_equivariant_minimalization_trace_descent.py
```

The focused suite has nine passing tests.  The machine certificate verifies
the frozen Gamma16 inputs, all nineteen
integral reflection matrices, their isometry and involution equations,
canonical-class fixing, the explicit correspondence formula, trace and
exceptional pairings, signatures, complement ranks, the standard-basis
defect matrices, and the pinned Type-VI/MII blowdown-frame predicates.  It
certifies the sharpness family and the stated calibrations; the general
descent and finite-orbit equivalence are the symbolic arguments audited
above.

## Remaining global gap

Nothing here constructs an intrinsic full trace carrier from an arbitrary
maximal arithmetic Coxeter chamber.  Even after a carrier is supplied, no
theorem yet forces the existence of a geometric contraction frame with
`D_pi=0`, nor finiteness of the search over such frames.  This result replaces
an impossible complement-rank target with an exact descent invariant; it
does not yet close the universal polytope-to-minimal-surface bridge.
