# Independent hostile audit: corrected genus-one packet and the exact `P_n` boundary

**Verdict (2026-08-11): PASS after two scope repairs.**  The local packet is
indeed rank-blind, the five positive full-nef carriers replay, and the
marking-independent `P17`--`P19` obstruction is supported by the exhaustive
surface-MMP archives.  The hostile pass repaired two statements without
changing the theorem:

1. the putative plane cubic in the no-go proof need only be integral and
   smooth at the first centre, not globally smooth; and
2. the obstructing quartic is
   \[
   A_n=4H-E_1-\cdots-E_{\min(n,18)},
   \]
   so `A17` follows seventeen centres while `A18` and `A19` follow eighteen.

The second repair is literal, not cosmetic: the frozen `P18` and `P19` root
packets do not contain the `P17` quartic class.  The Abel argument survives
because all three actual quartics share exactly the first eleven centres
with the cubic.

## 1. The corrected local packet has no rank information

In `I_(1,n)`, for every `n>=12`, put

\[
 K=(-3,1^n),\qquad
 C=(3,-1^{11},0^{n-11}),\qquad
 Z=C+K=(0^{12},1^{n-11}).
\]

Direct calculation gives

\[
 C^2=-2,qquad K.C=2,qquad p_a(C)=1,qquad
 C.Z=0,qquad Z^2=11-n.
\]

Because `C` is primitive of square `-2`, its reflection is integral in every
rank.  The independent executable repeats these identities through `n=64`;
the invariant tuple remains exactly `(-2,2,1,integral)`, while `Z^2` absorbs
all rank dependence.  This proves only that the local NS/canonical data
cannot yield a bound.  It is not a construction of a full-nef carrier in
every rank.

## 2. The positive range is a complete-cone theorem

The source archive was replayed and SHA-pinned, including its exact CoxIter
relation matches and its absolute-irreducibility replay for the cubic,
quartic, and sextic.  The positive cases are not short-curve searches.

- For `P12` and `P13`, proximity and the two Bezout inequalities
  `A3<=d`, `A11<=3d` give a symbolic nonnegative decomposition of every
  possible irreducible plane-curve class in `L,C` and the exceptional chain.
- For `P14` and `P15`, the extra quartic inequality cuts out exact rational
  cones with 94 and 145 primitive rays.  Every ray has an exact nonnegative
  `L,C,Q` plus exceptional decomposition.
- For `P16`, the quartic and sextic inequalities cut out 325 primitive rays,
  each with an exact nonnegative `L,C,Q,S` plus exceptional decomposition.

Consequently the displayed negative curves generate the complete Mori
cones, and duality gives literal full-nef chambers `P12` through `P16`.
The unique arithmetic-genus-one facet is the smooth strict cubic `C` in
these explicit positive constructions.

## 3. Why the MMP enumeration removes the marking hypothesis

Suppose a smooth projective rational surface had complete Mori cone equal
to one of the `P17`, `P18`, or `P19` facet cones.  The abstract facet set
contains a unimodular basis, checked independently in all three dimensions;
there is therefore no unexamined finite-index lattice embedding.

On a nonminimal rational surface, a `K`-negative extremal ray is generated
by a smooth rational `(-1)`-curve `e`.  Completeness of the facet cone puts a
current irreducible facet image on that ray.  It cannot be a proper multiple:
a distinct irreducible curve numerically equal to `q e` would meet `e` in
`-q`, while `e` itself has primitive square `-1`.  Thus every geometric MMP
is among the exact numerical branches

\[
                     D\longmapsto D+(D.e)e.
\]

Branches using a square-`-1` image of positive genus merely overcount and
are harmless.  Every surviving branch reaches `P2` or `F2`.  Reconstructing
the canonical class independently from both endpoints gives the same
sets, of sizes

\[
                         8,qquad24,qquad120.
\]

The hostile replay checks every stored weighted-diagram permutation against
the full facet Gram, verifies that the first `n+1` facets form a unimodular
basis, and recomputes the canonical pairing-profile orbit.  Those orbits
are exactly the archived MMP canonical profiles, again of sizes 8, 24, and
120.  After a diagram re-marking, `K=(-3,1^n)`; adjunction then makes the
standard `B_n` chain smooth rational, so successive Castelnuovo contraction
recovers the required free plane cluster.

## 4. The corrected Abel obstruction

In the standard marking the required facets are

\[
 L=H-E_1-E_2-E_3,qquad
 C=3H-E_1-\cdots-E_{11},
\]

\[
 Q=4H-2E_1-E_2-\cdots-E_{15},qquad
 A_n=4H-E_1-\cdots-E_{\min(n,18)}.
\]

The exact intersections are

\[
 L.C=C.Q=Q.A_n=0,qquad C.A_n=1.
\]

The image of `C` is an integral plane cubic.  Its coefficient one at the
first centre makes it smooth there.  From `L.C=0`, Bezout gives

\[
                \mathcal O_C(1)=\mathcal O_C(3p).
\]

The common multiplicities of `C` and `Q` already contribute
`2+10=12`, so the selected chain leaves `C` at the twelfth centre.  Each
actual `A_n` follows that chain and is smooth at `p`; hence

\[
                       A_n|_C=11p+r,qquad r\ne p,
\]

where `r` is an effective degree-one Cartier point.  But

\[
 \mathcal O_C(A_n)=\mathcal O_C(4)=\mathcal O_C(12p)
\]

forces `O_C(r)=O_C(p)`.  The degree-one Abel map is injective on the smooth
locus of an integral plane cubic, including the nodal and cuspidal cases;
therefore `r=p`, a contradiction.  This excludes every rational-surface
marking of `P17`, `P18`, and `P19` in characteristic zero.

## 5. What property 4 must remember

The stopped family identifies a concrete invariant beyond the Gram lattice:

> the oriented facet-contraction order must compile to an actual proximity
> cluster, and all required curve restrictions must be compatible in the
> Picard groups of positive-genus facets.

This is a genuine property-4/effectivity condition.  It explains both sides:
the Hilbert branches for `P12`--`P16` are nonempty, while the `P17`--`P19`
branches fail in `Pic^1(C)`.

The conclusion is deliberately scoped.  It is an iff theorem only for the
exact family `P12`--`P19`; it is not a universal `n<=16` theorem for every
maximal chamber with a genus-one facet, and it says nothing against the
classical abstract hyperbolic polytopes.

## 6. Replay

The independent executable is
[`audit_rational_pfirst_corrected_genus_one_pn_boundary_independent_hostile.py`](../scripts/audit_rational_pfirst_corrected_genus_one_pn_boundary_independent_hostile.py),
with separate tests in
[`test_rational_pfirst_corrected_genus_one_pn_boundary_independent_hostile.py`](../tests/test_rational_pfirst_corrected_genus_one_pn_boundary_independent_hostile.py).

```bash
/Users/sasha/.venvs/global/bin/python \
  scripts/audit_rational_pfirst_corrected_genus_one_pn_boundary_independent_hostile.py

/Users/sasha/.venvs/global/bin/python -m pytest -q \
  tests/test_rational_pfirst_corrected_genus_one_pn_boundary_independent_hostile.py \
  tests/test_rational_p17_mmp.py \
  tests/test_rational_p18_p19_mmp.py
```
