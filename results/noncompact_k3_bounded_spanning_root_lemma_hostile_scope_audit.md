# Hostile scope audit: the cusped K3 bounded spanning-root lemma

**Verdict (2026-08-11): PASS as a theorem in the finite-sided exact
full-Picard/full-nef scope.** The new argument closes the precise global
gate left by the cusp-divisibility family. For every Picard rank
\(3\le\rho\le20\), a finite-volume exact K3 root chamber contains
\(\rho\) literal physical square-\((-2)\) facet vectors spanning
\(\operatorname {NS}_{\mathbf Q}\), with all off-diagonal intersections at
most \(13\).

This audit is a hostile scope check by the theorem owner, not an independent
second proof. It isolates each place where compactness might have been
silently reintroduced and records the counterfactual if that input is
removed. An independent audit should replay the same gates before the result
is used as a frozen Property 4 dependency.

## 1. The visual bound does not use compactness

The minimal-aperture calculation needs only:

- a finite set of facet supports, so a minimum aperture is attained;
- an interior observation point;
- projective closed-facet incidence with one chosen facet \(E\); and
- the elementary sine quotient bound.

If \(F\) and \(G\) each meet \(E\), split according to their supporting
hyperplanes. Intersecting or asymptotic supports have \(b<2\) or \(b=2\)
directly. Only for ultraparallel supports are the two visual cones disjoint.
In that case their positive angular gap is at most the aperture of \(E\),
even when the two facets meet \(E\) only at different ideal points.
Minimality of \(E\) then puts both sine quotients below two and gives

\[
                 b(F,G)<4\cdot2\cdot2-2=14.             \tag{1.1}
\]

No vertex link, compact facet, simple-polytope identity, or cusp type occurs
in this calculation. The use of the **closed star** is essential: affine
ridge adjacency would omit facets which meet only at infinity. It would be
incorrect to set the gap to zero for overlapping visual cones and retain
the ultraparallel sine identity; the direct \(b\le2\) branch is logically
separate.

## 2. Finite volume is used exactly once

At an ideal point \(u\), move \(u\) to infinity. Above a high horosphere,
all nonincident hemispheres disappear and the chamber is the product of a
Euclidean incident-wall section \(C\) with a vertical ray. Finite hyperbolic
volume forces \(C\) to be bounded. Boundedness of \(C\) is equivalent to
the Euclidean wall normals positively spanning the horosphere; lifting them
back to Lorentz space gives

\[
          \operatorname {span}\{r_F:u\in\overline F\}=u^\perp. \tag{2.1}
\]

This proof uses the complete incident set and works for a product of affine
Coxeter simplices or any other nonsimple bounded Euclidean link. It never
classifies that link.

If finite volume is dropped, \(C\) may be unbounded and (2.1) fails. If one
keeps only a displayed affine subdiagram, the chosen normals may also span a
proper subspace of \(u^\perp\). These are genuine failures, not cosmetic
proof gaps.

## 3. The two-point bridge is exact

The closure of the positive-dimensional facet \(E\) in the closed Klein
ball is a compact convex set and therefore has two distinct extreme points
\(x,y\). At a finite extreme point, active facet normals span \(x^\perp\):
otherwise a small timelike displacement in their common kernel produces a
segment through the alleged extreme point. At an ideal extreme point,
finite volume makes the full horospherical section bounded, so it is an
actual ideal vertex and (2.1) gives the same conclusion. This remains valid
when the cusp link is nonsimple.

Every facet through either point lies in the closed star of \(E\). Moreover,

\[
 (x^\perp+y^\perp)^\perp
   =\mathbf Rx\cap\mathbf Ry=0,                         \tag{3.1}
\]

because \(x,y\) are two distinct projective points. Hence their incident
packets together span the full rank-\(\rho\) Lorentz space. Since both
packets lie in the same closed star, (1.1) also bounds every cross-pairing.

The one-cusp control is sharp: one complete cusp packet has rank
\(\rho-1\) and spans only \(u^\perp\). Duplicating that packet does not add
the transverse direction. The second distinct point, not merely more roots
at the first point, is what repairs the cusp-divisibility escape.

## 4. The K3 conversion uses literal same-vector return

The normalized bound (1.1) becomes an integral Gram bound only because each
selected normal is the same physical square-\((-2)\) vector in the full
Picard lattice. For two distinct irreducible K3 curves,

\[
             b(E_i,E_j)=E_i.E_j\in\mathbf Z_{\ge0}.       \tag{4.1}
\]

Thus strict \(b<14\) means \(E_i.E_j\le13\). Real spanning by integral
vectors is rational spanning, so one can extract exactly \(\rho\) of them.

This step is automatic for exact full-Picard/full-nef return. It is not
automatic for a primitive trace block inside a larger Picard lattice, where
an ambient physical curve can merely project to the native wall line. The
projected-root scaling family remains a valid countercontrol to that weaker
formulation.

## 5. Host and complete-chamber finiteness need no new geometric gate

For the extracted root lattice \(R\), the Gram diagonal is fixed at \(-2\)
and each of the \(\binom\rho2\) off-diagonal entries has fourteen choices.
After discarding singular and wrong-signature matrices, only finitely many
root lattices remain. Full spanning and Picard integrality give

\[
                  R\subseteq\operatorname {NS}(X)\subseteq R^*. \tag{5.1}
\]

Every Picard overlattice is represented by a subgroup of the finite group
\(R^*/R\). The rank ceiling \(\rho\le20\) makes the union finite. No
discriminant bound or primitive-fibre divisibility bound is missing between
the spanning lemma and host finiteness.

For each oriented lattice, all chambers of the locally finite complete
square-minus-two mirror arrangement form one Weyl-group orbit: a generic
path crosses finitely many mirrors and successive root reflections transport
one chamber to the next. Therefore literal complete K3 root chambers also
form one orbit per host. Reflection maximality is not needed in this complete
root-return lane.

## 6. Exact nonclaims

The theorem and certificate do **not** claim any of the following:

- that every raw maximal arithmetic chamber admits the required exact K3
  carrier (Property 4);
- that an infinitely-faced convex chamber has an attained minimal aperture;
- that a local affine ADE packet alone controls fibre divisibility;
- that projected physical roots equal native trace-module roots;
- that reflection maximality is needed for the bounded packet; or
- that cusp or reflective-lattice classifications have been reproduced.

The first item is the remaining global inverse-carrier problem from the
pasted four-property program. The bounded spanning theorem is a consequence
after Properties 1--2 supply the finite-sided exact full-Picard root chamber;
it is not a construction of that carrier.

## 7. Executable pressure tests

The companion certifier checks, over exact rational arithmetic and in every
rank \(3\) through \(20\):

- finite/finite, finite/ideal, ideal/finite, and ideal/ideal pairs of distinct
  complete orthogonal packets have union rank \(\rho\);
- one packet, or a duplicated packet, has rank \(\rho-1\);
- strict normalized \(<14\) is frozen as integral \(\le13\);
- the pre-filter Gram box has size \(14^{\binom\rho2}\); and
- every hostile scope flag above remains literal.

Run the certifier with

    python3 scripts/certify_noncompact_k3_bounded_spanning_root_lemma.py \
      --check data/enumerations/noncompact_k3_bounded_spanning_root_lemma.json

and run the regression suite with

    PYTHONPATH=. pytest -q \
      tests/test_noncompact_k3_bounded_spanning_root_lemma.py

The executable checks the linear and discrete interfaces. The visual-cone
identity, the bounded-horosphere implication, and the compact-convex
two-extreme-point fact are proved symbolically in the owner report rather
than inferred from samples.
