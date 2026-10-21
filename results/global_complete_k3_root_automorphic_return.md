# Global inverse return for complete integral norm-two chambers

**Status (2026-08-09): proved exact inverse prototype/closed denominator,
by synthesis of published classification theorems, with a new project-level
quantifier cleanup and cross-classification.**

This is the strongest exact noncircular complete-root inverse prototype
currently available in the project.  It starts with no surface and no
K3-embedding hypothesis.  Its input is the whole finite-volume class of
complete square-minus-two reflection chambers in even integral hyperbolic
lattices.  Nikulin--Vinberg classification forces every input lattice onto
the published 118-element list; global Torelli then supplies an exact-Picard
K3 carrier, and the K3 ample-cone theorem returns the same chamber.  Haowu
Wang's independent classification of automorphic corrections gives an
additional exact split: in hyperbolic dimensions 5 through 18, 51 marked
chamber types have a Borcherds/2-reflective automorphic correction and 18 do
not.

This is not yet the flagship AG-leverage theorem.  The `118` classification
and the bound `n<=18` enter first from the Nikulin--Vinberg classification of
the Coxeter/lattice input itself; they are not independently derived at the
last arrow from a softer algebraic-geometric boundedness theorem.  The K3,
nef-cone, automorphism, and Cox-ring equivalences are genuine structural
information, but they do not sharpen that already-used classification.

The 118 lattices and their Coxeter chambers, and Wang's `51+18` theorem, are
published.  They are not new abstract polytopes.  New project content is the
surface-free inverse formulation with K3 embeddability proved as an output,
the exact `P -> finite list of K3 carrier-component types -> same P`
quantifiers, and the intersection of Wang's obstruction list with the
already certified nineteen rank-11--19 rows (`7` corrected and `12`
uncorrected).  The word *finite* applies to discrete primitive-embedding,
orientation, and moduli-component package types; each such type can contain
a positive-dimensional, hence infinite, family of K3 period points.

## 1. The surface-free class

Let `Root2FV` consist of oriented marked pairs `(S,C)`, modulo oriented
integral isometry, with the following properties.

1. `S` is an arbitrary even integral lattice of signature `(1,n)`, `n>=2`.
2. With

   \[
   \Delta_2(S)=\{r\in S:r^2=-2\},\qquad
   W_2(S)=\langle s_r:r\in\Delta_2(S)\rangle,
   \]

   `C` is the closure of a component of the chosen positive cone after all
   mirrors `r^perp` are removed.
3. The hyperboloid projectivization `P(C)` has finite hyperbolic volume.

No algebraic surface, K3 lattice, primitive embedding, period, nef cone,
automorphism group, Cox ring, Weyl vector, or automorphic form occurs in this
definition.  Equivalently, condition 3 says that `W_2(S)` has finite index in
`O^+(S)`, but finite volume is retained as a literal Coxeter/Gram property of
the input chamber.

The compact and noncompact finite-volume subloci are denoted `Root2Comp` and
`Root2FV^nc`; they are kept separate throughout.

## 2. Global classification and automatic K3 carrier

### Theorem 2.1 -- 118-type exact return theorem

Every `(S,C)` in `Root2FV` has a nonempty finite exhaustive list of marked
exact-Picard K3 primitive-embedding/orientation/moduli-component package
types, and every type returns the same marked integral Coxeter chamber `C`.
The period locus inside a package type is generally continuous and may be
positive-dimensional; the theorem does not claim that there are only
finitely many K3 surfaces.  The rank distribution is

| `rank S` | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | `>=20` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| types | 26 | 14 | 9 | 10 | 9 | 12 | 10 | 9 | 4 | 4 | 3 | 3 | 1 | 1 | 1 | 1 | 1 | 0 |

Consequently

\[
 \boxed{|\mathrm{Root2FV}|=118,\qquad \dim P\le18.}
\]

Exactly eight marked types are compact: six in `H2` and two in `H3`.
The other 110 types are finite-volume noncompact.  In particular there is no
finite-volume complete integral norm-two chamber in `H19` or above.

This counts marked integral lattice/chamber types.  It does not assert,
without a separate forgetful-map audit, that all 118 remain distinct after
forgetting the integral lattice and retaining only the unmarked real metric
polytope.

More intrinsically, let `S` be any even integral lattice of signature
`(1,n)`, `n>=2`, and let `C` be one chamber of its complete `(-2)`-root
arrangement.  Then the following are equivalent:

1. `P(C)` has finite hyperbolic volume;
2. after the harmless global sign change used by the cited sources, `S` is
   one of the 118 Nikulin--Vinberg lattices;
3. there exists an exact-Picard projective K3 surface `X` with `NS(X)=S`,
   and every such marked K3 has finite automorphism group;
4. the nef/effective cone of such an `X` is rational polyhedral; and
5. `Cox(X)` is finitely generated.

The implications between 3--5 use genuinely algebraic input: global Torelli
and Sterk's cone theorem, together with Artebani--Hausen--Laface, Theorems
2.7 and 2.11.  The equivalence with 1 and the list in 2 is the published
Nikulin--Vinberg arithmetic/AG classification itself.  It is not presented
here as a new Chern-inequality shortcut or a newly sharpened worldwide
dimension bound.

### Proof

For a complete root arrangement, `W_2(S)` acts simply transitively on the
chambers in a chosen positive-cone component.  Thus finite volume of `C` is
equivalent to finite index of `W_2(S)` in the arithmetic lattice `O^+(S)`:
`S` is a hyperbolic 2-reflective even lattice in the terminology of Nikulin,
Vinberg, Gritsenko--Nikulin, and Wang.

Nikulin classified these lattices in every rank except four, and Vinberg
classified rank four.  Wang's Table 3 restates the complete rank counts above
and records that the models occur in Gritsenko--Nikulin, Section 3.2 (with the
published rank-ten correction `U+D4+4A1`).  The counts sum to 118 and vanish
from rank 20 onward.

The same published classification is the K3 finite-automorphism
classification.  Gritsenko--Nikulin explicitly state, by global Torelli, that
the elliptically 2-reflective lattices on this list give all Picard lattices
of complex K3 surfaces of Picard rank at least three and finite automorphism
group.  Roulleau's atlas gives geometric K3 models and the complete
configurations of smooth rational curves for all 118 lattices.  Therefore a
primitive embedding `S -> Lambda_K3` is a theorem output for every input; it
is not an admissibility condition smuggled into `Root2FV`.

For completeness, the carrier compiler for a marked input has finite
discrete output at the level of primitive-embedding, orientation, and
moduli-component package types.  It enumerates primitive embeddings of the
fixed lattice `S` into `Lambda_K3` through complement genera, complement
isometry classes, and finite discriminant anti-isometry double cosets,
retaining orientations and the marked chamber.  Published classification
makes this list nonempty.  It does not enumerate individual K3 surfaces:
for each embedding type, the complement contains a generally continuous
exact-Picard period locus obtained by avoiding all Noether--Lefschetz
hyperplanes.  Choosing any point of that locus gives a projective K3 surface
with `NS(X)=S` exactly.  The K3 Weyl group acts simply transitively on the
complete `(-2)`-root chambers, so one Weyl re-framing sends the nef cone to
`C`.
Because the input already used every square-minus-two mirror, there is no
hidden-root alternative.  Hence

\[
             \mathbf P(\operatorname{Nef}X)=\mathbf P(C)=P
\]

with the same saturated lattice, normalized wall Gram, labels, face
arrangement, and all walls.  Conversely every exact-Picard finite-automorphism
K3 surface gives a member of `Root2FV`.  This proves both essential
surjectivity and exact return.

The compact sublocus is the independent Nikulin rank-three/Vinberg rank-four
classification already replayed in
[`compact_complete_k3_root_return_classification.md`](compact_complete_k3_root_return_classification.md):
`S1,...,S6,L24,L27`.  Every rank at least five rational quadratic space is
isotropic by Meyer, so no further member can be compact.  Thus the remaining
`118-8=110` are finite-volume noncompact.

## 3. Borcherds' idea becomes an exact 51/18 criterion

For a hyperbolic 2-reflective lattice `S`, Wang calls a positive integer `m`
an **automorphic correction** when `U(m)+S` admits a 2-reflective modular
form.  The associated orthogonal modular variety carries a reflective
divisor, and its Fourier expansion at the chosen cusp supplies the
Borcherds/Lorentzian Kac--Moody correction of the same norm-two Weyl group.

### Theorem 3.1 -- exact automorphic-correction split

Among the 69 members of `Root2FV` of rank at least six, equivalently in
hyperbolic dimensions `5,...,18`, exactly

\[
                         \boxed{51+18=69}
\]

have, respectively do not have, an automorphic correction.  The 18 negative
lattices (negative-definite root summands are written without signs) are

```text
U+E8+E7,       U+E8+D6,       U+E8+D4+A1,
U+E8+D4,       U+D8+D4,       U+E8+4A1,
U+E8+3A1,      U+D8+3A1,      U+E8+A3,
U+D8+2A1,      U+2D4+2A1,     U+E8+A2,
U+E6+A2,       U+D4+A3,       U+D5+A2,
U+D4+A2,       U+A4+A2,       U+A3+A2.
```

Every other rank-at-least-six member has at least one correction.  This is
necessary and sufficient, not a search through a bounded modular-form
weight.  It is Wang's Theorem 7.2, proved from the classification of
2-reflective modular forms using Jacobi forms and explicit Borcherds
quasi-pullbacks.

This sharpens the role of vector-valued modular forms in the inverse program:
they do not automatically correct every finite-covolume reflection group.
They give an exact additional algebraic structure on 51 types, while the
remaining 18 are published counterexamples to automatic correction.  The K3
carrier and the same-polytope return remain valid on both sides of the split.

Wang's theorem deliberately starts at signature `(1,n)` with `n>=5`.
Nothing here classifies automorphic corrections for the 49 rank-three,
rank-four, and rank-five types.

## 4. Exact intersection with the nineteen high-rank rows

Intersecting Wang's 18-element obstruction list with the independently
certified rank-11--19 K3 denominator gives the following new project-level
ledger.

| `H` dimension | lattice `S` | facets | correction? |
|---:|---|---:|---|
| 10 | `U+E8+A1` | 12 | yes |
| 10 | `U+D8+A1` | 14 | yes |
| 10 | `U+2D4+A1` | 22 | yes |
| 10 | `U+D4+5A1` | 90 | yes |
| 11 | `U+E8+A2` | 13 | **no** |
| 11 | `U+E8+2A1` | 14 | yes |
| 11 | `U+D8+2A1` | 19 | **no** |
| 11 | `U+2D4+2A1` | 59 | **no** |
| 12 | `U+E8+A3` | 14 | **no** |
| 12 | `U+E8+3A1` | 17 | **no** |
| 12 | `U+D8+3A1` | 39 | **no** |
| 13 | `U+E8+D4` | 15 | **no** |
| 13 | `U+D8+D4` | 20 | **no** |
| 13 | `U+E8+4A1` | 27 | **no** |
| 14 | `U+E8+D4+A1` | 21 | **no** |
| 15 | `U+E8+D6` | 19 | **no** |
| 16 | `U+E8+E7` | 19 | **no** |
| 17 | `U+2E8` | 19 | yes |
| 18 | `U+2E8+A1` | 24 | yes |

Thus exactly seven of the nineteen have a correction and twelve do not.  The
corrected high-rank dimensions are `H10`, `H11`, `H17`, and `H18`; the whole
published strip `H12`--`H16` is uncorrected in this complete-nef denominator.
This is a positive classification of extra automorphic structures on known
polytopes, not a claim of seven new abstract polytopes.

## 5. What has and has not been achieved

This theorem closes essential surjectivity and exact return on a broad
non-singleton surface-free class:

```text
arbitrary P in Root2FV
  -> finite exhaustive list of K3 embedding/orientation/moduli-component types
     (nonempty by classification; continuous period families lie inside)
  -> exact full-nef return of the same P
  -> reuse of the prior lattice classification: 118 marked types, dimension <=18
  -> exact automorphic perk in H5--H18: 51 corrected / 18 uncorrected.
```

It strictly strengthens the previous high-rank theorem by removing primitive
K3 embeddability from the input and covering every rank from three through
nineteen.  It also answers the Borcherds side question exactly on the
rank-at-least-six part of this class.

It does **not** complete the non-negotiable flagship chain.  Essential
surjectivity here is proved by first invoking the exhaustive
Nikulin--Vinberg classification of `P`/lattices and only then interpreting
every entry as a K3 Picard lattice.  Accordingly, the numerical conclusion
`118` and the boundary `n<=18` have already been supplied on the input side;
Sterk and Artebani--Hausen--Laface give a genuine native-AG dictionary but no
new terminal theorem about `P`.  This result is therefore the exact closed
denominator against which a broader mixed-root or totally-real theorem must
be measured.  The flagship remains to find an exhaustive `P ->` carrier
transition for such a broader class and use independent AG
classification/boundedness to move a polytope boundary or construct new
survivors.

It does **not** classify mixed root lengths, incomplete reflection systems,
folded or invariant chambers, totally-real eigensections, or arbitrary
compact arithmetic Coxeter polytopes.  In particular Bugaenko's compact
`H8` chamber is an irrational golden eigensection of a rank-18 K3 divisor
lattice, not a direct rational complete norm-two chamber, and lies outside
`Root2FV`.  Therefore the global compact example frontier and the universal
compact nonexistence boundary do not move.

The dimension-18 bound and the 118 lattices are classical classification
results, not new worldwide theorems.  The honest advance is an exact
inverse-carrier prototype and a positive/negative automorphic criterion in
the largest presently classifiable direct-root denominator, not a new
AG-derived polytope boundary.

## 6. Sources and replay

Primary sources:

- H. Wang, [*The classification of 2-reflective modular forms*](https://doi.org/10.4171/JEMS/1358), especially Table 3 and Theorem 7.2 ([arXiv](https://arxiv.org/abs/1906.10459));
- V. Gritsenko and V. Nikulin, [*Lorentzian Kac--Moody algebras with Weyl groups of 2-reflections*](https://doi.org/10.1112/plms.12084), especially Section 3.2 and the global-Torelli K3 interpretation ([arXiv](https://arxiv.org/abs/1602.08359));
- X. Roulleau, [*An atlas of K3 surfaces with finite automorphism group*](https://doi.org/10.46298/epiga.2022.6286), for the 118 K3 families and their complete curve configurations ([arXiv](https://arxiv.org/abs/2003.08985)); and
- the original Nikulin--Vinberg hyperbolic 2-reflective classifications cited in those papers.

For the algebraic equivalences, see H. Sterk, *Finiteness results for
algebraic K3 surfaces*, and M. Artebani, J. Hausen, and A. Laface,
[*On Cox rings of K3 surfaces*](https://doi.org/10.1112/S0010437X09004576),
Theorems 2.7 and 2.11.

The exact rank arithmetic and the `7/12` high-rank intersection are replayed
by
[`certify_global_k3_root_automorphic_return.py`](../scripts/certify_global_k3_root_automorphic_return.py)
and
[`test_global_k3_root_automorphic_return.py`](../tests/test_global_k3_root_automorphic_return.py).

The gate-by-gate source and scope check is frozen in the
[`hostile audit`](global_complete_k3_root_automorphic_return_hostile_audit.md).
