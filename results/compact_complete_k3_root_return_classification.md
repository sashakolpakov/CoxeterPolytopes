# Exact return classification for compact complete K3 root chambers

**Status (2026-08-09): proved for the stated surface-free class.**

This is a literal

```text
compact complete Q-integral (-2)-root chamber
  -> exhaustive K3 period/nef carrier
  -> Nikulin--Vinberg surface classification
  -> exact classification of the original polytopes
```

theorem.  It removes both the rational Weyl-vector hypothesis of `C_K3W`
and any K3-embeddability hypothesis from the input.  The output is exactly
eight compact marked lattice/chamber types: six polygons in `H2` and two
polyhedra in `H3`.  Three of them, `S4`, `S6`, and `L27`, have no rational
lattice Weyl vector and are therefore genuine positive strata beyond the
previous five-type compact return theorem.

The eight lattices and polytopes are classical, not new worldwide examples.
What is new in the project is the quantifier-clean inverse theorem, the
identification of the three additional non-Weyl-vector strata, and the exact
matrix replay of all eight chambers.

## 1. A class defined without surfaces

Let `N` be an even integral lattice of signature `(1,n)`, with `n>=2`, and
choose one component of its positive cone.  Put

\[
 \Delta_2(N)=\{\delta\in N:\delta^2=-2\},\qquad
 W_2(N)=\langle s_\delta:\delta\in\Delta_2(N)\rangle .
\]

Let `C` be the closure of a component of the complement of all mirrors
`delta^perp` in the chosen positive cone.  Define `K3RootComp` to be the
class of oriented marked pairs `(N,C)`, modulo oriented integral isometry,
for which the hyperboloid projectivization

\[
                         P(N,C)=\mathbf P(C)
\]

is compact.  Equivalently, the complete square-minus-two root chamber is a
compact hyperbolic Coxeter polytope.  No surface, K3 lattice embedding,
period, nef cone, Weyl vector, finite automorphism group, or dimension bound
occurs in this definition.

There is no rootless/vacuous exception.  If `Delta_2(N)` is empty, its
chamber is the whole `H^n`, hence is not compact.  Compactness also makes the
wall list finite: the arithmetic root-mirror arrangement is locally finite,
so only finitely many mirrors meet a compact chamber.

Every facet reflection is automatically integral:

\[
              s_\delta(x)=x-\frac{2(x,\delta)}{\delta^2}\delta
                          =x+(x,\delta)\delta\in N.
\]

Thus `K3RootComp` is a non-singleton class of compact arithmetic Coxeter
polytopes of simplest type over `Q`, specified entirely by lattice and
complete-root data.

## 2. The return theorem

### Theorem 2.1 -- automatic K3 carrier and exact classification

For every `(N,C)` in `K3RootComp`:

1. `rank(N)<=4`, hence `dim P<=3`;
2. `N` has a unique primitive embedding orbit, under `O(Lambda_K3)`, in

   \[
       \Lambda_{K3}=U^3\oplus E_8(-1)^2;
   \]

3. the exact-Picard marked K3 period locus is nonempty, of complex dimension
   `20-rank(N)`;
4. after orienting a marking by the chosen positive-cone component and then
   applying its unique `W_2(N)` root-chamber re-framing, every surface in
   that locus has

   \[
                       \operatorname{Nef}(X)=C;
   \]

   so extraction returns the same normalized labelled polytope, including
   every wall; and
5. `(N,C)` is one of exactly eight marked lattice/chamber types:
   `S1,...,S6,L24,L27` in the table below.

Conversely, the exact-Picard K3 stratum attached to each of those eight
lattices is nonempty, its complete nef chamber is compact, and it returns the
listed member of `K3RootComp`.  Therefore

\[
 \boxed{
  K3RootComp=\{S1,S2,S3,S4,S5,S6,L24,L27\},\qquad \dim P\le3 .
 }
\]

This is essential surjectivity of the exact-K3-carrier return onto the
**whole stated polytope class**, not recognition of eight surfaces selected
in advance.

### Proof

The complete chamber `P` is a compact fundamental polyhedron for `W_2(N)`;
the standard reflection-arrangement argument also says that `W_2(N)` acts
simply transitively on the chambers in the chosen positive-cone component.
Hence `W_2(N)` is a cocompact lattice in `Isom(H^n)`.  It lies in the
arithmetic lattice `O^+(N)`, and a lattice subgroup of a lattice has finite
index.  The arithmetic quotient is therefore compact, so Godement's
compactness criterion makes the rational quadratic space `N_Q` anisotropic.
Meyer's theorem makes every
indefinite rational form in at least five variables isotropic.  Consequently

\[
                         \operatorname{rank}N\le4.       \tag{2.1}
\]

The stable-range primitive-embedding theorem used by Morrison now applies
without any additional discriminant gate: every even hyperbolic lattice of
rank at most ten has a primitive embedding in the K3 lattice, unique under
`O(Lambda_K3)`, and occurs as the entire Neron--Severi lattice of a
projective complex K3 surface.  A generic period in
`P(N^perp_C)` avoids all Noether--Lefschetz hyperplanes, so `NS(X)=N`
exactly and the locus has dimension `20-rank(N)`.

For a K3 surface, the nef cone is the closure of one component of the
positive cone cut by **all** square-minus-two mirrors.  Riemann--Roch chooses
an effective sign for every root, and the K3 Weyl group acts simply
transitively on the root chambers.  The roots supporting facets are the
simple effective roots and hence classes of irreducible smooth rational
curves.  Every square-minus-two reflection in `N` extends by the same
formula to an integral isometry of `Lambda_K3` and fixes `N^perp`
pointwise.  First choose the sign/component of the marking so that the ample
positive cone maps to the input's chosen component.  Composing it with an
element of `W_2(N)` then sends
the nef chamber to the prescribed `C`; that element is unique once the
initial marking is fixed.  Since the input arrangement was complete, there
are neither missing facets nor hidden root cuts.  This proves the exact
`P -> X -> same P` return.

It remains to classify `N`.  Compactness excludes a nef isotropic class,
hence excludes an elliptic fibration.  Moreover, `W_2(N)` has finite index
in `O^+(N)`, so the nef-chamber stabilizer in the integral orthogonal group
is finite; the K3 Torelli description (and the finite kernel of the action
on `NS`) gives finite automorphism group.  Nikulin's rank-three compact
classification gives exactly the six lattices `S1,...,S6`.  Vinberg's
rank-four compact classification gives exactly `L24,L27`.  Roulleau states
both exhaustions, proves geometric models, and records the complete finite
sets and Gram matrices of their irreducible `(-2)`-curves.  All full root
chambers for a fixed lattice are `W_2(N)`-conjugate, so a lattice row cannot
split into additional marked chamber types.  This gives the eight rows
below.  Conversely, those complete curve Grams have compact chambers,
closing both implications.  `square`

## 3. The eight returned polytopes

The displayed `NS` matrices use Roulleau's bases.  The exact complete
facet-wall Grams are frozen in
[`certify_compact_k3_root_return.py`](../scripts/certify_compact_k3_root_return.py).
For roots of square `-2`, intersection `0` means angle `pi/2`, intersection
`1` means angle `pi/3`, and intersection at least `3` is ultraparallel.  None
of the eight wall Grams has a parallel pair.

The displayed `f`-vector omits the terminal polytope entry: it is
`(f0,f1)` in dimension two and `(f0,f1,f2)` in dimension three.

| type | `NS` Gram | `H`-dimension | facets | `f`-vector | finite pairs `(m2,m3)` | Weyl vector? |
|---|---|---:|---:|---|---:|---|
| `S1` | `diag(6,-2,-2)` | 2 | 6 | `(6,6)` | `(6,0)` | yes |
| `S2` | `diag(36) + A2(-1)` | 2 | 6 | `(6,6)` | `(0,6)` | yes |
| `S3` | `diag(12) + A2(-1)` | 2 | 4 | `(4,4)` | `(0,4)` | yes |
| `S4` | `[[2,1,2],[1,-2,1],[2,1,-2]]` | 2 | 4 | `(4,4)` | `(0,4)` | **no** |
| `S5` | `diag(4) + A2(-1)` | 2 | 4 | `(4,4)` | `(2,2)` | yes |
| `S6` | `[[2,2,3],[2,-2,1],[3,1,-2]]` | 2 | 6 | `(6,6)` | `(2,4)` | **no** |
| `L24` | `[[2,1,1,1],[1,-2,0,0],[1,0,-2,0],[1,0,0,-2]]` | 3 | 6 | `(8,12,6)` | `(6,6)` | yes |
| `L27` | `[[12,2,0,0],[2,-2,1,0],[0,1,-2,1],[0,0,1,-2]]` | 3 | 8 | `(12,18,8)` | `(6,12)` | **no** |

The facet-wall classes span `NS` integrally in every row.  The certifier
checks an explicit unimodular change of basis from the displayed `NS` Gram
to a basis selected from the walls; equality of determinants alone is not
used.  The exact wall Grams have inertia

```text
S1,S2,S6: (1,2,3);   S3,S4,S5: (1,2,1);
L24:      (1,3,2);   L27:      (1,3,4).
```

Exact cone enumeration gives respectively

```text
6, 6, 4, 4, 4, 6, 8, 12
```

extreme rays, all timelike and none ideal or hyperideal.  Every wall is
active, so the lists are irredundant.  This is a literal compactness replay,
not an inference from the word "finite" in a surface classification.

The normalized wall Grams are pairwise distinct.  In particular `S3` and
`S4` have the same four `m=3` boundary angles but different opposite dotted
distances, so the theorem classifies metric Coxeter polytopes rather than
only nerves or angle sequences.

## 4. Why this strictly enlarges the Weyl-vector prototype

For equal norm-two walls `delta_i`, a rational lattice Weyl vector is a
solution of

\[
                         (w,\delta_i)=1\quad\text{for all }i. \tag{4.1}
\]

If `G` is the complete wall Gram, (4.1) is solvable exactly when

\[
                  \operatorname{rank}G
               =  \operatorname{rank}[G\mid\mathbf 1].       \tag{4.2}
\]

The exact rank test gives

```text
Weyl-vector rows:       S1, S2, S3, S5, L24  (five);
non-Weyl-vector rows:   S4, S6, L27           (three).
```

The corresponding vector squares are `1/2,1/4,1,2,2`, respectively, so
the five solutions are timelike, not merely formal solutions in a singular
wall-coordinate presentation.

Thus the earlier five compact `C_K3W` classes are exactly the Weyl-vector
part of this eight-class theorem.  Removing that hypothesis adds three
actual compact `P -> K3 -> P` survivors, not just hypothetical carrier
packages.

Roulleau's geometric models make the three new strata concrete:

- `S4`: a double plane branched over a smooth sextic with one tritangent line
  and one 6-tangent conic;
- `S6`: a double plane with a 6-tangent conic and two specially tangent
  cuspidal cubics; and
- `L27`: a double plane with one tritangent line and three 6-tangent conics.

Their exact-Picard moduli dimensions are respectively `17,17,16`.

## 5. Finite exhaustive compiler

For marked input `(N,C)`, the carrier list has one primitive K3 embedding
orbit: equation (2.1) puts the rank below Morrison's uniqueness cutoff.  The
remaining package is a chosen type-IV period component (together with its
conjugate when complex orientation is not fixed), its exact-Picard hyperplane
complement, and the unique `W_2(N)` chamber re-framing.  This is a finite list
of moduli packages, although each nonempty package is of course a positive-
dimensional family.

If the input is supplied only as a normalized integral norm-two wall Gram,
the marked compiler is still finite:

1. form the nondegenerate wall span;
2. enumerate the finite set of even reflection-stable integral overlattices
   between the wall span and its dual;
3. on each overlattice enumerate every norm-two mirror meeting the proposed
   compact chamber--an interior timelike vector and the finitely many exact
   vertices give a computable pairing bound, whose orthogonal negative-
   definite form gives an effective coordinate bound--and retain precisely
   the complete-root branches;
4. test the resulting rank-three/four lattice against the eight exact
   isometry classes; and
5. output its unique primitive K3 embedding orbit and period package.

This also handles saturation rather than assuming that the wall span is the
whole Picard lattice.  The classification says that every surviving branch
ends in one of the eight rows; in all eight terminal rows the walls in fact
generate `NS`.

The intrinsic nonemptiness criterion before using compactness is the usual
one: an even hyperbolic `N` has an exact-Picard K3 carrier if and only if it
embeds primitively in `Lambda_K3`; periods must avoid the rational
Noether--Lefschetz hyperplanes.  Compactness makes that criterion automatic
here by first forcing `rank(N)<=4`.  Full-root completeness is then exactly
the no-hidden-wall condition, not an extra surface guess.

## 6. What transfers back, and what does not

The theorem returns both a negative and a positive conclusion to every input
in `K3RootComp`:

```text
negative: no member exists in H4 or above;
positive: exactly six H2 and two H3 members survive;
classification: no ninth compact complete norm-two integral root chamber.
```

It is stronger than a carrier-relative rank estimate and stronger than the
five-type rational-Weyl-vector theorem.  It is not a bound for mixed-norm
integral chambers, totally-real eigen-nef sections, finite invariant folds,
or arbitrary compact arithmetic Coxeter polytopes.  Bugaenko's compact `H8`
polytope is outside this direct rational root class and remains untouched.
The theorem therefore also sharpens the design lesson: any compact surface
carrier in dimension at least four must be a proper irrational eigencone,
folded/invariant section not carrying this full rational root lattice, or a
different automorphic object--never an entire K3 norm-two nef chamber.

## 7. Sources and replay

Primary carrier, classification, and geometric sources:

- D. R. Morrison, *On K3 surfaces with large Picard number*, Invent. Math.
  75 (1984), 105--121
  ([DOI](https://doi.org/10.1007/BF01403093));
- V. V. Nikulin, *Surfaces of type K3 with a finite automorphism group and a
  Picard group of rank three*, Proc. Steklov Inst. Math. 165 (1985),
  131--155 ([bibliographic/full-text record](https://www.mathnet.ru/eng/tm2277));
- E. B. Vinberg, *Classification of 2-reflective hyperbolic lattices of rank
  4*, Trans. Moscow Math. Soc. 68 (2007), 39--66; and
- X. Roulleau, *On the geometry of K3 surfaces with finite automorphism group
  and no elliptic fibrations*, Int. J. Math. 33 (2022), 2250040
  ([DOI](https://doi.org/10.1142/S0129167X22500409),
  [arXiv](https://arxiv.org/abs/1909.01909)).

The pinned arXiv-v5 source hashes are

```text
source tar  ec64645999be6cc0e0ad3ef03fbf046a5e641f2c180e370ef2c91d052ec805a7
main TeX    effc17df3a96e4ead5fa8cee247e0fdffe178b970e6c76586592577ce9b62589
```

The source states the six-type rank-three and two-type rank-four exhaustions
and supplies the complete curve and lattice Grams.  The certifier checks the
eight exact ranks, inertias, explicit unimodular identifications of wall span
with the displayed `NS`, Coxeter entries, same-time-component extreme rays,
facet-pair incidences, `f`-vectors, irredundancy, pairwise distinction under
wall relabelling, and Weyl-vector rank tests.  It does not pretend to
recompute the Nikulin--Vinberg classification.  In particular, the finite
matrix replay proves compactness once the facet lists are accepted; the fact
that these are the complete irreducible `(-2)`-curve/facet lists is imported
from the cited classification.  Other effective roots can meet lower-
dimensional faces, but their inequalities are consequences of the simple
curve walls and do not create hidden cuts.

Run:

```bash
python3 scripts/certify_compact_k3_root_return.py --check-source
python3 -m pytest -q tests/test_compact_k3_root_return.py
```
