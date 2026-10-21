# Oriented markings and arithmetic period groups for the 19 high-rank K3 root chambers

**Frozen scope:** the nineteen finite-volume rows of `K3RootFV[11,20]`.  
**Result:** the positive-cone, primitive-embedding, complete-chamber, and
oriented arithmetic period/change-of-marking image questions are uniform and
exact for all nineteen rows.  
**Boundary:** no presentation of an infinite orthogonal group is claimed, and
surface automorphisms, diagram symmetries, and geometric family monodromy are
not relabelled as one another.

## The oriented arithmetic theorem

Fix one of the nineteen lattices (N), a primitive embedding
(j:N\hookrightarrow\Lambda_{K3}), its orthogonal complement (T=j(N)^\perp),
a positive-cone component (P_N), and a complete square-minus-two chamber
(C\subset P_N).  Write

\[
\widetilde O^+(L)=\ker\bigl(O^+(L)\longrightarrow O(q_L)\bigr).
\]

Then:

1. The primitive embeddings of (N) in (Lambda_{K3}) form one orbit under
   (O(\Lambda_{K3})), and still form **one orbit under**
   (O^+(\Lambda_{K3})).  Thus passing to the oriented K3 lattice does not
   double any of the nineteen embedding types.
2. The two components of the positive cone of (N) form one (O(N))-orbit
   and two (O^+(N))-orbits.  This is a statement about a fixed lattice.  If
   one regards a chosen component as part of an oriented object and allows an
   isometry to carry chosen component to chosen component, there is only one
   oriented integral isomorphism type.
3. In a fixed component, (W_2(N)) acts simply transitively on the complete
   root chambers.  Consequently the raw fixed-marking chamber set is an
   infinite (W_2(N))-torsor, but there is one chamber orbit under (W_2(N))
   and hence one under (O^+(N)).
4. Both discriminant representations are onto:

   \[
   O^+(N)\twoheadrightarrow O(q_N),\qquad
   O^+(T)\twoheadrightarrow O(q_T).
   \]

5. For **pointwise (N)-polarization**, let

   \[
   \Gamma^{+,\mathrm{pt}}_j=
   \{g\in O^+(\Lambda_{K3}):g|_{jN}=\mathrm{id}\}.
   \]

   Its restriction to `T`, hence its oriented arithmetic period/change-of-
   marking **image** on a chosen type-IV component, is exactly

   \[
   \operatorname{res}_T(\Gamma_j^{+,\mathrm{pt}})=\widetilde O^+(T),
   \]

   so its discriminant image is trivial.  This is the component-preserving
   version of Dolgachev's arithmetic change-of-marking calculation.  Calling
   it the geometric monodromy of a family is conditional on the standard
   moduli-stack/deck-group identification; that identification is not
   independently certified here.
6. If the marking of (N) is forgotten but its isometry type and nef chamber
   are retained, the compatible ambient normalizer is

   \[
   \Gamma^+_{j,C}=\{g\in O^+(\Lambda_{K3}):g(jN)=jN,\ g(jC)=jC\}.
   \]

   Restriction to (N) gives the exact sequence

   \[
   1\longrightarrow\widetilde O^+(T)
   \longrightarrow\Gamma^+_{j,C}
   \longrightarrow\operatorname{Stab}_{O^+(N)}(C)
   \longrightarrow1.
   \]

   Its action on `T` has image all of `O^+(T)`, and its common discriminant
   image is all of `O(q_T)\simeq O(q_N)`.  Thus `O^+(T)` is the
   **T/period image of** `Gamma^+_{j,C}`; it is not an equality
   `Gamma^+_{j,C}=O^+(T)`.

This distinguishes two arithmetic images that the unmarked classification
report had left bundled together: the pointwise-polarized T-image is the
stable kernel, whereas the unmarked lattice/chamber normalizer has the full
component-preserving orthogonal T-image.

## Why orientation does not split the embeddings

Every displayed (N) visibly splits a unimodular (U).  Roulleau's proof
also gives a `U` in every `T`: for ranks 11 through 18, the lattice
is embedded in (U^2\oplus E_8(-1)^2), leaving the third (U); the rank-19
case is checked separately.  A copy of the unimodular (U) is primitive and
therefore an orthogonal summand.

For a basis (e,f) of (U), let (h=e+f), so (h^2=2).  Its reflection is

\[
s_h=\begin{pmatrix}0&-1\\-1&0\end{pmatrix}.
\]

It is integral, has determinant (-1), reverses one positive direction, and
has trivial discriminant action because (U) is unimodular.  Applied in (T),
it extends while fixing (N) pointwise.  Hence a full-(O(\Lambda_{K3}))
transporter with the wrong positive-three-space orientation can be corrected
without changing the embedded (N).  This proves the (O^+)-orbit statement.

The same reflection in the visible (U\subset N) shows that (O(N)) exchanges
the two positive-cone components and that restricting a discriminant lift to
(O^+(N)) loses no discriminant image.  Applied in (T), it gives the
corresponding orientation correction for the type-IV components.

## Why the discriminant images are exact

Nikulin's Theorem 1.14.2 implies that an even indefinite lattice splitting
(U) maps onto the orthogonal group of its discriminant form.  We use this
strong, easily checked special case; none of the delicate weaker local
hypotheses is being inferred merely from a rank inequality.  The explicit
positive reflection above corrects the orientation while acting trivially on
the discriminant module, so surjectivity remains true after replacing (O)
by (O^+).

Every square-minus-two reflection acts trivially on (A_N): for
(x\in N^\vee),

\[
s_r(x)-x=(x,r)r\in N.
\]

Since (W_2(N)) is simply transitive on chambers, every element of (O^+(N))
can be moved into the chamber stabilizer by a unique Weyl element without
changing its discriminant action.  Therefore

\[
\operatorname{Stab}_{O^+(N)}(C)\twoheadrightarrow O(q_N).
\]

The primitive unimodular gluing anti-isometry (q_T\simeq-q_N) then gives the
normalizer exact sequence above.

## The nineteen finite discriminant quotients

The certifier reconstructs `N=U+R(-1)` and independently computes its Smith
form.  The actual ADE discriminant quadratic forms are taken from the standard
ADE table encoded in the certifier; they are **not** independently derived
from inverse Gram matrices.  Conditional on those table entries, the
certifier enumerates the full finite orthogonal groups.  The frozen expected
orders are regression values, not a second derivation.  The last column is
simultaneously

\[
[O^+(T):\widetilde O^+(T)]
=[O^+(N):\widetilde O^+(N)].
\]

| case | `A_N` invariant factors | order of `O(q_N)` |
|---|---:|---:|
| `R11_UE8A1` | `2` | 1 |
| `R11_UD8A1` | `2,2,2` | 2 |
| `R11_UD4D4A1` | `2^5` | 72 |
| `R11_UD4A1x5` | `2^7` | 40,320 |
| `R12_UA2E8` | `3` | 2 |
| `R12_UE8A1x2` | `2^2` | 2 |
| `R12_UD8A1x2` | `2^4` | 12 |
| `R12_UD4x2A1x2` | `2^6` | 1,440 |
| `R13_UE8A3` | `4` | 2 |
| `R13_UE8A1x3` | `2^3` | 6 |
| `R13_UD8A1x3` | `2^5` | 120 |
| `R14_UE8D4` | `2^2` | 6 |
| `R14_UD8D4` | `2^4` | 120 |
| `R14_UE8A1x4` | `2^4` | 24 |
| `R15_UE8D4A1` | `2^3` | 6 |
| `R16_UE8D6` | `2^2` | 2 |
| `R17_UE8E7` | `2` | 1 |
| `R18_UE8E8` | trivial | 1 |
| `R19_UE8E8A1` | `2` | 1 |

These are arithmetic group indices.  They are **not** asserted to be degrees
of coarse period maps: central isometries, notably (-1), can act trivially on
the projectivized period domain.

## Four groups that must remain separate

- A weighted Coxeter-diagram automorphism is only a permutation preserving the
  displayed wall data; it need not extend to the lattice.
- The integral chamber stabilizer `Stab_{O^+(N)}(C)` is a lattice group.  It
  maps onto `O(q_N)`, but it is not thereby a surface automorphism group.
- The automorphism group of a particular K3 surface is the subgroup satisfying
  the nef and Hodge-extension conditions.  Roulleau's generic groups are
  recorded only as comparison labels.
- The pointwise lattice-polarized arithmetic T/period image is
  `tilde O^+(T)`; the unmarked chamber normalizer's T/period image is
  `O^+(T)`.  Neither formulation independently identifies geometric family
  monodromy.

## Exact residual boundary

No generators-and-relations presentation for any of the nineteen infinite
groups `tilde O^+(T)` has been derived from the local pinned sources,
Nikulin's bibliographic input, or the table-driven finite discriminant
enumeration.  Giving such presentations would require explicit complement
Gram bases followed by a separate arithmetic-group presentation computation.
The present theorem classifies the arithmetic T-images and their finite
discriminant images exactly, but it does not disguise a group name as a finite
presentation or an arithmetic change-of-marking group as independently
certified geometric monodromy.

The separate [rank-11--14 literal Vinberg theorem](k3_rank11_14_literal_vinberg_volume.md)
now combines with the five rank-15--19 audits to close literal endpoint and
full parabolic-neighbourhood `V` for **all 19** rows.  There is no remaining
literal-V gap here.

Strict M3c nevertheless remains open for two substantive packages: the
ambient class is not the bounded raw `IC0` specification, and there are no
explicit presentations of the nineteen infinite arithmetic groups or an
independently certified geometric-monodromy packet.  The present work closes
the orientation/embedding/chamber and arithmetic discriminant-image caveat,
not those residual requirements.

## Replayable evidence

- [rowwise manifest](../data/enumerations/k3_high_rank_oriented_monodromy.json)
- [certifier](../scripts/certify_k3_high_rank_oriented_monodromy.py)
- [tests](../tests/test_k3_high_rank_oriented_monodromy.py)
- [upstream 19-stock inventory](../data/enumerations/roulleau_rank11_14_source_stock_inventory.json)
- [upstream rank-15--19 chambers](../data/known_examples/k3_high_rank_nef_chambers.json)
- [upstream rank-15--19 literal audits](../data/known_examples/k3_high_rank_nef_vinberg_volume.json)
- [upstream rank-11--14 literal theorem](k3_rank11_14_literal_vinberg_volume.md)
- [downstream 19-type synthesis, context only](k3_high_rank_full_root_classification_transfer.md)

The machine dependency graph is deliberately acyclic: the certifier neither
reads nor hashes the downstream 19-type synthesis.  It rebuilds from the
Roulleau 19-stock inventory, the two rank-15--19 chamber/volume archives, the
rank-11--14 literal archive, the embedded exact ADE factor table, and the
locally checked Roulleau/Dolgachev sources.

Source inputs and provenance:

- V. V. Nikulin, *Integral symmetric bilinear forms and some of their
  applications*, Math. USSR-Izv. 14 (1980), Theorem 1.14.2,
  DOI `10.1070/IM1980v014n01ABEH001060`.  This is a bibliographic/primary-web
  input at `https://www.mathnet.ru/eng/im1677`, **not** a locally hash-pinned
  source;
- I. Dolgachev, *Mirror symmetry for lattice polarized K3 surfaces*,
  Proposition 3.3, arXiv `alg-geom/9502005`, pinned source SHA-256
  `4bc625d0854b5c1c304d098e7247fd54597e0c781ea48c4320aa0d9f8b091baf`;
- X. Roulleau's pinned TeX, SHA-256
  `e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3`,
  especially lines 1142--1168, 1184--1195, and 1221--1227.

Run:

```bash
python3 scripts/certify_k3_high_rank_oriented_monodromy.py --check-archive --check-sources
PYTHONPATH=src:. pytest -q tests/test_k3_high_rank_oriented_monodromy.py
```
