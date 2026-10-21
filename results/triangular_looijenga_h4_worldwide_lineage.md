# Worldwide lineage of the four triangular Looijenga `H4` types

## Result

The lineage question is settled for `T6a`, `T6b`, and `T7`:

> **None of these three is a new abstract hyperbolic Coxeter polytope.**
> What is new in the project is their realization as complete nef chambers of
> triangular Looijenga pairs.

The exact matches are:

| type | facets | exact earlier item | status |
|---|---:|---|---|
| `T6a` | 6 | Tumarkin `H4` pyramid family `07_02`, `(k,l)=(2,4)` | known abstract polytope by 2004 |
| `T6b` | 6 | Tumarkin `H4` pyramid family `07_01`, `(k,l)=(2,4)` | known abstract polytope by 2004 |
| `T7` | 7 | Im Hof doubly truncated orthoscheme item 13; independently Ma--Zheng `P_{2,1}` | known abstract polytope by 1985/1990 |
| `T8` | 8 | no exact match in the pinned ordinary or doubling CoxIter corpora | worldwide priority unresolved |

The last line is deliberately weaker.  The published eight-facet
classifications in dimension four are **compact**; `T8` is noncompact.  No
complete finite-volume, noncompact `H4` eight-facet classification was located
as of 8 August 2026.  Corpus absence therefore cannot prove novelty.

There is, however, a significant exact structural result for all four types:

> **Congruence-preimage theorem.**  Let `P4` be Vinberg's five-wall chamber for
> the full root-reflection group of `I_(1,4)`, and let
> `rho: W(P4) -> W(A4)=S5` be reduction modulo two.  Then
>
> ```text
> W(T6a) = rho^(-1)(W(A2+A1)),   index 10,
> W(T6b) = rho^(-1)(W(A2)),      index 20,
> W(T7)  = rho^(-1)(W(A1+A1)),   index 30,
> W(T8)  = rho^(-1)(W(A1)),      index 60.
> ```

Thus `T8`, even without a resolved publication-priority claim, is not an
unstructured fishing find: it is the canonical index-60 `A1` congruence cover
of Vinberg's `P4`.  The quotient `W(A4)=S5` is exactly the Weyl structure
already visible in the surface compiler through `K^perp=A4(-1)`.  This is a
clean arithmetic bridge between the algebro-geometric carrier and the Coxeter
polytope.

The replayable certificate is
[`triangular_looijenga_h4_worldwide_lineage.json`](../data/enumerations/triangular_looijenga_h4_worldwide_lineage.json).

## 1. Exact diagram matches

The comparison forgets the integral wall norm `-1` versus `-2`, normalizes the
Gram matrix to Coxeter relations, and checks all facet permutations.  This is
the correct equivalence relation for the abstract hyperbolic Coxeter polytope.

### `T6a`

In a suitable order the diagram is the path

```text
4 -- 4 -- 3 -- 4 -- infinity,
```

where `infinity` is a bold, parallel edge.  The unique exact transcription in
the pinned CoxIter snapshot is

```text
graphs/Tumarkin pyramids n+2/4-tum04_07_02-k=2-l=4.coxiter.
```

The zero-based permutation from compiler order to this source order is

```text
(0,2,3,5,4,1).
```

This is the second left-hand family in Tumarkin's `H4` pyramid table, with
`(k,l)=(2,4)`.  Tumarkin's theorem states that the displayed tables exhaust
finite-volume Coxeter pyramids over products of two simplices; together with
Kaplinskaja and Esselmann, the paper completes the finite-volume `n+2`-facet
classification.  Hence this is a complete-classification match, not a bounded
database resemblance.

### `T6b`

In compiler order its non-right relations are

```text
1--2 : 3,     1--6 : 4,
2--4 : 4,     2--5 : 4,
3--6 : infinity.
```

The unique exact Tumarkin transcription is

```text
graphs/Tumarkin pyramids n+2/4-tum04_07_01-k=2-l=4.coxiter,
```

with two diagram isomorphisms; the lexicographically first zero-based
permutation is

```text
(2,1,3,0,5,4).
```

This is the first left-hand family in Tumarkin's `H4` table at the same
parameter value `(2,4)`.  The complete `n+2` theorem again settles the abstract
lineage.

Primary source: [Tumarkin, *Hyperbolic Coxeter n-polytopes with n+2 facets*](https://arxiv.org/abs/math/0301133), submitted in 2003 and published in
*Mathematical Notes* 75 (2004), 848--854.

### `T7`

The diagram is a seven-cycle with cyclic edge labels

```text
infinity, infinity, 4, 4, 4, 4, infinity.
```

It has two independently certified earlier occurrences.

1. It is item 13 in the CoxIter transcription of Im Hof's list of doubly
   truncated orthoschemes:

   ```text
   graphs/doubly truncated orthoscheme/4-imohf85_13.coxiter.
   ```

   In the displayed Im Hof list it is the dimension-four diagram in row 3,
   column 3.  The lexicographically first exact compiler-to-source permutation
   is `(3,5,1,0,6,2,4)`.

2. Parsing all 331 rows of Ma--Zheng's complete finite-volume `H4` seven-facet
   census gives exactly one unlabelled match, `P_{2,1}`.  Their table row is

   ```text
   P_{2,1}: 220202 | 42224 | 2202 | 422 | 24 | 2,
   ```

   in the lexicographic pair convention.  The first exact permutation is
   `(1,4,5,0,3,2,6)`.  The published invariants are three cusps, arithmetic,
   and volume `pi^2/48`.

Sources:

- H.-C. Im Hof, *A class of hyperbolic Coxeter groups*, *Expositiones
  Mathematicae* 3 (1985), 179--186.
- H.-C. Im Hof, *Napier cycles and hyperbolic Coxeter groups*, *Bulletin de la
  Societe Mathematique de Belgique*, Serie A 42 (1990), 523--545; see the
  [digitized diagram list](https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/poly/napier.pdf).
- [Ma--Zheng, complete 331-polytope census](https://arxiv.org/abs/2401.13698),
  v3, 23 December 2024.

Therefore `T7` is an Im Hof orthoscheme with a new algebro-geometric carrier,
not a new abstract seven-facet chamber.

## 2. Exact covolumes and indices

For every diagram, the certificate enumerates all finite special subgroups and
uses Solomon's formula

```text
chi(W) = sum_(S finite special) (-1)^|S| / |W_S|.
```

In dimension four, Gauss--Bonnet gives

```text
Vol(P)/pi^2 = (4/3) chi(W).
```

The results are:

| chamber | `chi(W)` | `Vol/pi^2` | index in `P4` |
|---|---:|---:|---:|
| Vinberg `P4` | `1/1920` | `1/1440` | 1 |
| `T6a` | `1/192` | `1/144` | 10 |
| `T6b` | `1/96` | `1/72` | 20 |
| `T7` | `1/64` | `1/48` | 30 |
| `T8` | `1/32` | `1/24` | 60 |

The `T7` result independently reproduces Ma--Zheng's published volume.  All
candidate wall reflections are integral reflections of `I_(1,4)`, so their
groups lie in `W(P4)`.  The volume ratios are therefore actual subgroup
indices, not merely numerical coincidences.

## 3. Why the congruence-preimage statement follows

Use the standard simple roots

```text
e1-e2, e2-e3, e3-e4, e4, e0-e1-e2-e3
```

for Vinberg's `P4`.  The four square-minus-two roots form `A4`; the
square-minus-one reflection is the identity modulo two.  Exact enumeration of
the generated matrix groups over `F2` gives

| group | mod-two image | image order | image index in `S5` |
|---|---|---:|---:|
| `W(P4)` | `W(A4)=S5` | 120 | 1 |
| `W(T6a)` | `W(A2+A1)` | 12 | 10 |
| `W(T6b)` | `W(A2)` | 6 | 20 |
| `W(T7)` | `W(A1+A1)` | 4 | 30 |
| `W(T8)` | `W(A1)` | 2 | 60 |

For a subgroup `H < W(P4)`,

```text
[W(P4):H]
  = [rho(W(P4)):rho(H)] * [ker(rho):ker(rho) intersect H].
```

For each of the four types, the first factor already equals the independently
computed hyperbolic covolume index.  The second factor is therefore one, which
proves

```text
H = rho^(-1)(rho(H)).
```

This also proves arithmeticity of all four groups and replaces four isolated
examples by one exact algebraic hierarchy.

## 4. The honest boundary for `T8`

The pinned [CoxIter repository](https://github.com/rgugl/CoxIter) snapshot has
1,214 ordinary files.  Its `H4` facet histogram is

```text
facets: 5   6   7   8
files: 64  75  42   1.
```

It has the exact `T6a`, `T6b`, and `T7` matches above and no `T8` match.  The
separate 2,959-file reflection-doubling archive contains five eight-facet `H4`
records, none matching `T8`.

This is useful negative evidence but not a worldwide theorem.  The current
[HCPdm classification repository](https://github.com/GeoTopChristy/HCPdm)
states completeness for all seven-facet `H4` polytopes and compact
eight-facet `H4` polytopes.  The July 2026 paper
[Guo--Ma--Zang--Zheng](https://arxiv.org/abs/2607.14715) completes the
eight-facet census in `H5`, not the noncompact eight-facet census in `H4`.

Accordingly the correct ledger entry is:

```text
T8 abstract priority: unresolved worldwide;
T8 exact structure:   canonical index-60 A1 congruence cover of P4;
T8 AG lineage:        previously frozen full-nef Looijenga carrier.
```

## 5. What is genuinely positive

The lineage correction does not erase the positive result.  It sharpens it:

- `T6a`, `T6b`, and `T7` give genuinely new **full-nef algebraic carriers**
  for three known hyperbolic Coxeter chambers.
- The surface compiler simultaneously realizes the four standard parabolic
  congruence preimages `A2+A1`, `A2`, `A1+A1`, and `A1` inside the `A4/S5`
  quotient of `W(I_(1,4))`.
- `T8` gains a new exact arithmetic interpretation even though publication
  priority for its literal eight-wall diagram remains open.

That is a stronger and more useful algebro-geometric statement than an
unsupported claim that the underlying six- or seven-facet diagrams were new.
