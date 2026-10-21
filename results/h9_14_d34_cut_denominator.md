# TOOT: the variable surface-free D34 cut denominator

Status: **exact variable-class `IC0/U_weak/R/V` theorem PASS**
(2026-08-08).  Formal strong surface universality `U` is **not** claimed.

## Result

The former D34 compiler began with six archived Gram matrices.  This theorem
does not.  Its only inputs are two explicit unlabelled 14-vertex Coxeter
support templates, called `S16` and `S17`.  On either support:

- choose exactly eight support edges to have relation `m4`;
- give every remaining support edge relation `m3` and every nonedge `m2`;
- permit primitive wall norm squares only in `{2,4}`; and
- ask, rather than assume, whether the resulting integral Gram has rank and
  negative-root inertia `(10;1,9,4)`.

No Gram, surface, polytope, Vinberg transcript, or archived example is an
input.  The support edge lists are frozen in the production certificate and
are legitimate candidate-nerve/support data in the raw surface-free datum.

The exhaustive result is:

| stage | `S16` | `S17` | total |
|---|---:|---:|---:|
| raw labelled eight-`m4` assignments | 12,870 | 24,310 | **37,180** |
| raw support-automorphism orbits | 376 | 3,234 | **3,610** |
| primitive `{2,4}`-compatible labelled assignments | 3,270 | 1,519 | **4,789** |
| compatible diagram orbits | 136 | 238 | **374** |
| scaled-dual oriented integral Gram lanes | 272 | 476 | **748** |
| rank-ten Lorentzian corank-four diagram orbits | **1** | **1** | **2** |

Thus all 3,610 raw diagram orbits and all 748 integral Gram lanes have been
decided.  The unique survivor on `S16` has relation histogram
`75m2+8m3+8m4`; the unique survivor on `S17` has
`74m2+9m3+8m4`.  Their derived edge lists identify them, after the
classification rather than by input, with the known T2 and T1 diagrams.

This replaces the corpus-relative D34 statement by a genuinely variable,
surface-free denominator theorem.  It is a new classification and a new
reconstruction of the two polytopes, not a claim that T1 or T2 are previously
unknown abstract diagrams.

## The cut criterion

The main reduction is elementary but decisive.  Encode norm square `2` or
`4` at a vertex by a bit.  Crystallographic integrality gives

```text
m3 edge  => equal endpoint bits,
m4 edge  => unequal endpoint bits.
```

Consequently a primitive `{2,4}` normalization exists exactly when the `m4`
edges form a graph cut `delta(X)`.  Conversely, every such cut produces an
integral symmetric Gram: an `m3` pairing is half the common norm and an `m4`
pairing is half the larger norm.  Since each support is connected, `X` and
its complement give exactly the two scaled-dual norm orientations.

This proves exhaustiveness without searching arbitrary integer matrices.  It
also gives a reusable criterion: in the minimal one-/two-sheet alphabet,
failure of the `m4`-cut condition is an immediate arithmetic obstruction.

## Exact Gram census

The 748 exact inertia outcomes are:

| support | exact oriented-lane histogram |
|---|---|
| `S16` | `216 rank14/(2,12,0)`, `46 rank13/(2,11,1)`, `8 rank12/(2,10,2)`, **`2 rank10/(1,9,4)`** |
| `S17` | `326 rank14/(2,12,0)`, `102 rank14/(3,11,0)`, `36 rank13/(2,11,1)`, `4 rank13/(3,10,1)`, `6 rank12/(2,10,2)`, **`2 rank10/(1,9,4)`** |

The two Lorentzian lanes on each support are the scaled-dual orientations of
one diagram.  In both cases the passing `m4` cut is fixed by every support
automorphism, giving diagram automorphism orders `48` for `S16/T2` and `8`
for `S17/T1`.

## Exhaustive weak K3 packages

Only after the diagram census does the existing generic lattice compiler run.
It takes the HNF quotient of all fourteen wall roots, enumerates every even
reflection-stable fixed-lattice overclosure, embeds the distinct long-root
classes into the discriminant quadratic space of `E8(-2)`, constructs the
rank-18 graph glue, and checks primitive K3 complement glue.

| derived diagram | scaled-dual candidate counts | total | unique weak package |
|---|---:|---:|---|
| `S16/T2` | `16+5` | 21 | `(18,4,0)`, complement `U(2)+U(2)`, root-span index 1 |
| `S17/T1` | `5+2` | 7 | `(18,2,0)`, complement `U+U(2)`, root-span index 2 |

All **28** stable even overclosures are therefore part of the variable
compiler, rather than imported from the old D34 corpus.  Each survivor has
fixed lattice `U+E8(-2)`, 24 square-minus-two roots, an exact fourteen-orbit
fold, and a primitive gluing to the K3 lattice.  Generic periods in the
rank-four complement give a nonempty full-Picard weak K3 stratum.  This proves
weak surface existence; it does not prove that the folded cone is a complete
invariant nef cone.

## Fresh polytope reconstruction and Vinberg gate

The certificate next rebuilds a rank-ten root-span lattice and fourteen
inequalities directly from each generated Gram.  Exact `cddlib` enumeration
gives a strict positive interior and, for each of the fourteen walls, a
positive-square relative-interior witness with that wall pairing zero and all
other wall pairings strictly positive.  Hence the cone is full-dimensional
and every advertised facet is irredundant before the volume routine is called.

The finite-volume proof is then a literal, fresh Vinberg endpoint and complete
parabolic-neighbourhood audit:

| derived diagram | exact rays | finite | ideal | elliptic rank-8 edges | endpoint types `(0,2)/(1,1)/(2,0)` | bad edges | bad cusp neighbourhoods |
|---|---:|---:|---:|---:|---|---:|---:|
| `S16/T2` | 110 | 90 | 20 | 683 | `97/362/224` | 0 | 0 |
| `S17/T1` | 127 | 112 | 15 | 707 | `53/300/354` | 0 | 0 |

The independently enumerated finite and ideal rays agree exactly with the
Vinberg finite endpoint and full affine-neighbourhood sets.  Both generated
diagrams therefore define finite-volume, noncompact Coxeter polytopes in
`H9`.  This `V` proof does not import the old factory transcript.

## Exact universality boundary

The result discharges more than corpus-relative `IC0`, but not the formal
strong `U` of the inverse-carrier program.

- **Arithmetic `U0`: closed.**  Every object satisfying the stated raw
  D34-cut datum lies among the 3,610 diagram orbits, and every primitive
  `{2,4}` realization, scaled dual, stable even overclosure, and minimal
  `E8(-2)` lift lies in the certificate.
- **Weak-surface `U_weak`: closed.**  Each of the two Lorentzian survivors has
  a primitive full-Picard K3 lattice/period package.
- **`R`: closed at the weak marked level.**  All 24 roots and fourteen orbit
  sums fold back to the generated labelled integral Gram.
- **`V`: closed afresh.**  Cone existence, irredundancy, all edge endpoints,
  and every full parabolic neighbourhood are exact.
- **Formal strong `U/C`: not closed.**  T2 has a separate complete
  invariant-nef realization.  The unique compiled T1 package is cut by an
  exchanged-root wall, so it cannot be its complete invariant nef section.
  The present denominator excludes neither a different K3 lift alphabet nor
  another kind of algebraic carrier.

In particular, this theorem converts the previous T1 failure into a precise
criterion: the minimal `E8(-2)` carrier construction is exhaustive for the
stated D34-cut alphabet, and it succeeds strongly on T2 but fails strongly on
T1.  A broader formal `U` theorem must enlarge the carrier alphabet or shrink
the raw class by a surface-free condition that separates these two diagrams.

## Reproduction

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_d34_cut_denominator.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_h9_14_d34_cut_denominator.py
python3 -m py_compile scripts/certify_h9_14_d34_cut_denominator.py
```

Frozen artifacts:

- production script SHA-256:
  `c85b362148df3ea554b03999daa75dd866b25304c7ae4f5e13ce350dd33c98a0`;
- JSON certificate SHA-256:
  `fa89247c46f7063cba4b26904b1ff9c749647bf3f001354eb107237391cad1f9`;
- test SHA-256:
  `db512a82030cefa4875f0fb61e5dd1290bcc4643deb3fa1149fe82b89ed9e785`.

Machine-readable certificate:
[`h9_14_d34_cut_denominator.json`](../data/enumerations/h9_14_d34_cut_denominator.json).

Related scope certificates:
[`h9_14_t1_t2_d34_unfolding_census.md`](h9_14_t1_t2_d34_unfolding_census.md),
[`h9_14_s16_strong_transfer_and_s17_countercalibration.md`](h9_14_s16_strong_transfer_and_s17_countercalibration.md),
[`h9_14_t1_strong_k3_nef_obstruction.md`](h9_14_t1_strong_k3_nef_obstruction.md),
and
[`invariant_nef_exchanged_root_obstruction.md`](invariant_nef_exchanged_root_obstruction.md).
