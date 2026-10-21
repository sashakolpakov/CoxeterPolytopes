# Polytope-first Coble recognition: hostile `H10` audit

## Result

Starting only with the pre-existing oriented wall lattice of
`D16-D5-H10-30`, one can reconstruct its complete one-boundary Coble
decoration before invoking a surface:

```text
unique characteristic norm -1 facet b:       facet 26
norm -2 walls orthogonal to b:                12
other norm -1 walls E with b.E=1:             17
extreme rays of the b-perp root cone:         38 finite + 17 ideal + 0 hyperideal
complete cusp identity:                       {E+b} = all 17 ideal rays
root diagram:                                 Kondo Gamma_II
```

Thus the polytope itself recovers

\[
 b=-K_S,\qquad \operatorname{CM}=b^\perp,\qquad
 \{R_i\},\qquad f_E=E+b,
\]

and regenerates every one of its 30 facets as

\[
 \{b\}\ \sqcup\ \{R_1,\ldots,R_{12}\}\ \sqcup\
 \{f-b:f\text{ an ideal ray of the CM chamber}\}.
\]

Only after this arithmetic reconstruction do we use Kondō: his
classification and explicit complex Type-II example prove that the extracted
package has nonempty complex Coble moduli.  This is a genuine `P -> surface`
recognition theorem in the fixed one-boundary finite-automorphism Coble
class.  It is not a theorem that every arithmetic Coxeter polytope has a
surface.

The exact archive is
[`polytope_first_coble_h10_hostile_audit.json`](../data/enumerations/polytope_first_coble_h10_hostile_audit.json),
generated independently by
[`certify_polytope_first_coble_h10_hostile.py`](../scripts/certify_polytope_first_coble_h10_hostile.py).

## The integral uniqueness gate

An equality of real root Grams would not by itself identify the integral
decoration.  The hostile replay therefore computes the missing saturation
gate.

The twelve extracted Type-II roots span an index-two sublattice of
`b^perp=E10`.  Its Gram determinant is `-4`, with Smith diagonal

```text
(1,1,1,1,1,1,1,1,2,2),
```

so its discriminant group is `(Z/2)^2`.  There are two nonzero isotropic
order-two glue classes that give even unimodular index-two extensions.  The
exact weighted-diagram automorphism group has order 48 and exchanges those
two classes.  Consequently the `E10` saturation, and hence the integral
Type-II root embedding, is unique up to a diagram automorphism.

This closes the integral `U` gate that a graph-only comparison would leave
open.

## What is polytope-only, and where geometry enters

The following steps use only the archived oriented lattice and walls:

1. divide the common form by two to obtain the forced odd-unimodular
   normalization `I_(1,10)`;
2. primitive-normalize the walls and enumerate characteristic norm-`-1`
   facets;
3. identify `b`, prove `b^perp` is even unimodular of signature `(1,9)`, and
   split `I_(1,10)=E10+<-1>` integrally;
4. separate the walls by the numerical identities
   `b.R=0, R^2=-2` and `b.E=1, E^2=-1`;
5. enumerate every extreme ray of the root cone exactly and verify
   `{E+b}` is its complete primitive-isotropic cusp set;
6. certify the integral saturation uniqueness above.

The first non-arithmetic input is Kondō's theorem that the complex
finite-automorphism one-boundary Coble surfaces are precisely Types I and II,
both with twelve effective roots, together with his explicit existence of
Type II.  General Coble geometry then interprets `b=-K`, `B=2b`, and
`E=f-b`; it is not used to discover those vectors.

No surface-first wall list is imported by the hostile certifier.  It imports
no project scripts at all, so its replay is also independent of the main
Coble certifier and has no circular import path.

## Exhaustive frozen local scan

The same test was run on every pre-existing local `H10` chamber artifact
having explicit odd-unimodular `I_(1,10)` wall coordinates or a
reflection-stable odd-unimodular overclosure.  This is six chamber types,
24 integral presentations, and 154 characteristic norm-`-1` facet
candidates.

| chamber | integral presentations | characteristic candidates | raw Coble incidence survivors | complete cusp regeneration | decision |
|---|---:|---:|---:|---:|---|
| `P10`, archived as `R12-E8A1x2-C2-H10-12` | 1 | 1 | 1 | 1 | Coble obstruction: CM chamber has 10 roots |
| `R12-D8A1x2-C2a-H10-15` | 3 | 6 | 0 | 0 | numerical incidence obstruction |
| `R12-D8A1x2-C2b-H10-16` | 1 | 2 | 0 | 0 | numerical incidence obstruction |
| `R12-D4x2A1x2-p67-H10-38` | 15 | 120 | 0 | 0 | numerical incidence obstruction |
| `R12-D4x2A1x2-S6transposition-H10-45` | 3 | 24 | 0 | 0 | numerical incidence obstruction |
| `D16-D5-H10-30` | 1 | 1 | 1 | 1 | **nonempty complex Type-II Coble carrier** |

The early failures are exact: a putative boundary has a forbidden pairing
with some wall, rather than the required alternatives `(b.v,v^2)=(0,-2)`
or `(1,-1)`.

The other complete arithmetic survivor is the classical twelve-wall `P10`
chamber.  Its unique candidate produces ten roots and one cusp, and the
single remaining norm-`-1` facet is exactly `E=f-b`.  If this were the full
nef cone of a one-boundary Coble surface, its finite-volume CM reflection
chamber would force finite automorphism group.  Kondō's complex
classification then obstructs it because both possible one-boundary types
have twelve roots, not ten.

That is only a **one-boundary finite-automorphism Coble-carrier obstruction**.
It is not a nonexistence theorem for `P10`, which is of course a genuine
classical finite-volume Coxeter polytope, and it says nothing about carriers
by other kinds of varieties.

## Recognition criterion in this class

For an oriented odd-unimodular `H10` wall datum `(L;v_i)`, the following is
necessary and sufficient for a complex finite-automorphism one-boundary
Coble carrier:

1. there is a characteristic norm-`-1` facet `b`;
2. every other facet is either a root `(v^2,b.v)=(-2,0)` or an exceptional
   wall `(-1,1)`;
3. `b^perp` is `E10` and the extracted root cone has finite volume;
4. the exceptional walls are exactly `f-b` for all primitive isotropic cusp
   rays `f` of that root cone;
5. the complete integral root decoration, including its saturation, is
   Kondō Type I or II.

Necessity is the divisor/adjunction/Halphen dictionary.  Sufficiency at the
last step is the nonemptiness of Kondō's corresponding surface type.  Outside
this explicitly stated class the test is only a filter.

Primary geometric source: Shigeyuki Kondō,
[*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591),
especially the Main Theorem/Table 2, Example Type II, and the results labelled
`CME10`, `NefCone`, `Vinberg2`, `Halphen`, and `fibration-isotropic`.

## Replay

```bash
PYTHONPATH=src:. python3 scripts/certify_polytope_first_coble_h10_hostile.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_polytope_first_coble_h10_hostile.py
```
