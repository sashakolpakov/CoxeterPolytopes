# TOOT: the raw `H14/60` chamber uniquely reconstructs the Type-VI Coble surface package

**Status (2026-08-09): exact `POLYTOPE -> SURFACE` recognizer for the
five-boundary Type-VI row; broader nine-class decision architecture stated,
but template-free essential surjectivity remains open.**

The inverse calculation starts with only an even arithmetic lattice `M` and
sixty oriented primitive wall vectors.  It does **not** read a surface, a
Néron--Severi marking, boundary labels, Coble--Mukai labels, or a stored cusp
lift table.  Nevertheless it recovers, uniquely,

```text
M and 60 raw mirrors
  -> one five-boundary de-scaling N=I_(1,14)
  -> K=-1/2 sum B_i, characteristic, K^2=-5
  -> 5 B(-4) + 10 R(-2) + 10 crossing E(-1) + 35 cusp E(-1)
  -> CM=E10 and the labelled 20-root Gamma_VI chamber
  -> every one of 57 CM cusps and all 285 cusp/boundary lift tests
  -> the identical set of 60 oriented mirrors
  -> Kondō Type VI, characteristic 3, one surface isomorphism class.
```

The exact archive is
[`polytope_first_coble_type_vi_h14.json`](../data/enumerations/polytope_first_coble_type_vi_h14.json),
generated independently by
[`certify_polytope_first_coble_type_vi_h14.py`](../scripts/certify_polytope_first_coble_type_vi_h14.py).
The companion forward construction remains
[`coble_type_vi_char3_five_boundary_full_ns_h14.md`](coble_type_vi_char3_five_boundary_full_ns_h14.md).

This is the first high-dimensional Coble row in the project for which the
boundary frame is discovered from the arithmetic Coxeter datum rather than
carried over from a surface construction.  It is not a new abstract
polytope: the chamber is Alexeev's known `(15,7,1)` chamber.

## 1. Raw input and the normalization that matters

Let `(M,q_M; r_1,...,r_60)` be the archived raw datum.  Exact arithmetic gives

\[
 \operatorname{sig}(M)=(1,14),\qquad |\det M|=2^7,
 \qquad (r,a,\delta)=(15,7,1).
\]

The primitive walls have square/divisibility distribution

```text
15 * (-2,div 1) + 35 * (-2,div 2) + 10 * (-4,div 2).
```

Every reflection is an integral involution of `M`, the complete wall Gram is
acute, and every pair has a literal Coxeter, parallel, or ultraparallel
relation.  The fifteen `(-2,div 1)` roots form a `Z`-basis of `M`.  Call them
the *main roots*.  Nothing geometric has entered at this point.

The factor of two cannot be ignored.  The intended Picard form is not
`q_M`; it is `q_N=q_M/2` on a proper index-sixteen sublattice `N`.  A boundary
wall is represented arithmetically by its half-class `h_i=B_i/2`, so

\[
 q_M(h_i)=-2,qquad q_N(B_i)=q_N(2h_i)=-4.
\]

This explains why simply calling all short roots `(-2)` curves would destroy
the divisor interpretation.

## 2. Exhaustive boundary-frame theorem

For every five-subset `H={h_1,...,h_5}` of the fifteen main roots, require
the `h_i` to be pairwise orthogonal and define, entirely inside `M`,

\[
 K_H=-\sum_{i=1}^5h_i,
 \qquad
 N_H=\left\langle
       \{a:a\text{ is a main root outside }H\},
       2h_1,\ldots,2h_5,K_H
      \right\rangle .                                    \tag{2.1}
\]

Put `q_H=q_M/2` on `N_H`.  The recognizer applies the following gates in
order:

1. `q_H` is integral;
2. `(N_H,q_H)` is odd unimodular of signature `(1,14)`;
3. `K_H^2=-5` and `K_H` is characteristic;
4. every raw wall descends to `N_H`, replacing `h_i` by the physical divisor
   `2h_i` on the five selected mirrors;
5. every non-boundary wall satisfies adjunction
   `K_H.C=-2-C^2` and has one of the allowed boundary profiles.

The exhaustive counts are

| prefix gate | survivors |
|---|---:|
| all five-subsets of the main roots | 3,003 |
| pairwise orthogonal | 338 |
| `q_M/2` integral on `N_H` | **1** |
| odd unimodular, characteristic `K`, all walls descend | **1** |
| complete adjunction alphabet | **1** |

Thus the five-boundary package is unique as an unordered frame.  Its
existence is not inferred from the visible `K5^(2)` graph: 338 orthogonal
frames survive that weak combinatorial condition.  The decisive gate is the
integral de-scaling.

For the unique frame,

\[
 [M:N]=2^4,
 \qquad
 \det(q_N)=\frac{\det(q_M)[M:N]^2}{2^{15}}
            =\frac{2^7 2^8}{2^{15}}=1.                  \tag{2.2}
\]

The computed lattice is therefore `I_(1,14)`.  Moreover

\[
 -2K=B_1+\cdots+B_5,qquad B_i^2=-4,qquad B_iB_j=0.
\]

Adjoining the five half-boundaries to `N` recovers `M` exactly, and doubling
the Picard pairing recovers the raw form entry by entry.  This closes both
the scaling and saturation gates.

### Necessity of the search denominator

On a five-boundary Type-VI surface, each `h_i=B_i/2` is a short root of
divisibility one because it pairs oddly with the crossing exceptional
curves.  The fifteen short divisibility-one roots consist of the five
boundary halves and ten crossing exceptional curves.  Hence every genuine
Type-VI boundary frame occurs among the 3,003 enumerated subsets.  The search
is exhaustive, not a lucky choice of five vertices.

## 3. The divisor alphabet is forced

With `N`, `K`, and the `B_i` recovered, each raw mirror receives its role
from its square and its five boundary intersections:

| role | numerical test in `N` | count |
|---|---|---:|
| boundary | `B_i^2=-4`, other boundary pairings zero | 5 |
| actual `(-2)` curve `R` | `R^2=-2`, `R.B_i=0` for all `i` | 10 |
| crossing `(-1)` curve `E` | `E^2=-1`, profile `(1,1,0,0,0)` | 10 |
| single-boundary `(-1)` curve `E` | `E^2=-1`, profile `(2,0,0,0,0)` | 35 |

All sixty satisfy `K.C=-2-C^2`.  This gives the exact

```text
5 + 10 + 10 + 35 = 60
```

partition without reading the forward certificate's labels.

## 4. Fold to `CM`, including saturation

Regard `M` as the doubled extended Picard module

\[
 \widetilde N=\langle N,B_1/2,\ldots,B_5/2\rangle,
 \qquad q_M=2q_N.
\]

The Coble--Mukai lattice is computed independently as the saturated Smith
kernel

\[
 \operatorname{CM}=\{x\in M:x.h_i=0\text{ for }i=1,\ldots,5\},
 \qquad q_{CM}=q_M/2.                                    \tag{4.1}
\]

It is even unimodular of signature `(1,9)`, hence `E10=U+E8(-1)`.  The raw
divisors fold by

\[
 R\longmapsto R,
 \qquad
 E_{ij}\longmapsto \alpha_{ij}=2E_{ij}+h_i+h_j.          \tag{4.2}
\]

The ten images of each kind are integral `(-2)` roots and span the whole
Smith kernel, not a proper root sublattice.  Their labelled Gram is exactly
Kondō's `Gamma_VI`: the actual `(-2)` roots form the Petersen half, the
`(-1)` roots form its triangular-graph half, and the ten matched pairs have
intersection two.  The exact labelled isometry set has size `120`, recovering
the `S5` symmetry.

This is the classification-template gate, but it is reached only after the
boundary frame, Picard lattice, divisor marking, and integral CM saturation
have already been reconstructed from `P`.

## 5. All cusps, not only the roots

Exact rational double description of the twenty-root CM chamber gives

```text
657 rays = 600 finite + 57 primitive ideal + 0 hyperideal.
```

For every primitive cusp `f` and every recovered boundary half `h_i`, the
recognizer tests

\[
                        E=f-h_i\in N.                    \tag{5.1}
\]

All `57*5=285` pairs are tested.  The result is

```text
35 integral lifts: exactly 7 for each boundary;
22 cusps: no integral lift at any boundary;
0 cusps: lifts at two or more boundaries.
```

Every integral lift has square `-1`, boundary profile `(2,0,0,0,0)`, and is
one of the raw walls.  Conversely every one of the 35 raw walls of that
profile is obtained once.  This is the crucial full-Picard step: recognizing
only `Gamma_VI` would recover a codimension-five CM section, not the full nef
cone.

Finally the recognizer rebuilds the mirror system as

\[
 \{h_i\}_{5}
 \sqcup \{R\}_{10}
 \sqcup \{(\alpha_{ij}-h_i-h_j)/2\}_{10}
 \sqcup \{f-h_i\text{ integral}\}_{35}.                 \tag{5.2}
\]

This set is literally equal, with orientation, to the sixty input vectors.
It is therefore the same closed chamber, not merely an isomorphic nerve or a
CM face.

An independent full-cone replay gives

```text
28,872 rays = 27,845 finite + 1,027 ideal + 0 hyperideal;
all 60 facet faces have rank 14.
```

Thus the returned chamber is irredundant, finite-volume, and noncompact.
No inference of finite volume is made merely from finite sidedness.

## 6. Necessary-and-sufficient Type-VI recognizer

Let `R_VI` be the surface-free class of oriented primitive wall systems on an
even integral lattice of signature `(1,14)` for which the preliminary raw
arithmetic operations above are defined: primitive normalization,
reflection integrality, the Coxeter-angle gate, and finite rational cone
enumeration.  Do not assume a surface or a boundary marking.

For `(M,{r_a}) in R_VI`, define `VI-RETURN(M,{r_a})` by the finite tests:

1. the root-type and main-basis gates of Section 1;
2. exhaustive five-frame de-scaling and Coble conditions of Section 2;
3. the adjunction partition of Section 3;
4. saturated `CM=E10` and labelled `Gamma_VI` of Section 4;
5. exhaustive integral cusp lifting and literal same-wall return of Section
   5.

Then, over an algebraically closed field of characteristic different from
two,

\[
\boxed{
 \begin{array}{c}
 (M,\{r_a\})\text{ is the full-nef arithmetic chamber of a finite-Aut}\cr
 \text{strong Coble surface with five boundary components}
 \end{array}
 \iff
 \operatorname{VI\!\!-RETURN}(M,\{r_a\})\text{ accepts}.}
                                                               \tag{6.1}
\]

Necessity follows from adjunction, the definition of the extended Picard and
CM lattices, Kondō's description of the twenty effective irreducible roots,
and the Halphen/cusp dictionary.  Sufficiency uses the final labelled package
gate: Kondō proves existence and uniqueness of precisely this Type-VI
five-boundary surface; the same-wall return transports its complete nef cone
to the raw input.

Only at this last implication does algebraic classification enter.  It
forces

```text
characteristic 3;
one surface isomorphism class;
Aut(S)=S5;
R-invariant (E6 + D5, Z/2Z).
```

These are consequences extracted from the recovered polytope package, not
metadata supplied to the inverse search.  The primary source is Shigeyuki
Kondō, [*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591),
especially the Main Theorem, the Type-VI example, the root-exhaustion and
Halphen results, and the uniqueness section.

## 7. What extends to all nine classes

The same finite compiler has a uniform mathematical form.  Given a raw
oriented arithmetic wall datum, it must:

1. enumerate **all** possible boundary frames and all admissible de-scalings
   to an odd unimodular rational-surface Picard lattice;
2. reconstruct `K=-sum(B_i)/2`, verify characteristic parity, saturation,
   and adjunction;
3. fold actual `(-2)` and crossing `(-1)` walls into the saturated rank-ten
   CM lattice;
4. enumerate every primitive CM cusp and every integral boundary lift;
5. compare the resulting labelled root/glue package with the finite Kondō
   catalogue; and
6. demand exact same-oriented-wall return and a literal finite-volume gate.

On the independently defined raw screening class `C_Cob^raw` consisting of
systems that pass steps 1--4, completion of the nine finite package archives
would give the non-tautological decision theorem

```text
P has a finite-Aut strong Coble full-nef carrier
  iff
one of the nine classified package gates accepts and returns exactly P.
```

The domain is not defined as “the nine accepted answers”: it is the larger
surface-free C0--C2 numerical screening class.  The catalogue is a final
decision procedure on that domain.  Forward exact templates are presently
available in `H10`, `H11`, this `H14` row, and at the sharp `H17` endpoint;
see [`coble_nine_class_polytope_surface_transfer.md`](coble_nine_class_polytope_surface_transfer.md)
and [`coble_type_mii_full_ns_h17.md`](coble_type_mii_full_ns_h17.md).

Kondō's classification supplies an immediate **necessary carrier theorem**.
For every raw polytope that really has such a carrier,

\[
 b\in\{1,2,5,8\},\qquad
 \dim P=9+b\in\{10,11,14,17\}.                            \tag{7.1}
\]

Thus this carrier species has no members in `H12`, `H13`, `H15`, `H16`, or
above `H17`, while exact positive chambers occur at the `H14` and `H17`
boundaries.  This conclusion concerns Coble-carriability; it is not abstract
Coxeter-polytope nonexistence.

## 8. The still-missing L4 theorem

There is a crucial distinction between the finite nine-package **decision
theorem** above and the stronger transfer theorem the project ultimately
wants.

Numerical lattice, root, and cusp conditions do not by themselves prove that
the putative divisor classes are simultaneously effective on an algebraic
surface.  Consequently it would be circular to define `C_Cob` as “passes one
of the nine templates” and then advertise “there are nine types” as a new
classification of an independently existing polytope class.

The missing statement is an essential-surjectivity lemma of the following
kind.

> **Numerical Coble realization lemma (open).** Construct from every
> template-free raw package passing C0--C2 a finite algebraic realization
> stack, and prove that the raw conditions force this stack to be nonempty;
> its universal surface must have the reconstructed boundaries and roots,
> and the cusp compiler must return its entire nef chamber.

A concrete route is a *blowdown-incidence compiler*.  From `P`, enumerate all
finite contraction trees of recovered `(-1)` classes to `P2` or a Hirzebruch
surface.  For each tree construct a finite-type scheme whose variables are
the point cluster and curve equations, and whose open/closed conditions
enforce

```text
the divisor classes and incidences from P;
|-K| empty;
|-2K| the unique reconstructed boundary;
the required root curves effective and irreducible;
the prohibited incidences absent.
```

Three things would then have to be proved:

1. every genuine finite-Aut Coble carrier has a contraction tree enumerated
   by the compiler;
2. a template-free P-only criterion decides nonemptiness of at least a
   non-singleton family of these schemes; and
3. Kondō--Vinberg root exhaustion plus the cusp-lift calculation excludes
   hidden walls, giving exact same-P return.

The K3 double-cover period stack is an alternative in characteristic zero,
but the Type-VI and Type-MII endpoints require a positive-characteristic
crystalline realization theorem.  Until one of these bridges is proved, the
Coble program is accurately rated as:

```text
Type-VI inverse calibration:                    proved;
nine-package finite carrier decision theorem:  architecture, templates incomplete;
template-free automatic-carrier class (L4):    open essential-surjectivity lemma.
```

## 9. Reproduction

Run from the repository root:

```bash
PYTHONPATH=src:. python3 scripts/certify_polytope_first_coble_type_vi_h14.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_polytope_first_coble_type_vi_h14.py
```

The focused test checks the frozen exhaustive frame census, odd-unimodular
de-scaling, labelled CM saturation, all 285 cusp lifts, exact same-P return,
and the explicit scope boundary.
