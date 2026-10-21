# Compact automorphic pushback: the trace-reflection barrier and the exact L4 gap

**Frozen status (2026-08-09): no strict L4 theorem.**  Starting with an
arbitrary compact arithmetic Coxeter polytope, the canonical determinant line
does give a genuine projective Hilbert-orthogonal carrier and an exact marked
real-slice return.  It does **not** force a pure reflective automorphic divisor,
and the classical reflective-form theorems of Ma, Wang, and
Gritsenko--Hulek do not classify the resulting compact carrier.

The hostile audit gives a sharper reason than merely saying “the field is
wrong.”  If the ground field has degree `d>1`, an `F`-reflection has a
`d`-dimensional `(-1)`-eigenspace after restriction of scalars.  It is not a
rational pseudoreflection, and its fixed locus in the rational trace type-IV
domain has codimension `d`, not one.  The native Hilbert branch divisor
therefore has no **automatic** interpretation as a branch divisor for the full
rational trace orthogonal group.  A coincident rational-root hyperplane can
occur only as an additional root-by-root arithmetic condition; it is not the
same reflection and supplies no universal transfer.

Thus there is no noncircular, surface-free, non-singleton compact class in the
current determinant-line package on which existing classical AG results yield
a new dimension, finiteness, or nonexistence theorem for every original
polytope.  The exact completion package consists of two independent new
theorems: automatic bounded **pure** support on a raw Coxeter class, and a
Hilbert-reflective finiteness theorem.  Compensation is explicitly excluded.

## 1. Acceptance test

The requested implication has the form

```text
arbitrary P in a surface-free compact arithmetic class C
  -> every canonical determinant-line gluing
  -> a pure globally reflective divisor with exact same-P return
  -> independent AG classification or boundedness
  -> a theorem about every original P in C.
```

The following do not pass:

1. defining `C` by existence of the desired modular form, vanishing of its
   obstruction, bounded slope, or membership in a finite reflective list;
2. adding effective special divisors which miss the marked real slice;
3. recovering only the facet mirrors while allowing an untested mirror to cut
   the interior of `P`;
4. proving finiteness of type-IV carrier lattices without proving that `P` is
   the complete root chamber and exhausting determinant decompositions and
   gluings; or
5. applying a theorem for a finite-index subgroup of a rational orthogonal
   group to an infinite-index restriction-of-scalars subgroup.

The canonical determinant line is a fixed carrier species.  Every even
reflection-stable gluing is finite and enumerable, but all of them have the
same rational `F`-quadratic space, so the obstruction below applies to every
gluing simultaneously.

## 2. Setup

Let `P` be a compact arithmetic Coxeter polytope in `H^n`, with minimal totally
real ground field `F` of degree `d`.  After an even integral normalization let

\[
 L_P/\mathcal O_F
\]

have rank `n+1`, signature `(1,n)` at one real place `sigma_0`, and be negative
definite at every other real place.  Its determinant-line extension

\[
 \widehat L_P=L_P\perp K_P
\]

has rank `n+2` over `F` and signatures

\[
 (2,n),\qquad (0,n+2),\ldots,(0,n+2).
 \tag{2.1}
\]

Bruinier's construction associates to (2.1) a type-IV Shimura variety of
complex dimension `n`.  It is projective whenever the quadratic space is
anisotropic, and in this signature pattern anisotropy is automatic for
`d>1`: an `F`-isotropic vector would remain isotropic at a definite real
place.  This is the valid native algebraic carrier.

A facet root `r` defines an `F`-reflection `s_r`.  On the native type-IV
domain its fixed locus is a divisor, and its intersection with the canonical
real slice is the original hyperbolic facet mirror.  This establishes the
arrangement carrier; it says nothing yet about principality of the sum of
those divisors.

## 3. Compactness forces the non-rational branch

### Proposition 3.1

If `P` is compact and `n>=4`, then `d>1`.

### Proof

If `F=Q`, the Lorentzian quadratic space of `P` has dimension `n+1>=5`.
Meyer's theorem makes every indefinite rational quadratic form in at least
five variables isotropic.  The associated rational arithmetic hyperbolic
group is then nonuniform.  An arithmetic reflection subgroup of finite index
cannot have a compact fundamental polytope.  Hence a compact arithmetic `P`
in this range is not rational.  `square`

Compact dimensions two and three are deliberately excluded from this
proposition: anisotropic rational ternary and quaternary forms exist there.
The rational noncompact branch is treated separately in Section 7.

## 4. The trace-reflection barrier

Let `m=n+2`, let

\[
 T=\operatorname{Res}_{F/\mathbf Q}\widehat V_P
\]

be the underlying rational vector space, and use the trace bilinear form.
After an integral scaling it contains the trace lattice of `widehat L_P`.

### Theorem 4.1 (trace-domain and eigenspace barrier)

For `d=[F:Q]`, the trace quadratic space has

\[
 \operatorname{sig}(T)=\bigl(2,\,d(n+2)-2\bigr).
 \tag{4.1}
\]

Put `N=d(n+2)-2`.  Then:

1. the full rational type-IV domain `D_T` has complex dimension `N`;
2. the native Hilbert type-IV domain has dimension `n` and embeds in `D_T`
   with codimension

   \[
   N-n=(d-1)(n+2);
   \tag{4.2}
   \]

3. an `F`-reflection `s_r` has a `d`-dimensional `(-1)`-eigenspace over
   `Q`, namely the rational space underlying `F r`; hence it is a rational
   pseudoreflection if and only if `d=1`;
4. `F r` is negative definite for the trace form, so the fixed type-IV
   subdomain of `s_r` in `D_T` has codimension `d`; and
5. the image of a native facet divisor has codimension at least

   \[
   (d-1)(n+2)+1
   \tag{4.3}
   \]

   in the rational trace quotient.

### Proof

At the distinguished real place there are two positive and `n` negative
directions.  Each of the other `d-1` places contributes `n+2` negative
directions.  This proves (4.1), and (4.2) follows by subtracting dimensions.

Over `F`, the `(-1)`-eigenspace of `s_r` is the line `F r`.  Forgetting `F`
makes it a rational vector space of dimension `d`.  At every real embedding
the vector `r` is negative, so the trace form on that eigenspace is negative
definite.  Its orthogonal complement therefore has signature `(2,N-d)` and
its type-IV domain has dimension `N-d`.  This proves the fixed-locus claim.
Finally a native divisor has dimension `n-1`, while `D_T` has dimension `N`,
which gives (4.3).  `square`

The theorem does **not** say that the same set-theoretic hyperplane can never
be the pullback of a rational root divisor.  A particular trace vector may
coincidentally satisfy the rational integrality condition for another
rank-one reflection.  That is a separate exact root test.  The theorem says
that the original arithmetic reflection and its branch divisor do not pass
functorially to the classical rational branch.

### Proposition 4.2 (infinite-index barrier)

The natural rational algebraic subgroup

\[
 G=\operatorname{Res}_{F/\mathbf Q}SO(\widehat V_P)
   \longrightarrow SO(T)
\]

is proper for `d>1`.  More precisely,

\[
 \dim SO(T)-\dim G
   ={d(d-1)(n+2)^2\over2}>0.
 \tag{4.4}
\]

Consequently the native arithmetic group has infinite index in every full
trace orthogonal arithmetic group containing it.

### Proof

The two dimensions are

\[
 {d(n+2)(n+1)\over2},\qquad
 {d(n+2)(d(n+2)-1)\over2},
\]

whose difference is (4.4).  An arithmetic subgroup of `G` is Zariski dense
in `G`, while a finite-index subgroup of a full trace orthogonal arithmetic
group is Zariski dense in `SO(T)`.  Properness therefore rules out finite
index.  `square`

This also kills the obvious norm construction: multiplying a native Hilbert
form over cosets of the full trace group would require an infinite product.

## 5. Why the published classical theorems do not close the loop

Ma starts with an integral `Z`-lattice of signature `(2,N)` and the full
rational orthogonal modular variety.  His bounded-slope theorem gives
finiteness of such lattices carrying reflective forms of bounded vanishing
order.  It does not assert that a compact Hilbert-orthogonal branch divisor
extends to that rational quotient.

Wang's strongest current statements are likewise classical:

- there is no even `Z`-lattice of signature `(l,2)` with a reflective modular
  form for `l=21` or `l>=23`, except the dimension `26` case;
- the more structural root-system and finiteness theorems use hypotheses such
  as `2U+L` or `U+K`; and
- the complete `2`-reflective classification contains forty-two classical
  lattices for dimension at least eight.

None of those hypotheses follows from (2.1).  The determinant extension is
anisotropic over `F` and need not split even one hyperbolic plane over `F`.
After tracing to `Q`, the native group is infinite index and the facet
reflection is not the rational pseudoreflection in Theorem 4.1.

There are two valid one-way operations, neither of which gives the desired
pushback:

1. a classical form on the full trace quotient may be restricted to the
   Hilbert subdomain, but its other rational-root divisors can produce hidden
   slice mirrors and exact same-`P` return must be checked independently; and
2. nonexistence of a form for the full trace group does not imply nonexistence
   for its infinite-index Hilbert subgroup.

For illustration only, if one illegally fed a quadratic-field compact input
to Wang through the trace dimension, then `N=2n+2`: `n=11` would give the
excluded dimension `24`, `n=12` the exceptional dimension `26`, and all
`n>=13` would again be excluded.  This attractive pseudo-bound is **not a
theorem**, precisely because the trace-reflection and infinite-index gates
fail.

## 6. Exact support on the native Hilbert carrier

The native carrier still has a rigorous finite Borcherds gate.  For a fixed
finite principal part, Bruinier's cusp pairing is necessary and sufficient
for a weakly holomorphic Whittaker input.  A deep congruence level can
separate the selected facet sheets, and root integrality of every selected
global component can be checked exactly.

This gives **decidability for a fixed Borcherds support**, not automatic pure
support for every compact reflective chamber.  The Bugaenko `H8/11` control
is already decisive at that logical level:

- its determinant-line gluing is unique;
- the pure norm-two Whittaker column is obstructed; and
- the complete first slice-disjoint compensator shell is also obstructed.

The general safe-column theorem forces a correction at some larger support,
but that correction is precisely what the strict L4 test forbids.  It may
miss the real slice and still make the global divisor non-reflective.

Nor does the Whittaker obstruction rule out every arbitrary automorphic form
with the same special support.  Existing converse theorems for Borcherds
products require a hyperbolic split; Metzler's 2025 theorem explicitly keeps
that hypothesis.  The compact anisotropic Hilbert lattice has no such split.
Thus “arbitrary pure automorphic support exists” is not presently reduced to
the finite Whittaker matrix in this branch.

Unbounded compensation cannot be substituted.  On any projective carrier,
for an ample class `A` and effective target `T`, Serre vanishing gives

\[
 T+C_k\sim kA,\qquad C_k\geq0
\]

for all sufficiently large `k`.  Therefore existence of an effective
automorphic divisor containing the facets, with uncontrolled extra support
or slope, is algebro-geometrically tautological and cannot bound `P`.

## 7. Compact and noncompact branches

The obstruction above is specific to the compact `n>=4` branch.

For `F=Q`, the determinant carrier is a genuine classical integral
orthogonal variety.  If one **assumes** a globally reflective form with
bounded multiplicity, Ma's finiteness theorem and Wang's classifications can
be applied to the carrier lattice.  But:

1. rational arithmetic `P` is noncompact in `n>=4` by Proposition 3.1;
2. pure support is still an additional condition, as the known rational
   determinant-line positive and negative controls show;
3. the original `P` must be proved to be the complete root chamber; and
4. determinant decompositions, gluing alternatives, hidden mirrors, and
   cusp finite-volume tests must still be exhausted.

Consequently even a finite list of classical carrier lattices is not, by
itself, a finite list of the original polytopes.  No compact/noncompact
inference is mixed in this audit.

## 8. Hostile audit of candidate classes

| proposed class | raw and surface-free? | purity forced? | independent compact AG theorem? | verdict |
|---|---:|---:|---:|---|
| all compact marked arithmetic `P`, `n>=4` | yes | no | no | fail |
| compact complete root chambers | yes | no | no | fail |
| vanishing Hilbert Whittaker obstruction | lattice-computable, but modular-carrier-defined rather than independently Coxeter-defined | Borcherds-only | no | fails independence and AG gates |
| bounded-slope or reflective-form class | no | assumed | Ma only in classical branch | circular |
| trace lattice belongs to a Ma/Wang list | list membership is raw but encodes the terminal classification | no functorial facet reflection; exact restriction unproved | classical only | fail |
| compensated effective divisor | yes | explicitly no | ample completion is tautological | fail |

No row passes all four gates.  In particular, there is no presently known
non-singleton compact high-dimensional subclass meeting the user's strict
acceptance test.

## 9. Kodaira consequences are outputs, not contradictions

Gritsenko--Hulek show that a strongly reflective modular form of sufficiently
large weight can force a classical orthogonal modular variety to be
uniruled; related formulas determine zero or negative Kodaira dimension in
special weights.  These results begin **after** the pure strongly reflective
form is supplied.

They do not force that form from `P`, bound its weight, or contradict the
existence of a compact Coxeter slice.  A neat compact type-IV cover has ample
automorphic canonical class and is already of general type, while a branched
finite quotient can have lower Kodaira dimension.  General type of the
carrier also places no bound on the dimension of its totally real Coxeter
slice.  Bruinier constructs compact Hilbert-orthogonal varieties in
arbitrarily large dimensions divisible by four over suitable real quadratic
fields, so projectivity or generic Kodaira type alone cannot be the missing
restriction.

## 10. The minimum completion package

There is not one missing technical estimate.  Two logically independent
theorems are necessary.

### Lemma A: automatic bounded purity

Find an independently defined, non-singleton compact raw Coxeter class `C`
such that, for every `P in C` and every canonical determinant-line gluing,
either:

- a pure globally reflective divisor exists, with multiplicity or normalized
  slope bounded independently of the number of facets, and its zero-auxiliary
  mirrors return exactly the complete chamber `P`; or
- that gluing is excluded by an exact finite obstruction.

This lemma must include a compact Hilbert converse if arbitrary automorphic
forms, rather than only Borcherds products, are allowed.

### Lemma B: Hilbert-reflective finiteness

For the anisotropic varieties

\[
 \operatorname{Res}_{F/\mathbf Q}SO(2,n),
\]

prove that pure reflective forms of the bounded multiplicity/slope in Lemma A
occur in only finitely many effectively enumerable `F`-lattices and
dimensions.  The theorem must be native over `F`; Theorem 4.1 proves that a
restriction-of-scalars appeal to the classical list is not automatic.

Only `A+B`, followed by complete determinant-decomposition and Vinberg return
tests, would produce the requested theorem on every original `P`.  Positive
terminal rows could then produce examples as well.

The frozen verdict is

```text
C        = PASS for the broad raw compact arithmetic input class
Carrier  = PASS for the projective determinant-line Hilbert arrangement
Purity   = FAIL: fixed Borcherds support is decidable but not forced;
           arbitrary support lacks a compact Hilbert converse
AG       = FAIL: no native Hilbert-reflective finiteness/classification
Trace    = FAIL as a universal bridge: d-dimensional (-1)-eigenspace,
           codimension-d fixed locus, infinite-index subgroup
Return   = conditional exact marked-slice return only
Pushback = FAIL; no new polytope bound or nonexistence theorem claimed
Rating   = L3 carrier + theorem-shaped L4 scope obstruction, not L4
```

## 11. Reproducibility

The exact symbolic certificate is
[`automorphic_pushback_l4_hostile_audit.json`](../data/enumerations/automorphic_pushback_l4_hostile_audit.json),
generated by
[`certify_automorphic_pushback_l4_hostile_audit.py`](../scripts/certify_automorphic_pushback_l4_hostile_audit.py)
and checked by
[`test_automorphic_pushback_l4_hostile_audit.py`](../tests/test_automorphic_pushback_l4_hostile_audit.py).

```bash
PYTHONPATH=src:. python3 scripts/certify_automorphic_pushback_l4_hostile_audit.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_automorphic_pushback_l4_hostile_audit.py
```

## 12. Primary sources

1. J. H. Bruinier,
   [*Regularized theta lifts for orthogonal groups over totally real fields*](https://arxiv.org/abs/0908.3076):
   the native Hilbert-orthogonal Shimura variety, projectivity for the
   anisotropic `d>1` branch, and Whittaker/Borcherds special divisors.
2. S. Ma,
   [*On the Kodaira dimension of orthogonal modular varieties*](https://arxiv.org/abs/1701.03225):
   rational `Z`-lattice scope and bounded-slope/vanishing-order finiteness.
3. H. Wang,
   [*2-reflective lattices of signature `(n,2)` with `n>=8`*](https://arxiv.org/abs/2301.11536):
   the forty-two-lattice classical classification.
4. H. Wang,
   [*On the classification of reflective modular forms*](https://arxiv.org/abs/2301.12606):
   the dimensions `21` and `>=23` exclusion, the dimension `26` exception,
   and the precise `2U+L` / `U+K` finiteness hypotheses.
5. V. Gritsenko and K. Hulek,
   [*Uniruledness of orthogonal modular varieties*](https://arxiv.org/abs/1202.3361):
   the Kodaira/uniruledness consequence after a strongly reflective form is
   already present.
6. I. Metzler,
   [*A new converse theorem for Borcherds products*](https://arxiv.org/abs/2512.20701):
   the current converse theorem and its explicit hyperbolic-split hypothesis.

Related project controls are the
[`strict compact inverse leverage audit`](borcherds_compact_inverse_leverage_audit.md),
the
[`determinant-line transfer theorem`](polytope_to_automorphic_divisor_transfer_theorem.md),
and the
[`Bugaenko determinant-line calibration`](bugaenko_h8_determinant_line_ag.md).
