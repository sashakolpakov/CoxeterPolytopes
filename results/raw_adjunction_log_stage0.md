# Raw adjunction/log Stage 0

**Status (2026-08-10): proved and executable on arbitrary raw input.**  The
compiler starts with only an integral unimodular Lorentzian lattice Gram
matrix and primitive oriented full-span negative wall vectors.  It accepts
no canonical class, physical scale, boundary marking, surface label,
characteristic, volume assertion, or carrier verdict.  It exhausts every
smooth-rational numerical scaling and every canonical vector allowed by
adjunction, characteristic parity, and the rational-surface square identity;
then it detects exactly the direct weak branch or a disjoint standard log
branch whose wall reflections preserve the inferred arithmetic module.

This is a genuinely raw **numerical** Stage 0, not a surface-existence
theorem.  It does not infer effectivity, irreducibility, complete Mori
generation, nonempty moduli, or same-polytope return.  No
Nikulin--Vinberg, Kond\=o, ABSW, Prasad, arithmetic-volume, or surface
classification theorem is used.

## 1. Closed raw interface

The input is

\[
       (L,G;r_1,\ldots,r_N),
\]

where `G` is an integral unimodular Gram matrix of inertia
\((1,\rho-1,0)\), and the `r_i` are distinct primitive oriented integral
vectors with negative square.  They must span \(L_{\mathbf Q}\), and their
pairings must be nonnegative off the diagonal.  The last condition is the
necessary intersection condition for distinct effective irreducible curves;
it is not treated as proof that such curves exist.

The machine-readable schema has exactly three keys:

```text
schema, picard_gram, wall_vectors.
```

Additional keys are rejected.  In particular a purported `canonical_class`
or boundary list cannot leak into the computation.  Full span is recomputed,
not supplied as a trusted flag.  No finite-volume flag is needed.

The exact implementation is
[`certify_raw_adjunction_log_stage0.py`](../scripts/certify_raw_adjunction_log_stage0.py),
the closed schema is
[`raw_adjunction_log_stage0_input.schema.json`](../data/fixtures/raw_adjunction_log_stage0_input.schema.json),
and the six frozen raw arrays are
[`raw_adjunction_log_stage0_inputs.json`](../data/fixtures/raw_adjunction_log_stage0_inputs.json).

## 2. Finite scaling theorem

Suppose a wall ray is the ray of a smooth rational curve and write its class
as

\[
                       C=d r,\qquad d\in\mathbf Z_{>0},
\]

where \(r\) is the supplied primitive vector and \(r^2=-q<0\).  If \(K\)
is the canonical class, adjunction says

\[
 -2=C^2+K.C=-qd^2+d(K.r).
\]

Consequently

\[
                         K.r=qd-\frac2d.                 \tag{2.1}
\]

Since the left side is integral, \(d\mid2\).  Moreover \(K\) is
characteristic, so

\[
                        K.r\equiv r^2=-q\pmod2.           \tag{2.2}
\]

The complete local table is therefore

| primitive square | scale | forced pairing | parity condition |
|---:|---:|---:|---|
| \(-q\) | `1` | \(q-2\) | always |
| \(-q\) | `2` | \(2q-1\) | exactly when \(q\) is odd |

There are no other physical scales compatible with smooth-rational
adjunction and characteristic parity.

Choose any \(\rho\) independent wall vectors.  Their forced pairings
determine \(K\) uniquely over \(\mathbf Q\).  Each basis wall has at most two
pairing choices, hence at most

\[
                              2^\rho                       \tag{2.3}
\]

vectors need be tested.  The implementation uses the matroid greedy
algorithm with weight zero on even-\(q\) walls and weight one on odd-\(q\)
walls, so it chooses a basis having the smallest possible number of binary
choices.  Every candidate is then tested against **all** walls, required to
be integral and globally characteristic, and required to satisfy

\[
                              K^2=10-\rho.                  \tag{2.4}
\]

Equation (2.4) is the Noether identity for a rational surface.  It is only a
necessary numerical gate here; the raw lattice does not prove rational
surface existence.

This proves termination and exhaustiveness.  Full span is essential: without
it, wall pairings leave an affine family of possible canonical vectors and
the stated finite enumeration no longer follows.

## 3. The two reflection modules must not be conflated

Strict reflection integrality on \(L=\operatorname{NS}(S)\) cannot be
imposed on every raw Coble wall.  For a primitive vector in a unimodular
lattice, the reflection in \(r\) preserves \(L\) only if \(-r^2\mid2\).
Thus it covers the weak roots and the twice-primitive Type-II boundary, but
it deliberately does **not** cover a primitive physical \((-4)\)-curve.

The compiler consequently tests two modules, inferred from the scaling
rather than supplied in the input.

### 3.1 Direct weak module

If every physical square is `-1` or `-2`, the candidate is accepted exactly
when every wall reflection preserves \(L\).  In fact the local table then
forces

\[
 (q,d)=(1,1)\quad\hbox{or}\quad(2,1),
\]

so this is precisely the weak `(-1)/(-2)` alphabet.  If these walls are
later proved to generate the complete Mori cone, \(-K\) is nef and
\(\rho\leq10\), giving full-nef chamber dimension at most `9`.

### 3.2 Inferred standard boundary module

Every physical wall outside the weak alphabet must have even square

\[
                          C_i^2=-2m_i,\qquad m_i\geq2.       \tag{3.1}
\]

The assignment itself determines the boundary set and the integers \(m_i\).
The compiler forms

\[
 M_\Delta=L+\sum_i\mathbf Z\frac{C_i}{m_i},\qquad
 \Delta=\sum_i\left(1-\frac1{m_i}\right)C_i.               \tag{3.2}
\]

It then verifies, by exact Hermite-normal-form membership calculations,
that every wall reflection preserves \(M_\Delta\).  It also verifies that the
form multiplied by \(\ell=\operatorname{lcm}(m_i)\) is integral.  The
standard log branch is accepted exactly when

\[
 \begin{gathered}
 C_i.C_j=0\quad(i\ne j),\\
 K+\Delta=0\quad\hbox{in }L\otimes\mathbf Q,\\
 s_{r_j}(M_\Delta)=M_\Delta\quad\hbox{for every wall }j.
 \end{gathered}                                               \tag{3.3}

No boundary indices or weights enter the input.

This formulation distinguishes the two ways a physical `-4` curve occurs:

- `C=2r`, with `r^2=-1`: the physical class is twice primitive and its
  reflection is already integral on `NS`; this is the `H10/30` Type-II row;
- `C=r`, with `r^2=-4`: the physical boundary is primitive and its reflection
  is integral only after adjoining `C/2`; this is what occurs in the
  `H11/26`, `H14/60`, and `H17/48` rows.

Treating these cases as the same integral-lattice statement would be a real
error.

The independent adapted-cover theorem in
[`rational_nef_reflection_integrality_gate.md`](rational_nef_reflection_integrality_gate.md)
shows that an actual tame rational-surface realization of a nonempty branch
(3.3) forces every \(m_i=2\).  The canonical adapted cover is then K3.  Thus
the only tame geometric output of the standard-coefficient branch is the
Coble index-two branch, with the conditional Picard ceilings `20` in
characteristic zero and `22` in odd positive characteristic.  Stage 0 itself
records generalized numerical branches before applying that geometric
theorem; it does not silently assume a characteristic.

## 4. Exact calibrations

The frozen fixtures are projections of existing physical lattice arrays,
but the solver sees only the Gram and wall vectors.  Source canonical
classes occur only in the separate expected-output field used after the
solve; source boundary markings are omitted entirely.  The raw key-set check
proves that neither enters the mathematical input.

| raw case | \(\rho\) | walls | basis assignments | numerical \(K\)'s | accepted output |
|---|---:|---:|---:|---:|---|
| toric `H6/9` | 7 | 9 | 2 | 2 | one weak branch, `6*(-2)+3*(-1)` |
| Coble Type II `H10/30` | 11 | 30 | 2 | 2 | one index-two log branch, `1*(-4)+12*(-2)+17*(-1)` |
| Coble Type I `H11/26` | 12 | 26 | 1 | 1 | one index-two log branch, `2*(-4)+10*(-2)+14*(-1)` |
| Coble Type VI `H14/60` | 15 | 60 | 1 | 1 | one index-two log branch, `5*(-4)+10*(-2)+45*(-1)` |
| Coble Type MII `H17/48` | 18 | 48 | 1 | 1 | one index-two log branch, `8*(-4)+24*(-2)+16*(-1)` |
| `Bl_p(F_4)` maximal triangle | 3 | 3 | 4 | 1 | none |

For every accepted Coble row, the boundary is inferred as the full set of
physical `-4` walls, it is pairwise disjoint, and the remaining incidence
table is recovered without labels:

\[
              (C^2,B.C)=(-2,0)\quad\hbox{or}\quad(-1,2).     \tag{4.1}

\]

The `H10` datum has a second integral characteristic vector of the required
square, but it scales seventeen primitive `-1` normals to physical `-4`
classes and fails \(K+\Delta=0\).  The toric datum similarly has a false
three-boundary alternative.  Retaining and explicitly rejecting these
solutions is important: merely finding *a* familiar canonical vector would
not prove that the raw inverse is unique.

### The maximal `Bl_p(F_4)` control

In the basis \((B,F,E)\), Stage 0 receives only

\[
 G=\begin{pmatrix}-4&1&0\\1&0&0\\0&0&-1\end{pmatrix},
 \qquad r=(B,E,F-E).
\]

It uniquely recovers

\[
                    K=(-2,-6,1),\qquad
                    C^2=(-4,-1,-1).
\]

It infers `B` as the sole index-two boundary candidate.  All three
reflections preserve the inferred half-boundary module, and the boundary is
vacuously disjoint.  Nevertheless

\[
                         K+\frac12B=(-3/2,-6,1)\ne0.          \tag{4.2}

\]

The case is therefore rejected exactly at the log equation after passing the
shape, disjointness, and adapted-module tests.  The fact that its Coxeter
triangle is reflection-maximal supplies no missing implication.  This is the
requested executable obstruction to deriving the Coble equation from
reflection maximality.

## 5. Invariance and verification

A unimodular basis change \(P\) sends

\[
             G\mapsto P^tGP,\qquad r_i\mapsto P^{-1}r_i,
             \qquad K\mapsto P^{-1}K.                         \tag{5.1}

\]

Every condition used by the compiler is preserved: primitivity, wall rank,
pairings, adjunction, characteristic parity, canonical square, physical
scales, disjointness, the log equation, and reflection preservation of the
transported adapted module.  A wall permutation merely permutes the local
equations and inferred boundary indices.  Thus the candidate set is
intrinsic even though the printed coordinates and the chosen computational
basis are not.

The certificate applies a deterministic nontrivial unimodular basis mutation
and a wall-order reversal to every calibration.  It transports every
recovered canonical vector back and checks equality of the complete set, not
just of the accepted branch.  All six cases pass.

Executable verification:

```bash
/Users/sasha/.venvs/global/bin/python \
  scripts/certify_raw_adjunction_log_stage0.py \
  --check-fixtures --check-archive

/Users/sasha/.venvs/global/bin/python -m pytest -q \
  tests/test_raw_adjunction_log_stage0.py
```

The frozen exact output is
[`raw_adjunction_log_stage0.json`](../data/enumerations/raw_adjunction_log_stage0.json).
The focused suite has eleven tests, including malformed-input rejection,
the two false canonical alternatives, the primitive/twice-primitive boundary
distinction, exact `Bl_p(F_4)` rejection, archive replay, and both invariance
mutations.

Pinned artifact hashes after the terminal replay:

```text
script   a60d4780c9897daf2c60d3bf432c3dc41983449c531ffada1871ba5c118171b5
schema   f007d407886b729259796d85891bbce49662e4631231e9444cebb48b19db97d9
fixtures c197624eb2986e44c2a991f3433cfe81d5b16ca00593282b481844721dcc4113
archive  55700e272eaeb527b0bcc913ab401a45898b9d076cf9de50f5159986f6403730
tests    d80c35143395b7474ca07ea421a6ad1103dc4fa010da0940c29fe532efc38348
```

## 6. What this advances, and what remains

The raw arrows now proved are

```text
(unimodular Lorentzian Gram, primitive oriented full-span walls)
    -> every adjunction/characteristic/Noether scaling and K
    -> exact weak-NS or inferred-standard-log branch gate
    -> conditional rank return once geometric completeness exists.
```

The remaining theorem is not another lattice census.  One must still prove,
for the intended polytope class, that the selected wall classes are effective
smooth rational curves on a nonempty surface and exhaust the complete Mori
cone (or prove a finite exhaustive alternative of carrier species).  Only
then does the weak `rho<=10` or Coble/K3 `rho<=20,22` return become a
dimension theorem for the original polytopes.  Stage 0 removes canonical and
boundary data as hidden inputs; it does not manufacture the missing
geometric existence and completeness arrow.
