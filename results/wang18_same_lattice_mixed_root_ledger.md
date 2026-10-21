# Wang's eighteen negative rows: exact reflection-maximality ledger

**Status (2026-08-10): closed, 18/18.**  There are no undecided rows.

This audit gives a concrete answer to the proposed shortcut
"maximal norm-two reflection group implies automorphic correction."  The
shortcut is false, but it misses by exactly one lattice in the classified
rank-at-least-six denominator:

> **Theorem.**  Let `S` be one of the 69 complete norm-two reflective even
> hyperbolic lattices of rank at least six.  If `W_2(S)` is maximal as an
> arithmetic reflection group, then either `S` has a Wang automorphic
> correction or
>
> \[
> S\cong U\oplus E_8(-1)\oplus E_7(-1).
> \]
>
> The exceptional group is Tumarkin's cofinite `Gamma_16`.  It has no Wang
> automorphic correction and is maximal **as a reflection group**.

The theorem is one-way.  It does not say that every corrected row is
reflection-maximal.  It also does not call `Gamma_16` a maximal arithmetic
lattice: its diagram symmetry gives a non-reflective normalizer extension.

The proof has three exact pieces.  Wang's Theorem 7.2 says that precisely 18
of the 69 rows have no correction.  The uniform construction below produces
a strict same-lattice reflection enlargement for 17 of those 18.  The
separate [`Gamma_16` hostile audit](wang_h16_maximal_reflection_no_correction.md)
proves maximality of the remaining row, including commensurable reflection
overgroups.  Thus the result does not use a determinant table as a surrogate
for arbitrary arithmetic maximality.

## 1. The three notions kept separate

For an even lattice `S` of signature `(1,n)`, put

\[
 W_2(S)=\langle s_a:a\in S,\ a^2=-2\rangle,
 \qquad
 W_{\rm ref}(S)=\langle s_r:s_r(S)=S,\ r^2<0\rangle.
\]

The ledger distinguishes:

1. **same-lattice enlargement:** `W_2(S) < W_ref(S)`;
2. **integral overclosure:** a proper integral `S < M < S tensor Q` (even
   or odd); and
3. **arithmetic-reflection maximality:** no larger arithmetic reflection
   group in the same `Isom(H^n)`, allowing commensurable lattice changes.

A proper overclosure is neither necessary nor, without new roots, sufficient
for a strict reflection enlargement.  Indeed, rows 12 and 17 have squarefree
determinant and no integral overclosure, yet each already has a longer root
inside `S`.  Conversely, row 1 requires the separate normalizer/maximality
argument to pass from integral lattice arithmetic to arbitrary commensurable
reflection groups.

## 2. Uniform mixed-root lemma

Write `S=U+R(-1)` and choose a basis `e,f` of `U` with
`e^2=f^2=0` and `(e,f)=1`.  Suppose that `omega in R^*` has exact order `m`
modulo `R` and positive-definite square

\[
                         \omega^2=\frac{2}{m}.
\]

Then

\[
                         r=m e+m\omega\in S
\]

is primitive, has `r^2=-2m` and `div_S(r)=m`, and its reflection is integral:

\[
 s_r(x)=x-\frac{2(x,r)}{r^2}r
       =x+\frac{(x,r)}{m}r.
\]

Primitivity follows from the exact order of `omega`; pairing with `f` gives
the exact divisibility.  For `m=2,3`, this reflection is not in `W_2(S)`.
Indeed, every codimension-one reflection in the Coxeter group `W_2(S)` is
conjugate to a norm-two wall reflection, hence has a norm-two rational normal
on the same line.  The line of `r` cannot have such a rational normal: its
norm ratio to `-2` is `m`, and neither 2 nor 3 is a rational square.

Since Wang's rows are 2-reflective, `W_2(S)` already has finite covolume.
Consequently every strict `W_ref(S)` obtained here is itself an arithmetic
finite-covolume reflection overgroup.  No new Vinberg run is needed.

The required discriminant classes are elementary:

- `D_n`: the vector weight `omega_1`, order 2, square 1;
- `A_3`: the middle weight `omega_2`, order 2, square 1;
- `2A_1`: `omega=(a+b)/2`, order 2, square 1; and
- `A_2`: a minuscule weight, order 3, square `2/3`.

Thus the first three give a `(-4,div 2)` root and the last gives a
`(-6,div 3)` root.

## 3. Complete 18-row table

Negative signs on ADE summands are suppressed.  `D4/2` in the root column
means that the displayed `D4` class supplies a root of square `-4` and
divisibility 2.  The overclosure column records one explicit glue class, not
an assertion that it is the unique or maximal overclosure.

| # | `S` | `H^n` | `|det S|` | same-lattice extra root | proper integral overclosure | maximal arithmetic reflection? |
|---:|---|---:|---:|---|---|---|
| 1 | `U+E8+E7` | 16 | 2 | none | none (squarefree) | **yes**, `Gamma16` |
| 2 | `U+E8+D6` | 15 | 4 | `D6: (-4,div 2)` | index 2, odd, D6 vector glue | no |
| 3 | `U+E8+D4+A1` | 14 | 8 | `D4: (-4,div 2)` | index 2, odd, D4 vector glue | no |
| 4 | `U+E8+D4` | 13 | 4 | `D4: (-4,div 2)` | index 2, odd, D4 vector glue | no |
| 5 | `U+D8+D4` | 13 | 16 | `D8: (-4,div 2)` | index 2, odd, D8 vector glue | no |
| 6 | `U+E8+4A1` | 13 | 16 | `2A1: (-4,div 2)` | index 2, odd, `(a+b)/2` | no |
| 7 | `U+E8+3A1` | 12 | 8 | `2A1: (-4,div 2)` | index 2, odd, `(a+b)/2` | no |
| 8 | `U+D8+3A1` | 12 | 32 | `D8: (-4,div 2)` | index 2, odd, D8 vector glue | no |
| 9 | `U+E8+A3` | 12 | 4 | `A3: (-4,div 2)` | index 2, odd, A3 middle glue | no |
| 10 | `U+D8+2A1` | 11 | 16 | `D8: (-4,div 2)` | index 2, odd, D8 vector glue | no |
| 11 | `U+2D4+2A1` | 11 | 64 | `D4: (-4,div 2)` | index 2, odd, D4 vector glue | no |
| 12 | `U+E8+A2` | 11 | 3 | `A2: (-6,div 3)` | none (squarefree) | no |
| 13 | `U+E6+A2` | 9 | 9 | `A2: (-6,div 3)` | index 3, even, diagonal `E6+A2` glue | no |
| 14 | `U+D4+A3` | 8 | 16 | `D4: (-4,div 2)` | index 2, odd, D4 vector glue | no |
| 15 | `U+D5+A2` | 8 | 12 | `D5: (-4,div 2)` | index 2, odd, D5 vector glue | no |
| 16 | `U+D4+A2` | 7 | 12 | `D4: (-4,div 2)` | index 2, odd, D4 vector glue | no |
| 17 | `U+A4+A2` | 7 | 15 | `A2: (-6,div 3)` | none (squarefree) | no |
| 18 | `U+A3+A2` | 6 | 12 | `A3: (-4,div 2)` | index 2, odd, A3 middle glue | no |

Totals:

\[
 \boxed{17\text{ strict same-lattice reflection enlargements},\quad
        15\text{ proper integral overclosures},\quad
        1\text{ maximal negative row},\quad
        0\text{ undecided}.}
\]

Every integral overclosure in the table is automatically `W_2(S)`-stable:
if `x in M` and `a in S` has square `-2`, then
`s_a(x)=x+(x,a)a in M`.  For the three negative answers, squarefreeness is
decisive because an index-`k` integral overclosure would satisfy
`|det S|=k^2|det M|`.

## 4. Why `U+E8+E7` is the unique same-lattice exception

For `S=U+E8(-1)+E7(-1)`,

\[
 A_S\cong \mathbf Z/2,\qquad q_S(1)=+\tfrac12\pmod {2\mathbf Z}.
\]

Let `r` be a primitive negative vector whose reflection preserves `S`, write
`r^2=-2d`, and put `delta=div_S(r)`.  Reflection integrality gives
`d | delta`, while the class `r/delta` has exact order `delta` in `A_S`.
If `d>1`, the order-two discriminant group forces `d=delta=2`.  But then

\[
 q_S(r/2)=r^2/4=-1\pmod {2\mathbf Z},
\]

contradicting `q_S(1)=1/2`.  Hence every integral reflection of `S` is a
norm-two reflection and `W_ref(S)=W_2(S)`.  Determinant 2 also excludes every
proper integral overclosure.

Those two facts alone do **not** exclude a reflection overgroup after a
nonintegral commensurable change of lattice.  That last step is exactly what
the separate [`Gamma16` hostile
audit](wang_h16_maximal_reflection_no_correction.md) supplies; it also keeps
reflection maximality distinct from the larger non-reflective diagram
normalizer.

The companion argument is short and Prasad-free.  A discrete reflection
overgroup contains `Gamma16` with finite index, remains nonuniform, and
Felikson--Tumarkin's facet-rank theorem forces its chamber to have at most 19
facets.  In `H16`, 17 facets would be a finite-volume Coxeter simplex (none
exists), the complete `n+2` classification excludes 18 facets, and Tumarkin
proves the noncompact 19-facet polytope unique.  Hence the two chambers are
congruent, the index is one, and `Gamma16` is maximal among **all discrete
reflection overgroups**, not merely integral or arithmetic ones.

## 5. Independent cross-check and reproduction

Row 12 was already present independently in the divisor-fold stock.  Its
frozen `U+E8(-1)+A2(-1)` full chamber has twelve `(-2,div 1)` walls and one
`(-6,div 3)` wall; Scharlau--Kirschmer List 11, entry 4 identifies the same
13-facet full reflection chamber.  See
[`roulleau_ue8d4_all_subgroup_folds.md`](roulleau_ue8d4_all_subgroup_folds.md).
This independently matches the uniform `A2` construction above.

The exact certifier constructs the ADE Cartan matrices, solves for the
fundamental weights over `Q`, and checks coordinate gcd, norm, divisibility,
integrality of every reflection coefficient, glue index, glue square,
parity, and determinant for all rows:

```bash
python3 scripts/certify_wang18_same_lattice_mixed_root_ledger.py
python3 -m pytest -q tests/test_wang18_same_lattice_mixed_root_ledger.py
```

The replay reports `18` rows, `17` mixed-root enlargements, `15` integral
overclosures, and `0` undecided; all six focused tests pass.

## 6. Primary lineage and claim boundary

- Haowu Wang, [*The classification of 2-reflective modular
  forms*](https://doi.org/10.4171/JEMS/1358), JEMS 26 (2024), Theorem 7.2:
  published `51+18` automorphic-correction split and the exact eighteen-row
  input list.
- Ruth Kellerhals, [*On Minimal Covolume Hyperbolic
  Lattices*](https://doi.org/10.3390/math5030043), Example 6: Tumarkin's
  arithmetic cofinite `Gamma_16`, its 19-generator diagram, and its exact
  Euler characteristic/covolume.
- Anna Felikson and Pavel Tumarkin, [*On finite index reflection subgroups of
  discrete reflection groups*](https://arxiv.org/abs/math/0305093): the
  facet-rank inequality for finite-index reflection subgroups.
- Pavel Tumarkin, [the complete finite-volume `n+2`-facet
  classification](https://arxiv.org/abs/math/0301133) and [the unique
  noncompact 19-facet polytope in
  `H16`](https://arxiv.org/abs/math/0311272): the low-facet maximality gates.
- Rudolf Scharlau and Markus Kirschmer, [reflective Lorentzian lattice
  tables](https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/):
  independent full-reflection data used for the row-12 cross-check and in the
  separate `Gamma_16` audit.

The 69 lattices, Wang's `51+18` split, and Tumarkin's polytope are published.
The uniform 17-row mixed-root elimination and the resulting unique-exception
maximality criterion are new project-level deductions.  No claim of prior
publication or worldwide novelty is made without a literature audit.
