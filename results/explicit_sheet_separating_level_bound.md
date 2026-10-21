# An explicit per-polytope sheet-separating level

**Status (2026-08-09): proved.**  The deep level in the exact facet-support
theorem is computable from the marked compact polytope.  This removes an
existential ambiguity, but it does not give a uniform bound on that level for
an infinite polytope class.

## 1. Setup

Let `M` be the unglued determinant-line lattice, viewed through its underlying
integral trace lattice, and fix a `Z`-basis `e_1,...,e_R`.  Let
`r_1,...,r_s` be the marked facet vectors.  For the finite set of target
norms

\[
 T=\{Q(r_i):1\leq i\leq s\},
\]

let

\[
 \mathcal U_T=
 \{v\in M:Q(v)\in T,\;v^\perp\cap P\ne\varnothing\}/\{\pm1\}.
\tag{1.1}
\]

Compactness of `P`, definiteness at every nondistinguished real place, and
discreteness of the Minkowski embedding make this set finite.  The standard
positive-majorant argument also makes it effectively enumerable from the
marked Gram data and an exact compactness certificate.

Form the finite set `Delta(P)` of all nonzero differences which must remain
nonzero modulo the new level:

\[
 \begin{split}
 \Delta(P)=&\{v-r_i,\;v+r_i:
       v\in\mathcal U_T,\;Q(v)=Q(r_i),\;v^\perp\ne r_i^\perp\}\\
 &\cup\{r_i-r_j,r_i+r_j: i<j,\;r_i^\perp\ne r_j^\perp\}.
 \end{split}
\tag{1.2}
\]

Write `d=sum_k d_k e_k` for `d in Delta(P)` and define

\[
 N_{\rm sep}(P)=
 1+\max\{|d_k|:d\in\Delta(P),\;1\leq k\leq R\},
\tag{1.3}
\]

with `N_sep(P)=1` if `Delta(P)` is empty.

## 2. Explicit separation theorem

### Theorem 2.1

For `M_P=N_sep(P) M`, every target class

\[
 \mu_i=r_i+M_P\in M_P'/M_P
\]

is defined, distinct modulo sign from the other target classes, and its
same-norm special divisor meets the closed real chamber exactly in the facet
`r_i^perp cap P`.  It contains no other facet, lower-face contact, or
interior-cutting mirror.

#### Proof

Because `r_i in M` and `M` is integral,
`(r_i,N_sep(P)M)` is integral, so `r_i in M_P'`.

If a nonzero `d in Delta(P)` belonged to `N_sep(P)M`, every basis coordinate
of `d` would be divisible by `N_sep(P)`.  But (1.3) makes every coordinate
strictly smaller in absolute value than `N_sep(P)`, while at least one is
nonzero.  This is impossible.  Hence none of the forbidden differences
vanishes modulo `M_P`.

Any same-norm sheet in the class `mu_i` which meets `P` is represented by a
vector `v in U_T` satisfying `v-r_i in M_P` (or the corresponding sign
relation).  The preceding paragraph forces `v^perp=r_i^perp`.  The same
argument separates distinct target facets modulo sign.  This proves all
claims.  `square`

The formula is basis-dependent but the conclusion is not.  Minimizing (1.3)
over reduced bases is an optional computational improvement, not part of the
theorem.

## 3. Discriminant-level growth

Let `ell(M)` be the level of the finite quadratic module `M'/M`, using the
quadratic-form convention of the Weil representation.  Scaling the lattice
gives

\[
 (NM)'=N^{-1}M'.
\]

For `x=y/N` with `y in M'`,

\[
 N^2\ell(M)Q(x)=\ell(M)Q(y)\in\mathbf Z.
\]

Consequently

\[
 \boxed{\ell(NM)\mid N^2\ell(M)}.
\tag{3.1}
\]

In particular the exact-support construction may use

\[
 \ell_P\leq N_{\rm sep}(P)^2\ell(M)
\tag{3.2}
\]

after replacing the actual level by this explicit clearing multiple.  This
is the level parameter which enters any quantitative diagonal or
orbit-square sign-change estimate.

## 4. What remains uniform

Equations (1.3) and (3.2) give a terminating, exact `P`-first algorithm:

```text
marked compact P
  -> enumerate all target-norm mirrors meeting P
  -> compute N_sep(P)
  -> exact facet-only discriminant classes
  -> an explicit clearing level for the diagonal Weil representation.
```

They do **not** bound `N_sep(P)` using only the dimension, field, discriminant,
or facet count.  Such a class-uniform estimate requires a quantitative bound
for the coordinates of all target-norm mirrors meeting `P`, or a different
level construction.  Nor does (3.2) bound the Borcherds weight or the top
automorphic degree of the resulting cover.

Thus this theorem closes the per-input level algorithm while leaving the
uniform Missing-Lemma-A+ problem intact.

