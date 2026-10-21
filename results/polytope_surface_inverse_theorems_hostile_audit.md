# Hostile audit of the polytope-to-surface inverse theorems

**Date:** 2026-08-09  
**Verdict:** **QUALIFIED PASS after repair.**

This audit attacks the two inverse-program notes
[`polytope_to_eigencarrier_hodge_rank_return.md`](polytope_to_eigencarrier_hodge_rank_return.md)
and
[`universal_ns_realization_and_nef_leverage_audit.md`](universal_ns_realization_and_nef_leverage_audit.md).
The rank obstruction, low-rank K3 realization theorem, and exact same-cone
criterion survive.  Four qualifications were necessary and have been patched:

1. common doubling preserves the integrality status of a reflection; it does
   not turn a nonintegral reflection into an integral one;
2. the integral field action is by the specified order stabilizing the
   lattice, not automatically by the maximal order;
3. arithmetic preservation belongs to the marked order-module before taking
   one real eigenspace, which need not contain a discrete full-rank
   `Z`-lattice; and
4. the same-cone statement needs closed-cone and regular-slice conventions,
   excluding an eigenspace contained in a complete wall.

No dimension bound in either report is asserted for arbitrary Coxeter
polytopes without an exhaustive carrier compiler.

## 1. Morrison rank eleven and exact Picard realization

The primary-source statement is exact:

- Morrison, Corollary 2.9(i): for `rho<=10`, every even lattice of signature
  `(1,rho-1)` occurs as the Neron--Severi lattice of an algebraic K3 surface,
  and its primitive embedding into the K3 lattice is unique;
- Morrison, Remark 2.11: every even lattice of signature `(1,10)` also occurs
  as the Neron--Severi lattice of an algebraic K3 surface, but no uniqueness
  statement is supplied without extra hypotheses.

Thus the report's `r<=11` **existence** theorem and `r<=10` uniqueness cutoff
are correct.  If `B` is integral, `2B(x,x)` is even for every lattice vector,
so `L(2)` is even.  Morrison realizes `L(2)` as the whole `NS`, not merely as
a sublattice.  The primary reference is D. Morrison, *On K3 surfaces with
large Picard number*, Invent. Math. 75 (1984), 105--121,
[DOI](https://doi.org/10.1007/BF01403093), Corollary 2.9 and Remark 2.11.

The marked period dimension is `20-r`: the period domain in
`P(N^perp_C)` is a quadric of dimension `(22-r)-2`.  A generic period avoids
the countable rational Noether--Lefschetz hyperplanes, hence has `NS=N`.
Because `N` has signature `(1,r-1)`, the resulting K3 is projective.  The
same calculation for a fixed irreducible-holomorphic-symplectic deformation
type is `b2-r-2`.  At `r=b2-2`, `N^perp_R` is exactly the positive period
plane; primitiveness then forces the integral `(1,1)` lattice to be exactly
`N`, despite the period locus being zero-dimensional.

## 2. Scaling and reflection integrality

For every positive scalar `c`,

```text
s_v(x)=x-2 B(x,v)/B(v,v) v
      =x-2 cB(x,v)/(cB(v,v)) v.
```

Therefore common scaling preserves mirrors, oriented halfspaces, normalized
Gram products, angles, ultraparallel distances, and the real reflection map.
Since `L(2)` has the same underlying abelian group as `L`, it also preserves
the yes/no answer to `s_v(L)=L`.  The original wording risked suggesting that
doubling itself supplied integrality; the repaired theorem assumes it when
an integral-reflection conclusion is used.

Likewise, if a chosen order `A` acts on `L`, every `a in A` remains an
integral self-adjoint endomorphism after doubling.  An element of
`O_K\setminus A` need not preserve `L`.  The maximal-order version is valid
only for an `O_K`-stable lattice.

## 3. Hodge orthogonality and the rank return

The Hodge argument is sound under the now-explicit `A`-stability hypothesis.
The order action extends to `K` on `N_Q`, and

```text
K tensor_Q R = product over sigma of R
```

gives `N_R=direct_sum E_sigma`, each summand of dimension `n+1`.  If
`sigma!=tau`, self-adjointness and an element separating the two embeddings
give `(E_sigma,E_tau)=0`.  The distinguished `E_0` is already Lorentzian of
signature `(1,n)`.  Its orthogonal complement in a space of positive index
one is negative definite, so all conjugate eigenspaces and `N^perp` are
negative definite.  Hence

```text
rank_Z N = [K:Q](n+1) <= rho(X).
```

Trace duality then uniquely recovers the `K`-bilinear form `h` from
`B(ax,y)=Tr(a h(x,y))`; self-adjointness gives bilinearity in both variables.
This proves the claimed admissibility at the conjugate real places.  It is a
necessary condition for this faithful trace/eigencarrier mechanism, not for
every conceivable unrelated surface package.

## 4. Carrier-rank table

The table is arithmetically consistent:

| type | `b2` or numerical rank | divisor-rank ceiling | quadratic `n` ceiling |
|---|---:|---:|---:|
| K3 | `22` | `20` | `9` |
| Enriques | `Num` rank `10` | `10` | `4` |
| `K3^[k]`, `k>=2` | `23` | `21` | `9` |
| generalized Kummer, dimension `>=4` | `7` | `5` | `1` |
| OG6 | `8` | `6` | `2` |
| OG10 | `24` | `22` | `10` |

For the irreducible-symplectic rows, the ceiling is `rho<=b2-2`; Enriques is
separate because its numerical divisor lattice itself has rank ten.  Applying
`2(n+1)<=R` gives exactly the displayed `H9,H4,H9,H1,H2,H10` list.  These are
Hodge-slot bounds only; existence, wall completeness, and finite volume are
independent gates.

## 5. Exact same-cone criterion

After the repaired conventions, the criterion is a literal two-inclusion
argument.  Write the closed geometric chamber as the intersection of **all**
its oriented geometric wall halfspaces.  If every input facet is one of the
restricted complete walls, then `A_X intersect E_0` is contained in `P+`.
If every complete wall functional is nonnegative on `P+`, the reverse
inclusion holds.  Completeness of the wall set is indispensable.

The regularity hypothesis rules out a wall containing the entire slice.  The
closed positive sheet includes isotropic rays and the origin, so the equality
also covers noncompact finite-volume chambers.  Equality of cones alone does
not prove equality of Coxeter data: the restricted oriented Gram and
preservation of the marked order-module must still be checked.

For a full-Picard K3, all square-minus-two roots supply the complete mirror
arrangement and Weyl transitivity permits re-marking the chamber as nef.
For an irreducible symplectic carrier, the analogous assertion requires the
complete relevant MBM or prime-exceptional arrangement and the appropriate
monodromy/chamber theorem; bounded wall squares are not by themselves a
complete wall census.

## 6. Counterexamples and larger reservoirs

The three-point blowup test is correct.  Both marked surfaces have

```text
NS=I_(1,3),  K=-3H+E1+E2+E3,
K^2=6, c2=6, chi(O)=1, q=pg=0.
```

For noncollinear points `-K` is ample.  For collinear points the strict
transform `R=H-E1-E2-E3` is effective, has `R^2=-2`, and satisfies
`(-K).R=0`; it supplies a genuine extra nef/effective wall.  Thus the marked
pair `(NS,K)` and the displayed numerical invariants do not determine the
nef cone.

For a product of curves the two fibre classes are nef and isotropic and meet
once, so the full nef cone has ideal boundary rays and cannot be a compact
direct full-nef carrier.  The remaining Neron--Severi term is governed by
Jacobian homomorphisms, not an arbitrary negative lattice.

For a finite surjective map `f:X->Y`, projection formula gives
`(f^*D_1.f^*D_2)=deg(f)(D_1.D_2)`, and `f^*D` is nef exactly when `D` is nef:
the reverse direction follows by taking a component above each curve of
`Y`.  Hence

```text
Nef(X) intersect f^*NS(Y)_R = f^*Nef(Y).
```

Covers transport a controlled cone; additional divisor classes upstairs
create another wall-completeness problem.  Neither products nor covers repair
the inverse compiler gap.

## 7. Residual scope

The repaired reports prove genuine implications **inside stated carrier
classes**.  They do not prove that every compact arithmetic Coxeter polytope
has a K3, Enriques, or irreducible-symplectic eigencarrier, and failure of one
primitive embedding excludes only that carrier species.  The universal
rational-blowup construction supplies weak semialgebraic divisor slices for
its admissible lattice-form class, not effective facets or a complete nef
cone.  Consequently no universal compact-polytope dimension theorem follows
until a surface-free polytope class is given an exhaustive strong compiler.

With those quantifiers retained, the two inverse reports pass.
