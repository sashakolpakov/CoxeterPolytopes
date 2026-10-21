# Hostile audit: full-`H^2` coherent golden-Lanner returns

**Audit status (2026-08-09): PASS after arithmetic, quantifier, and scope repairs.**

Audited objects:

- [`golden_lanner_full_h2_coherent_returns.md`](golden_lanner_full_h2_coherent_returns.md);
- [`certify_golden_lanner_hodge_coherent_carriers.py`](../scripts/certify_golden_lanner_hodge_coherent_carriers.py); and
- [`test_golden_lanner_hodge_coherent_carriers.py`](../tests/test_golden_lanner_hodge_coherent_carriers.py).

The twelve input chambers and their complete-root same-polytope Picard
returns are imported from the upstream exhaustive golden-Lanner certificate.
This audit does not reclassify Lanner diagrams or infer root completeness
from a finite wall list.  It audits the new step from each fixed Picard
return to a coherent action on the full rational K3 `H^2`.

## Arithmetic verdict

Let `F=Q(sqrt(5))`, let `H` be a frozen golden field Gram of rank `r`, and
put `m=11-r` and `d=det(H)`.  For the listed totally positive twist `t`, the
certificate uses

```text
C_t = <-phi*t, d/t, 1, ..., 1>
q_t = Tr_{F/Q}((-phi/sqrt(5)) C_t).
```

The sign bookkeeping is correct.  The determinant `d` has signs `(-,+)` at
the two real places, whereas both conjugates of
`-phi/sqrt(5)=-(phi+2)/5` are negative.  Consequently `C_t` has signatures
`(m-2,2)` and `(m,0)`, and the two real eigenspaces of `q_t` have signatures
`(2,m-2)` and `(0,m)`.  Thus the distinguished eigenspace, and only that
eigenspace, can contain the K3 period plane.

For every row the code diagonalizes the rational form exactly and compares
rank, real signature, determinant square class, and Hasse invariant with
`U^3 + E8(-1)^2`.  The finite-place list is complete: it includes `2` and
every prime dividing a diagonal numerator or denominator; at every omitted
odd prime all coefficients are units and every Hilbert symbol is `+1`.
Hilbert reciprocity is now checked as an internal consistency condition.
Independent PARI/GP spot checks agreed with the implemented odd and dyadic
Hilbert-symbol formulas.  Hence the use of Hasse--Minkowski is exhaustive,
not a finite heuristic.

## Embedding, action, and period verdict

For any fixed primitive embedding of the integral Picard lattice `N` in the
K3 lattice, with orthogonal complement `T`, the two decompositions

```text
N_Q + T_Q  ~=  Lambda_K3,Q  ~=  N_Q + q_t
```

share the same abstract nonsingular rational summand `N_Q`.  Witt
cancellation therefore gives `T_Q ~= q_t`.  This does not require choosing
an isometry of the total spaces that is already the identity on `N_Q`, so
the claim really does hold for every fixed primitive embedding.

Multiplication by `phi` on the Picard trace block is checked against the
upstream frozen matrix.  Its transported complement action and the resulting
block-diagonal rank-22 action are checked exactly for self-adjointness and
the relation `tau^2=tau+I`.  The latter makes the `F`-action unital and
faithful.

The exact-Picard assertion is also valid.  A nonzero rational vector in the
transcendental block cannot lie wholly in either irrational real eigenspace:
that would make a rational matrix have `phi` or its conjugate as an
eigenvalue on a rational eigenvector.  Each unwanted rational Hodge class
therefore cuts out a proper Noether--Lefschetz locus in the distinguished
period domain.  Avoiding their countable union gives `NS(X)=N` exactly.
The already certified complete Picard root arrangement then returns the same
metric simplex; the Hodge specialization introduces no new divisor wall.

## Integral-order and Kuga--Satake boundaries

A rational `F`-action on the integral K3 lattice automatically stabilizes
some finite-index order.  If `D` clears the denominators of `rho(phi)`, then

```text
O_D = Z + D O_F = Z[D phi]
```

acts integrally, since `(D phi)^2=D(D phi)+D^2`.  This is an existence
statement only: the certificate does not construct the rational isometry,
find the least conductor, or prove integral action by the maximal order.

For `B=2 tau-I`, the exact identities

```text
B^2 = 5 I,                 B^T q B = 5 q
```

make the restriction to the transcendental K3-type Hodge structure a Hodge
similarity of multiplier `5`.  Varesco, Proposition 3.1, consequently gives
a compatible self-isogeny of the Kuga--Satake abelian variety.  It does
**not** by itself give an algebraic correspondence on `X x X`; that stronger
claim needs an additional Kuga--Satake Hodge-conjecture hypothesis.  The
isogeny is a by-product and is not used as a carrier or same-polytope gate.

The main report now also points to a **separate** quadratic-CM audit that
algebraicizes four specially chosen periods.  That result does not follow
from `B` or from Varesco's Proposition 3.1 and is outside the proof audited
here.  Accordingly, the base certifier's
`X_times_X_algebraic_correspondence_certified = false` means “not certified
by this full-`H^2`/Kuga--Satake certificate”; it does not negate the distinct
four-row CM certificate.

## Repairs made by this audit

- The optional bounded twist search previously never ran because every row
  already had a frozen twist.  It now independently searches the bounded
  candidate set and requires the discovered twist to equal the frozen one.
- The certifier now checks the imported Picard action, the complement action,
  and the full rank-22 action, including minimal polynomials and
  self-adjointness.
- Sign, rank, complete-place, and Hilbert-reciprocity assertions are now
  executable gates rather than explanatory prose.
- The aggregate success flag now requires every row to pass; it no longer
  succeeds through a default value when a row is incomplete.
- The period argument, fixed-embedding form of Witt cancellation,
  finite-conductor order, and maximal-order caveat are explicit.
- The `B` similarity is certified, while `X x X` algebraicity and any use of
  Kuga--Satake as a carrier are explicitly recorded as false.
- The report and ledger now distinguish a positive coherent-return upgrade
  from a new abstract polytope or an algebraic-geometric inverse bound.

## Replayed tests

```text
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_golden_lanner_hodge_coherent_carriers.py \
  tests/test_golden_lanner_exhaustive_census.py

10 passed
```

The independent `--search-bound 5` replay rediscovers the frozen twist in
all twelve rows.  Python bytecode compilation also passes.

## Scope and novelty verdict

This is a genuine positive theorem: all twelve members of the exhaustive
equal-root golden Lanner class admit coherent rational full-`H^2`
`Q(sqrt(5))` K3 returns.  It is not a source of new abstract Coxeter
polytopes: those simplices are classical.  Nor is the absence of an `H^5`
row an algebraic-geometric consequence of the carrier; it is imported from
the upstream Lanner classification.  The result is therefore a positive
coherent refinement of a complete known class, not yet the desired
polytope-to-surface obstruction yielding a new universal dimension bound.
