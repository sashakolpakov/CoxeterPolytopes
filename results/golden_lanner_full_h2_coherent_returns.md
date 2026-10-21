# Full-`H^2` coherent K3 returns for all twelve golden Lanner chambers

**Status (2026-08-09): exact rational theorem.**  Every member of the
exhaustive equal-root golden Lanner class now has a coherent

```text
P -> projective K3 with full-H2 Q(sqrt(5)) action -> the same P
```

return.  This upgrades all twelve previously certified divisor-only
eigen-nef returns, rather than only Bugaenko's `H8` endpoint.  The abstract
polytopes are the classical Lanner simplices; the new assertion is the
full-weight-two Hodge realization.  "New" here means new to this project;
no worldwide priority claim is made without a dedicated literature audit.

The action is rational on `H^2`.  It automatically preserves some
finite-conductor order after clearing denominators, but this report does
**not** assert that the maximal order acts integrally on the K3 lattice,
compute the conductor of the preserved order, or represent the Hodge
operator by an algebraic correspondence in all twelve rows.

## 1. Result and finite ledger

Let `F=Q(sqrt(5))`, `O=Z[phi]`, `phi^2=phi+1`, and

\[
 \alpha=-\frac{\phi}{\sqrt5}.
\]

For a certified field Gram `H` of rank `r=n+1`, write `d=det(H)` and
`m=11-r`.  The following diagonal `F`-form supplies the transcendental
block:

\[
 C_t=\langle-\phi t,\ d/t,\ 1^{m-2}\rangle . \tag{1.1}
\]

The exact successful twists are:

| chamber | dimension | `rk NS` | `m=dim_F C_t` | twist `t` | verdict |
|---|---:|---:|---:|---|---|
| `(2,5,5)` | 2 | 6 | 8 | `1` | coherent same-`P` return |
| `(3,3,5)` | 2 | 6 | 8 | `2` | coherent same-`P` return |
| `(3,5,5)` | 2 | 6 | 8 | `2+phi` | coherent same-`P` return |
| `(5,5,5)` | 2 | 6 | 8 | `2+phi` | coherent same-`P` return |
| `[3,5,3]` | 3 | 8 | 7 | `1` | coherent same-`P` return |
| `[5,3,5]` | 3 | 8 | 7 | `1` | coherent same-`P` return |
| `D4(5,3,3)` | 3 | 8 | 7 | `1` | coherent same-`P` return |
| `C4(one 5)` | 3 | 8 | 7 | `2+phi` | coherent same-`P` return |
| `C4(opposite 5s)` | 3 | 8 | 7 | `3` | coherent same-`P` return |
| `[5,3,3,3]` | 4 | 10 | 6 | `1` | coherent same-`P` return |
| `[5,3,3,5]` | 4 | 10 | 6 | `1` | coherent same-`P` return |
| `D5(terminal 5)` | 4 | 10 | 6 | `1` | coherent same-`P` return |

There are no undecided or failed rows in this class.

## 2. Exact rational complement certificate

The trace form

\[
 q_t=\operatorname{Tr}_{F/\mathbf Q}(\alpha C_t)
\]

has rational rank `2m`.  The field determinant `d` is negative at the
distinguished embedding and positive at the conjugate embedding.  Every
listed `t` is totally positive.  Moreover

\[
 \alpha=-\frac{\phi+2}{5},\qquad
 \alpha'=-\frac{3-\phi}{5},
\]

so the trace scalar is negative at both real places.  The field form `C_t`
has two negative directions at the distinguished place and is positive
definite at the conjugate place.  After multiplication by `alpha`, the two
real eigenspaces of multiplication by `phi` on `q_t` therefore have
signatures

\[
 (2,m-2)\quad\hbox{and}\quad(0,m). \tag{2.1}
\]

In particular the distinguished eigenspace contains positive oriented
two-planes and the conjugate eigenspace is negative definite.

For each row the certifier forms

\[
 N(H)_\mathbf Q\perp q_t
\]

and compares it with

\[
 \Lambda_{K3,\mathbf Q}=U^3\perp E_8(-1)^2
\]

by the complete Hasse--Minkowski invariants:

1. rational rank `22`;
2. real signature `(3,19)`;
3. determinant square class `-1`; and
4. the Hasse invariant at infinity and every finite prime occurring in an
   exact rational diagonalization.

This finite-prime list is complete: at an omitted odd prime every diagonal
coefficient is a unit, so every pairwise Hilbert symbol is `+1`; the prime
`2` is always included.  Hilbert reciprocity is also checked exactly.  All
invariants agree in all twelve rows.  This is a proof of rational isometry,
not a numerical Gram comparison.  Multiplication by `phi` on both the frozen
Picard trace block and `C_t` is checked exactly to be self-adjoint and to
satisfy

\[
 \tau^2=\tau+I. \tag{2.2}
\]

Now choose any primitive K3 embedding of the integral Picard trace lattice
`N(H)` supplied by the earlier stable embedding theorem, and let `T` be its
orthogonal complement.  We have

\[
 N(H)_\mathbf Q\perp T_\mathbf Q\simeq\Lambda_{K3,\mathbf Q}
 \simeq N(H)_\mathbf Q\perp q_t.
\]

Witt cancellation gives `T_Q isometric to q_t`.  This holds for **every**
fixed primitive embedding: cancellation uses the same abstract rational
quadratic summand `N(H)_Q` and does not require a chosen total-space isometry
to restrict to the identity on it.  Transporting multiplication by `phi`
therefore extends the already certified Picard action to a unital faithful
self-adjoint `F`-action on all of `H^2_Q`.  The certifier rechecks the
minimal-polynomial and self-adjoint equations on the full 22-dimensional
orthogonal-sum model.  No choice of a different Picard lattice is hidden in
this step.

### Automatic finite-conductor integral order

Although the rational isometry is existential, it also gives an exact
existence theorem at the integral level.  Write

\[
 \rho:F\longrightarrow\operatorname{End}_{\mathbf Q}(\Lambda_{K3,\mathbf Q})
\]

for the transported action.  In any integral basis of `Lambda_K3`, the
matrix `rho(phi)` has finitely many rational denominators.  Choose a positive
integer `D` clearing all of them.  Then `D rho(phi)` is integral and

\[
 \mathcal O_D=\mathbf Z+D\mathcal O_F
             =\mathbf Z[D\phi]
\]

acts unitally on the integral K3 lattice.  Closure is literal because

\[
 (D\phi)^2=D(D\phi)+D^2.
\]

Thus every row carries an integral action by **some** order of finite
conductor.  The certificate does not construct the rational isometry, so it
does not freeze the least `D`; and it does not promote this to the maximal
order `O_F`.

### Kuga--Satake perk, with a strict boundary

Put `B=2 tau-I`, where `tau=rho(phi)`.  Self-adjointness and (2.2) give

\[
 B^2=5I,\qquad (Bx,By)=5(x,y). 
\]

Thus `B` restricts to a Hodge similarity of multiplier `5` on the
transcendental K3-type Hodge structure.  Varesco's Proposition 3.1 in
[*Hodge similarities, algebraic classes, and Kuga--Satake varieties*](https://arxiv.org/abs/2304.02519)
therefore supplies a self-isogeny of the associated Kuga--Satake abelian
variety compatible with the Kuga--Satake embedding.  The script checks
`B^2=5I` and the multiplier equation exactly.

This conclusion is only a Kuga--Satake isogeny.  Without the additional
Kuga--Satake Hodge-conjecture hypothesis used in Varesco's later
algebraicity result, it does **not** prove an algebraic correspondence on
`X x X`.  The isogeny is also a by-product of the Hodge return, not a carrier
or same-`P` gate.

### Four special-period algebraic upgrades

There is nevertheless a separate, unconditional `X x X` upgrade in four
rows.  For `(2,5,5)`, `[5,3,3,3]`, `[5,3,3,5]`, and `D5(terminal 5)`, an
explicit `K`-linear skew-adjoint operator `J`, with `J^2=-c`, enlarges the
transcendental action to a quartic CM field.  A generic period in the
resulting complex ball has exact Picard lattice and exact CM endomorphism
field.  Ramón-Marí's theorem then algebraicizes the rational Hodge action as
a codimension-two correspondence on `X x X`.  The complete proof, matrices,
and the eight sharp quadratic-extension failure gates are in
`results/golden_lanner_cm_algebraic_upgrade_audit.md`.

Thus the Kuga--Satake statement holds uniformly in all twelve rows, while a
surface-square algebraic correspondence is presently certified in four.

## 3. Hodge period, exact Picard lattice, and return of the same chamber

Choose a generic oriented positive two-plane in the distinguished
eigenspace in (2.1).  It satisfies the K3 period equations and is scalar for
the `F`-action, so the transported action is a Hodge action.  The complement
has `F`-rank `6`, `7`, or `8`.  A nonzero rational vector of `T_Q` cannot lie
in only one real eigenspace: otherwise a rational matrix would act on a
nonzero rational vector by the irrational scalar `phi` or `phi'`.  Hence the
orthogonality condition for every rational class outside `N(H)_Q` is a
proper locus in the distinguished period domain.  Avoiding their countable
union leaves a period with

\[
 \operatorname{NS}(X)=N(H)
\]

exactly.  Since `N(H)` is hyperbolic, the resulting K3 is projective.

The special Hodge-coherent period changes neither the integral Picard
lattice nor its complete restricted `(-2)`-root arrangement.  The upstream
certificate already proves for every row that:

- the original simplex is compact by the finite-vertex form of Vinberg's
  criterion;
- no norm-two mirror cuts its interior;
- the complete restricted K3 root arrangement is the field norm-two
  arrangement; and
- K3 Weyl framing returns the original metric simplex as the eigen-nef
  section.

Therefore the new full-`H^2` action returns the **same** `P`, not a
subdivision or an unspecified chamber.

## 4. The theorem boundary supplied by the same class

The upstream exhaustive census proves that these twelve rows are the whole
connected equal-root golden Lanner class: four in `H2`, five in `H3`, three
in `H4`, and none in `H5` or above.  Combining that classification with the
present complement theorem gives a positive/negative boundary package in
one finite statement:

\[
 \boxed{
 \begin{array}{c}
 \text{exactly twelve coherent same-}P\text{ K3 returns in the class},\\
 \text{and no higher-dimensional member of the class.}
 \end{array}}
\]

This is not a bound for arbitrary compact arithmetic Coxeter polytopes.
It is a complete positive/negative criterion for a sharply defined,
surface-free arithmetic simplex class.  Crucially, the `H5+` nonexistence
statement comes from Lanner's Coxeter classification upstream; it is **not**
an algebraic-geometric inverse bound derived from these K3 carriers.  The
new theorem on this page is the coherent positive return for each of the
twelve survivors.

## 5. Scope and lineage

- **Positive content:** twelve projective K3 realizations carrying a
  rational self-adjoint `Q(sqrt(5))` Hodge action on all of `H^2`, with exact
  Picard lattice and exact same-polytope return.
- **Not new abstract polytopes:** all twelve Coxeter simplices are classical
  Lanner diagrams.
- **Not an automorphism claim:** multiplication by `phi` is not an isometry
  of the intersection form.
- **Kuga--Satake but not `X x X`:** `B=2tau-I` gives a multiplier-five Hodge
  similarity and hence, by Varesco, a compatible Kuga--Satake self-isogeny.
  This alone does not algebraicize the full rational Hodge endomorphism as a
  correspondence on `X x X` and plays no carrier role.  A separate CM
  refinement does algebraicize four rows.
- **Integral-order boundary:** clearing denominators proves existence of an
  integral action by some finite-conductor order `Z+D O_F`.  Neither the
  least conductor nor maximal-order integrality is certified.

Thus the project lineage gain is a coherent realization theorem for an
entire known polytope class, with honest rational/order, Hodge/algebraic, and
positive/negative-bound boundaries.

## 6. Reproduction

- Same-`P` input theorem and certificate:
  `results/golden_lanner_exhaustive_census.md` and
  `data/enumerations/golden_lanner_exhaustive_census.json`.
- Full-`H^2` certifier:
  `scripts/certify_golden_lanner_hodge_coherent_carriers.py`.
- Tests:
  `tests/test_golden_lanner_hodge_coherent_carriers.py`.
- Quadratic-CM/algebraic upgrade audit:
  `results/golden_lanner_cm_algebraic_upgrade_audit.md`, with certifier
  `scripts/certify_golden_lanner_cm_upgrade.py` and tests
  `tests/test_golden_lanner_cm_upgrade.py`.

Run

```bash
PYTHONPATH=src:. python3 scripts/certify_golden_lanner_hodge_coherent_carriers.py
PYTHONPATH=src:. python3 -m pytest -q tests/test_golden_lanner_hodge_coherent_carriers.py
```

The optional bounded discovery replay

```bash
PYTHONPATH=src:. python3 scripts/certify_golden_lanner_hodge_coherent_carriers.py --search-bound 5
```

independently rediscovers every frozen nontrivial twist.
