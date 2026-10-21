# Hostile audit of the maximal compact OG10 return compiler

## Verdict

**Qualified pass after five precision repairs.**  No theorem-breaking defect
was found in the determinant/gluing argument, the rank-two period return, the
four OG10 wall types, or the compact ellipsoid bound.  For every **fixed
marked** datum in the stated compact quadratic `H10` input class, existence of
an exact-Picard OG10 carrier whose Kähler-chamber closure returns the same
polytope is a finite decision problem.

The pass does not mean that the repository currently contains an end-to-end
implementation of every finite step.  The executable certifier is the
arithmetic front end.  Nor does the theorem enumerate all abstract compact
arithmetic `H10` polytopes, produce a compact `H10` survivor, or prove an
NS-only `n<=9` bound.

The audited theorem is

\[
 \mathbf P\bigl(\overline{\mathcal K_X}\cap E^+\bigr)=P,
 \qquad
 \mathbf P\bigl(\mathcal K_X\cap E^+\bigr)=\operatorname{int}P.
\]

Calling the open Kähler cone itself the closed polytope would be false.

## Gate-by-gate audit

| gate | verdict | audit conclusion |
|---|---|---|
| Surface-free quantifier | pass | The input fixes a marked real-quadratic order-module, its even rank-22 trace lattice, the distinguished `H10` eigenspace, and an exactly certified compact Coxeter polytope.  Failure applies to that integral realization, not automatically to its abstract diagram. |
| Primitive complement and determinant | pass | For primitive `N -> Lambda_OG10`, the positive complement has rank two and `det(T)=3k^2/|det(N)|`; injectivity of both gluing projections gives `k | |A_N|` and `k | |A_T|`.  Hence only finitely many `k` and complement determinants occur. |
| Marked embedding orbits | repaired/pass | Quotienting by all of `O(N) x O(T)` is valid only after forgetting the marked domain.  The compiler must retain every discriminant graph, or quotient only by `O(T)` and the stabilizer of all marked `N` data.  The finite-graph enumeration remains exhaustive. |
| Positive binary forms | pass | Every positive even binary lattice has a reduced Gram `[[2a,b],[b,2c]]` with `0<=b<=a<=c`.  The four frozen determinant lists are complete; the script also performs an exact finite `GL(2,Z)` isometry search and finds no boundary duplicates. |
| Local length filter | pass | Reduction modulo `p` bounds `dim Rad(N/pN)` by corank two plus `dim Rad(Lambda/pLambda)`.  This gives `ell_p(A_N)<=2` for `p!=3` and `ell_3(A_N)<=3`. |
| Three-coprime lane | pass | If `3` does not divide `det(N)`, valuation and projection injectivity force `k=|det(N)|` and `det(T)=3|det(N)|`.  Nondegeneracy of `N/3N` excludes primitive Picard vectors of ambient divisibility three, so only the `(-2,1)` and `(-4,1)` OG10 types remain. |
| Maximal-rank periods | pass | A positive rank-two complement has exactly two oriented period lines.  Their real and imaginary parts span all of `T_R`, so the integral `(1,1)` lattice is exactly the primitive `N`.  A positive integral class then gives projectivity.  Orientation/marked-moduli choices are finite. |
| OG10 wall alphabet | pass | Mongardi--Onorati Proposition 5.4 gives exactly `(-2,1),(-4,1),(-6,3),(-24,3)`.  The first and third are the stably-prime-exceptional types; all four are needed for Kähler chambers. |
| Compact ellipsoid | pass | If a wall of square `s` meets `P`, then `0<-delta_E^2<=-s`.  Combining the hyperplane-distance bound with the positive-definite form `Qhat_(p0)(delta)=2(p0,delta)^2/(p0,p0)-delta^2` gives `Qhat_(p0)(delta)<=(-s)(2cosh^2(R)-1)`.  Thus every wall meeting compact `P` lies in a bounded positive-definite ellipsoid and can be enumerated exactly. |
| Regular slice point | repaired/pass | One may start the ellipsoid proof from any exact strict interior point.  After the finite relevant-wall list is known, choose an exact interior point off those mirrors.  The faithful quadratic module has no nonzero rational vector in only one real eigenspace, so no nonzero wall vanishes identically on `E`; a generalized package where this happens must be rejected. |
| Same-chamber iff | pass | Facet survival gives the complete chamber inside `P`; nonnegativity of every retained wall on `P` gives the reverse inclusion.  A wall outside the ellipsoid has no zero on connected `P`, hence has the orientation sign throughout and cannot add a hidden cut. |
| Realization of every selected component | repaired/pass | Mongardi--Onorati Theorem 5.5 gives the numerical wall components.  Conversely, a very general point of any component lies in a Kähler-type chamber by Markman Lemma 5.11, and Markman Theorem 5.16 realizes every such chamber in the same period fibre.  Onorati's `Mon^2=O^+` theorem controls the remaining marked orientation component. |
| Movable statement | repaired/pass | Using only `(-2,1),(-6,3)` computes a prime-exceptional/birational-Kähler chamber and its closure, not a Kähler chamber.  In closed-cone language this is the movable closure; the `-4` and `-24` walls return at the Kähler/nef level. |
| Executable scope | repaired/pass | The current script certifies Smith data, local lengths, determinant lanes, binary-form censuses, and the four frozen replays.  Discriminant anti-isometry enumeration and the compact ellipsoid wall enumeration are proved terminating but are not implemented end to end in this script. |

## Independent arithmetic replay

The certifier returns the following exact front-end data.

| trace determinant | Smith nonunits | surviving `(k,det T,#forms)` lanes |
|---:|---|---|
| `44` | `(2,22)` | `(44,132,4)` |
| `20` | `(2,10)` | `(20,60,4)` |
| `36` | `(6,6)` | `(12,12,2)`, `(36,108,5)` |
| `4` | `(2,2)` | `(4,12,2)` |

The selected complements `<2>+<66>`, `<2>+<30>`, `<2>+<6>`, and `A2(2)`
occur in their listed lanes.  This replay does not claim that every candidate
form or determinant lane survives discriminant anti-isometry.

Focused verification:

```text
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_og10_h10_maximal_return_compiler.py

3 passed in 0.79s
```

The standalone certifier also completed, emitted valid JSON, and reported
all eight frozen checks true.  Its canonical payload hash in this run was

```text
029ec46d3667be76e006680d2eddfb7245ea556bb72446d58603aa212c8a133c
```

## Exact remaining scope

The audit does **not** establish any of the following.

1. A finite P-first list of all abstract compact arithmetic quadratic `H10`
   diagrams, normal scalings, order lattices, and saturations.
2. Essential surjectivity from all compact arithmetic `H10` polytopes into
   this OG10 carrier class.
3. An implemented end-to-end discriminant-glue plus ellipsoid wall compiler.
4. A compact or finite-volume new `H10` Coxeter polytope.  The two norm-two
   calibrations remain infinite-volume, the norm-six carrier has open
   complete-chamber/volume status, and the `D4` package has a hidden wall.
5. An NS-only universal `n<=9` theorem.  The `H9` ceiling applies only to the
   separately defined coherent full-`H2` branch.

Thus the exact achievement is a terminating **carrier-relative return
criterion for each fixed marked compact input**, which is the finite inner
engine needed once a P-first denominator/essential-surjectivity theorem is
supplied.

## Sources and audited files

The four wall types and Kähler-cone description are
[Mongardi--Onorati, Proposition 5.4 and Theorem 5.5](https://arxiv.org/abs/2010.12511).
The period-fibre/Kähler-type chamber realization is
[Markman, Lemma 5.11 and Theorem 5.16](https://arxiv.org/abs/1101.4606).
The maximal OG10 monodromy statement is
[Onorati, Theorem 5.4](https://arxiv.org/abs/2002.04129).

Audited local files:

- [maximal compiler theorem](og10_h10_maximal_return_compiler.md)
- [`certify_og10_h10_maximal_return_compiler.py`](../scripts/certify_og10_h10_maximal_return_compiler.py)
- [`test_og10_h10_maximal_return_compiler.py`](../tests/test_og10_h10_maximal_return_compiler.py)
