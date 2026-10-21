# Independent review of the height-charge dimension-165 candidate

**The new local estimates and their global accounting pass this independent
review.** No mathematical discrepancy was found. The resulting research proof
candidate gives dimension at most 165 for finite-sided, finite-volume real
hyperbolic Coxeter polytopes, conditional only on the explicitly inherited
geometric chain already reviewed in the preceding packages. This is an
internal mathematical and exact-arithmetic review, not a formal verification,
publication-priority claim, or independent reproduction of all historical
geometric inputs.

The [primary proof](fresh_nikulin_height_charge_20261015.md),
[independent checker](../scripts/fresh_nikulin_height_charge_20261015_review.py),
and [review certificate](fresh_nikulin_height_charge_20261015_review.json)
are separate artifacts. The checker uses standalone `Fraction` arithmetic,
grounded matrix inverses, rational convex-dual verification, and polynomial
coefficient arithmetic. It imports **no project or producer code, SymPy,
optimizer, or floating-point solver**. It binds the final primary PY/JSON,
not its report.

## Exceptional factors: exhaustive support coverage

The checker independently builds the literal affine A1/E6/E7/E8/F4/G2
matrices. For the last two it uses rationally scaled root bases; the physical
norm-two Gram is checked entry by entry through signs and squares. Grounding
one node gives a rational principal inverse, from which the positive kernel,
normalized kernel squares, and point resistances are recovered. Every
principal deletion is positive definite, and every retained inverse is
entrywise nonnegative. The ordinary and forbidden-node omission rules are
reconstructed explicitly, rather than copied from the new record list.

All 40 records agree with the frozen preset choices. The zero-retained-contact
case is checked in each record, including the coincident-omission A1 case,
whose target is zero.

Write \(T=2r/9\), or \(2(r-1)/9\) in the single forbidden-node factor.
If actual total cost \(C\ge T\), the unit-charge conclusion follows from
\(W\ge0\). Otherwise each retained nonzero coordinate has inverse diagonal
strictly below T. Because every nonzero contact is at least 1 and all inverse
entries are nonnegative, replacing the contacts in a fixed support by 1 is a
valid cost lower bound, including all same-row cross terms. The independent
checker enumerates **every subset** of that full cheap-coordinate set. There
is no support-size or finite-angle cutoff.

This reconstructs exactly 71 singleton and 16 multiple-contact occurrences;
none of the latter has more than two retained coordinates. For each singleton,
the resistance quadratic is reconstructed independently and minimized at the
endpoints and its stationary point. With \(B=31/30\), the bound
\[
\min(\sigma,\tau)^2\ge
\max\!\left\{\min(h,k)^2,\frac{\max(h,k)^2}{B^2}\right\}
\]
is valid regardless of which endpoint has the larger actual height. Separate
zero/nonzero omitted-coordinate cases give the claimed unit charge in all 71
occurrences. Multiplying a separately valid energy lower bound by a separately
valid nonnegative height lower bound is safe; no unproved correlation is
required.

The sole inexpensive multiple support is the ordinary E7 pattern with one
contact at node 4 of the first row and node 7 of the second. Polarization of
the grounded resistance data gives exactly
\[
E(u,v)=\tfrac12(1+u^2+v^2),\qquad C\ge1.
\]
If both omitted coordinates vanish, \(u=v=1\), so \(W+C\ge5/2\).
Otherwise one height is at least 3; the common ratio bound gives
\(W+C\ge9/(2B^2)+1\). Both exceed \(14/9\).

The other 15 occurrences have cost at least \(3/2\) and possible deficit at
most \(2/9\), just using \(W\ge0\). More precisely, the four ordinary E7
occurrences have cost \(3/2\) and deficit at most \(1/18\); the four ordinary
E8 occurrences have cost \(14/9\) and deficit at most \(2/9\); the seven
forbidden-node E8 records have cost \(3/2\) and deficit at most \(1/18\).
The seven records represent alternative forbidden-node choices, not seven
factors simultaneously present. The review does not assume every possible
defect is attained; these are safe upper bounds on actual failure.

## Small A cycles: complete domains and exact dual certificates

For each cycle size \(m=2,\ldots,17\), a separately computed inverse of the
grounded A-type path reconstructs all point inner products:
\[
K_{ij}=\min(i,j)-ij/m.
\]
These agree with the explicit Euclidean vectors
\(v_k=e_1+\cdots+e_k-(k/m)\sum_{a=1}^m e_a\). Their squared distances are
\(R(i-j)\). Thus the signed Gram used for the two point lists is positive
semidefinite, including singular and coincident-potential cases.

The independent checker lists all retained nodes whose individual cost is at
most 2, enumerates all their subsets, and computes every same-row cross term.
It then rebuilds both endpoint domains, all domain pairs with total minimum
cost below T, and both zero/nonzero possibilities for each omitted coordinate.
The complete list has **580 cases**, agreeing exactly with the primary.
Actual contact coefficients may be arbitrary real numbers at least 1; the
support enumeration therefore covers more than a finite label set.

For every case it verifies the saved rational barycentric vector x and
multipliers \(\lambda\):
\[
x\ge0,\quad Ax=(1,1),\quad
s=Hx+A^T\lambda\ge0,\quad x_i s_i=0.
\]
It also checks \(x^THx=-\sum\lambda_i\). For any other feasible z,
\[
z^THz-x^THx=(z-x)^TH(z-x)+2z^Ts\ge0.
\]
This independently establishes the global convex minimum; no trust in the
producer's active-set optimizer or its treatment of degenerate minimizers is
needed. One verified dual witness suffices for each completely enumerated
support case. Allowing zero barycentric coordinates relaxes the actual
positive-contact domain, which is valid for a lower bound.

The actual heights are at least the two full support sizes. Combining the
height estimate with the verified convex energy and support-cost lower bounds
proves unit charge in every case.

## Unbounded classical ranks

For a large A cycle, a retained coordinate at cyclic distance at least 3 costs
at least \(3(m-3)/m>2\), and is impossible under the individual endpoint
budget. The inverse diagonal increases up to the antipode, so this excludes
all larger distances, not merely distance 3. At distance \(d\le2\), the two
directions have linear losses \(d(m\pm\varepsilon)/m\), where
\(\varepsilon=m\bmod2\), while the cost is \(d(m-d)/m\). Therefore the
loss-to-cost ratio is at most \((m+\varepsilon)/(m-2)\).

The frozen nonnegative-quadratic-remainder identity gives the energy lower
bound after these linear losses. Since both heights are at least 1,
\(W\ge E\); since \(C\le4\),
\[
W+C\ge R(\lfloor m/2\rfloor)-\frac{4(2+\varepsilon)}{m-2}.
\]
The independent checker reconstructs the cleared excess numerators
\[
m^2+6m-304\quad(m\text{ even}),\qquad
m^3+6m^2-457m+18\quad(m\text{ odd}).
\]
Substitution \(m=18+2t\) or \(19+2t\) gives strictly positive constant
terms and nonnegative remaining coefficients, with positive denominators.
This proves every larger rank. It also verifies the cross-multiplied identities
with the point resistance and target, independently of symbolic parsing.
Rotation preserves the ordinary opposite-omission estimate and avoids any
prescribed node when \(m\ge3\). The remaining A1 avoid-node case has target
zero.

For B/C/D the review independently re-derives all 24 prefix-coordinate classes
using rational interval segments and affine polynomials in nonnegative rank
and node-index parameters. It verifies the point-energy surplus, and in each
class either the linear loss is at most its inverse diagonal or that diagonal
strictly exceeds T. The latter excludes the contact in the relevant range
\(C<T\). Modified B2 is checked separately by an exact grounded matrix.
The geometric realization, complete class partition, and omission symmetries
are the same explicitly inherited ones as in the
[independently reviewed four-thirds proof](fresh_nikulin_unit_charge_20261015_review.md).
No new rank cutoff or omitted classical family is introduced here.

## Shared budget and common-height split

Define an actual failing factor by \(T-W-C>0\). The complete local analysis
shows that it consumes actual cost at least \(3/2\), and loses at most
\(2/9\). Since the two endpoint budgets give \(\sum C\le4\), three failing
factors would cost at least \(9/2>4\). Thus at most two factors can fail.
This counts factors, not the number of candidate support records.

With \(N=n-1=\sum r\) and the single avoid-node target penalty \(2/9\),
\[
\sum W\ge\frac{2N}{9}-\frac29-4-\frac49
=\frac{2N-42}{9}.
\]
The unchanged Lorentz identity consequently gives
\(-(x,y)\ge M(N-39)/9\) for \(1\le M\le31/30\). In the other range the
reviewed four-thirds estimate, with the **same preset omissions**, gives
\(-(x,y)\ge M(N-43)/9\). At \(N=165+t\), the respective lower bounds are
at least
\[
14+t/9,\qquad
14+1/135+31t/270.
\]
Both contradict the strict narrow-facet pairing bound \(-(x,y)<14\).
The equality endpoint of the near range is handled by that strict bound.
The ideal-vertex branch therefore gives \(n\le165\); the inherited compact
facet branch gives \(n\le143\).

This step genuinely uses the shared cost budget and common height ratio. It
does not assert a uniform unit-charge lemma, which is false, or assume that
the balanced counterexamples disappear. The same omission choices in both
ranges avoid any circular dependence on the subsequently found walls.

## Dependencies and final status

The new checker and this review establish the changed local inequalities,
complete support coverage, exact dual certificates, all-rank algebra, and
global numerical accounting. The following remain explicit previously reviewed
geometric inputs: Nikulin's narrow facet and its strict constant; the
compact/ideal facet dichotomy; actual terminating walls and framing ridges;
nonzero contact rows with the Coxeter gap; a common height ratio across all
affine factors; the two endpoint Schur budgets; and the compact-facet 143
branch. They are bound through the frozen
[169](fresh_nikulin_unit_charge_20261015_dimension169.md) and
[166](fresh_nikulin_coupled_charge_20261015_dimension166.md) chains.
The present review does not silently replace those geometric premises by a
finite computation.

Independent generation and saved replay pass. Source hashes, syntax, and local
links pass. No primary changes were requested.

```sh
python3 scripts/fresh_nikulin_height_charge_20261015_review.py
python3 scripts/fresh_nikulin_height_charge_20261015_review.py --verify
```
