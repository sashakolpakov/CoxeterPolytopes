# Full source review: A18 isolated-cap domain and norm screen

21 October2026. **PASS before execution.** This is a complete mathematical
specification and worker-source review; no worker import or execution,
enumeration, scientific process, Lean or TeX job was performed by this
reviewer. No output counts or runtime PASS are asserted here.

| Reviewed input | SHA256 |
| --- | --- |
| scripts/check_a18_isolated_cap_domain_and_norm_screen_20261021.py | `29da7434634c80fb001784d9d5f0d69c3d6247a7b6de1e55e2a0eb6b7ff22932` |
| results/arithmetic_whole_a18_isolated_cap_domain_and_norm_screen_20261021.md | `06172cc3d1f1e218a29d7e20916f105366d10344f4e086a4aeb4d4eed45e75ac` |
| results/arithmetic_whole_a18_point52_actual_field_degree_reduction_20261021.md | `b57aa75adbde5c1a1b6f78354c58bba47d03b2047a923438b133be9aee5f2707` |
| manuscripts/compact_arithmetic_dimension_bound.tex | `1f09df2976085fd1ad1aa56ed35f1ab182ff95a8332dcfa0b0b80d4849b82203` |

The reviewer authored the inherited b57 field/degree source with root's
shared quartic-discriminant contribution and derived a separate pending
quadratic inventory. This is not independent review of b57. The new root
specification, polynomial enumeration logic, raw-time formulas, norm
screen and reflected-pair deduction were independently checked here.

## Mathematical coverage and exact domains

The fresh rank18 source supplies the actual isolated-dot normalization,
4<beta<40/3, all other conjugates in(0,19/35), actual K=Q(beta), and
degree<=3. The reviewed clarification explicitly uses rational Meyer
isotropy in19 variables to exclude degree1. Thus selecting degrees2,3
is justified; rational-field exclusion was not silently attributed to b57.

The complete boxes are respectively9*7=63 quadratic pairs and10*14*3=420
cubic triples. These are coefficient domains, not lists of realizable
geometric caps. `coefficients` visits each tuple once with the stated
ascending polynomial signs. Quadratic irreducibility is tested by a
positive nonsquare discriminant. Cubic irreducibility uses every signed
integer divisor of M, complete by the rational-root theorem for a monic
cubic. Negative roots are tested too; no reducible cubic is retained.

`trim`, Horner evaluation and Euclidean remainder use `Fraction` only.
The derivative sequence j*p[j], j1..d, has exactly the correct ascending
derivative coefficients. Negative remainders strictly decrease degree;
irreducibility in characteristic zero ensures the final nonzero term is
constant. Sturm variation omits zero evaluated entries and uses V(a)-V(b).
All endpoints are rational and retained degree>=2 polynomials irreducible,
so p never vanishes at an endpoint. The two disjoint root-count tests
retain exactly d-1 distinct roots in(0,19/35) and one in(4,40/3).

The ten index tests i5..14 then require every small root in(0,38/k),
k=i(19-i). This is precisely the strict full-Gram conjugate Schur
condition. No physical weak cap in the stated isolated domain is omitted.
No other quadratic helper/finite role is represented by this enumeration.

## Raw-time normalization and both independent norm calculations

For raw old square-two roots and an actual isolated cap, scaling by
-beta/a is over K and gives c^2=2beta, marked product-beta. Its projection
is beta*Omega_i, Omega_i^2=k/19. Thus F=beta(k beta-38)/19 exactly.
The physical F is positive and every other F conjugate negative.

Independently expanding its norm gives

    quadratic: M(Mk^2-38Tk+1444)/19^2;
    cubic: M(Mk^3-38Sk^2+1444Tk-54872)/19^3.

`norm_numerator` implements those formulas with the correct signs.
Separately `multiplication_norm` reduces F modulo p, computes all
columns beta^j F in the basis1,beta,... and takes the rational determinant.
The matrix orientation and Gaussian elimination, including row-swap
signs, are correct. Every retained cell requires exact equality of these
two computations, nonzero numerator, and the required negative/positive
sign in degree2/3. This is not a floating-point root/norm approximation.

## Necessary grouping and the reflected-pair lemma

An actual common field forces equal polynomial-discriminant square class
by the integral discriminant-index formula. Common K-valued time lines
force F1/F2 to be a K-square, so their rational norm ratio is a positive
rational square. The fixed19^d denominator cancels. Grouping by degree,
discriminant squarefree part and absolute norm-numerator squarefree part
therefore cannot split an actual pair.

The trial-division invariant proves the returned part squarefree and
the recorded identity n=part*square^2. The grouping is deliberately only
necessary: equal discriminant classes need not identify cubic fields,
and equal norm classes do not prove a K-square. The output never claims
otherwise. Full cells, all polynomial/k types and multi-type groups are
retained; the cubic singleton predicate is computed, not hardcoded.

For a same-polynomial reflected pair, the physical beta is the unique
root above4, hence equal at both indices. The positive time coefficients
are equal. Direct raw inverse products give
c_i.c_j=2beta-beta^2*i(19-2i)/19. In the actual Cartan cycle the two
isolated contacts and the new norms cancel to give plus or minus this
raw product, so it is integral. Its nonintegral remainder would have
norm denominator19^d and numerator M^2*|i(19-2i)|^d. For i5..9 the
factors45,42,35,24,9 and M<=7 are prime to19. This is a valid integral
norm contradiction. The specification's CONDITIONAL at-most-one cubic
index consequence follows only if every computed cubic group has one
polynomial and one k. It supplies no omitted-endpoint energy argument.

## Termination, bindings, storage and honest output boundary

The source uses only the standard library, finite outer loops, rational
polynomials of degree<=3, at most3-by-3 matrices and terminating integer
trial division. There is no eager corpus loader, scientific import or
unbounded root/Pell search. Before/after hashes bind domain source,
specification, rational obstruction source and worker; the first and
third have fixed expected source pins. The reviewed specification/worker
hashes above must also be bound by the external launch/receipt.

The worker requires the bounded-owner environment and a NEW direct
results destination, then exclusively creates that output. An owner
environment alone is not a resource admission: the unchanged guarded
runner supplies the actual lock, monitor and caps. No source-review
verdict waives its fresh40GiB/6GiB admission gate or59s/1GiB limit.

Count identities cover all coefficient decisions and all index decisions.
The final specification correctly says it records COUNTS for all visited
rows, exact rejection counts, complete retained polynomials/cells/groups
and squarefree witnesses; it does not promise a separate rejected-row
table. This documentary precision and the explicit rational-degree input
were corrected before this final binding. No remaining source blocker was
found. Actual execution and any later geometric/energy deduction require
their own separately source-bound evidence.
