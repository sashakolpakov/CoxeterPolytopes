# The exact K3-lift boundary for the `A16-DUAL-C2-H11-29` chamber

## Theorem

Let (N) be Alexeev's even 2-elementary lattice `(16,6,1)`, and let
(a\in O(N)) be the transposition whose fixed reflection-chamber section is
the certified finite-volume noncompact `H11/29` Coxeter polytope.  There is
no finite-order lift with `NS(X)=N` on a projective K3 surface.

This is a generic-Picard theorem, not a global special-period no-go.  The
certificate also constructs an explicit low-order special-period lift: an
equivariant order-four Hodge isometry of the K3 lattice with Picard rank 18
and with fixed Néron--Severi lattice exactly (N^a).  That lift is not
geometric, because its Picard coinvariant contains two negated
`(-2)`-classes.  Consequently it cannot preserve a Kähler chamber.

Thus the `H11/29` chamber remains a genuine abstract/divisor-lattice
positive example, but this particular transposition has not been promoted
to a surface automorphism.

## The six-bit obstruction

The action of (a) on the discriminant group

\[
A_N\cong (\mathbb Z/2)^6
\]

is a nontrivial involution satisfying

\[
\operatorname{rank}_{\mathbb F_2}(\bar a-1)=2.
\]

The K3-lattice complement is the unique 2-elementary type of signature
`(2,4)`, length six and parity one; we use

\[
T=U(2)\oplus U(2)\oplus\langle-2\rangle^2.
\]

Because every entry of its Gram matrix is even and its determinant is
(2^6), the natural half-vector map identifies (A_T) with (T/2T).
An equivariant K3 glue therefore requires the reduction of the
transcendental action $B$ modulo two to be conjugate to $\bar a$.  In
particular,

\[
\operatorname{rank}_{\mathbb F_2}(B-1)=2.
\]

Suppose now that `NS(X)=N` and (B) has finite order.  A rational
cyclotomic primary summand of (T_\mathbb Q) not containing the eigenperiod
is perpendicular to that period and would enlarge `NS(X)`.  Hence all six
dimensions of (T_\mathbb Q) must have one cyclotomic type.  The elementary
inverse-totient list leaves only

```text
m = 1,2,3,4,6,7,9,14,18.
```

The cases other than (m=4) are eliminated immediately: (m=1,2) act
trivially modulo two, while the characteristic polynomials for
(m=3,6,7,9,14,18) do not reduce to ((x+1)^6).

For (m=4), the characteristic polynomial is

\[
(x^2+1)^3
\]

and

\[
|\det(B-1)|=\Phi_4(1)^3=8.
\]

But rank two modulo two means that the Smith form of (B-1) has four even
invariant factors, forcing (16\mid\det(B-1)).  This contradiction closes
the entire finite-order denominator when the Picard lattice is exactly
(N).

This gives a reusable criterion: before searching periods or Kähler cones,
compare the required discriminant action with the reductions of the only
cyclotomic types allowed by the transcendental rank.  The Smith valuation
of (B-1) can eliminate the one characteristic polynomial whose reduction
looks correct.

## The exact special-period near miss

The generic obstruction does not prevent algebraic specialization.  On the
displayed (T), take

\[
B=\begin{pmatrix}
0&0&-1&0&0&0\\
0&0&0&-1&0&0\\
1&0&0&0&0&0\\
0&1&0&0&0&0\\
0&0&0&0&-1&0\\
0&0&0&0&0&-1
\end{pmatrix}.
\]

It has order four and characteristic polynomial
((x+1)^2(x^2+1)^2).  In the certificate's exact half-vector bases, the
matrix

\[
P=\begin{pmatrix}
1&0&1&0&1&1\\
0&1&0&1&1&1\\
0&1&1&1&1&1\\
0&0&0&0&1&1\\
0&1&1&0&0&1\\
0&0&0&1&1&1
\end{pmatrix}
\]

intertwines the two discriminant actions.  All 64 classes are replayed with
the full parity-one quadratic value (4q\pmod 8), and (P) is an
anti-isometry.  Graph gluing produces an even unimodular lattice of
signature `(3,19)` and an integral order-four isometry with characteristic
polynomial

\[
(x-1)^{12}(x+1)^6(x^2+1)^2.
\]

Choose the positive eigenperiod

\[
\omega=(1,\sqrt2,i,i\sqrt2,0,0),\qquad B\omega=i\omega.
\]

Its rational orthogonal in (T) is exactly the last two
(\langle-2\rangle) directions.  The resulting Picard lattice has rank 18,
signature `(1,17)`, determinant `-16`, and Smith diagonal
`1^14 2^4`.  Its fixed lattice is integrally equal to the archived rank-12
fixed lattice (N^a), so the fixed hyperbolic space would still have
dimension 11.

The Kähler gate nevertheless fails exactly.  The two added algebraic
vectors have Gram matrix

```text
-2  0
 0 -2
```

and the glued action negates both.  Every invariant positive class is
orthogonal to them.  A finite automorphism would have an invariant Kähler
class obtained by averaging, which is impossible on a root hyperplane;
equivalently, Riemann--Roch makes one sign of each root effective, and an
automorphism cannot send it to its negative.

## What is closed and what remains open

Closed exactly:

- no finite-order realization with generic Picard lattice
  `NS(X)=N=(16,6,1)`;
- a complete 64-class discriminant glue for the smallest order-four
  specialization;
- exact preservation of the archived rank-12 fixed `NS` in that
  specialization; and
- exact failure of its Kähler gate by two explicit roots.

This does not prove that every possible special-period, higher-order lift
fails.  Such a lift would need extra algebraic directions forming a
rootless coinvariant while carrying one of the remaining mixed cyclotomic
types.  The theorem above is the search boundary for that narrower problem.

The source 118-wall chamber is a complete mixed reflection chamber, not the
literal nef cone of the hypothetical special K3 surface.  No nef-cone or
quotient-surface claim is made.

## Hyperbolic handoff and reproduction

The separate Vinberg certificate remains unchanged: the `H11/29` chamber
has 29 facets, 28,798 elliptic edges, 3,876 finite vertices, 253 cusps, and
zero bad edge or affine-neighbourhood audits.  It is finite-volume and
noncompact.

Permanent files:

- [exact obstruction certificate](../data/enumerations/a16_dual_c2_h11_29_k3_lift_obstruction.json)
- [certificate generator](../scripts/certify_a16_dual_c2_h11_29_k3_lift_obstruction.py)
- [regression tests](../tests/test_a16_dual_c2_h11_29_k3_lift_obstruction.py)
- [abstract `H11/29` positive theorem](alexeev_118_dual_h11_29_positive_example.md)

Run from the repository root:

```bash
PYTHONPATH=.:src python3 scripts/certify_a16_dual_c2_h11_29_k3_lift_obstruction.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_a16_dual_c2_h11_29_k3_lift_obstruction.py
```

Frozen payload SHA-256:

```text
96d621be828afb0dcb29040fb85ae024ddde7dc807a8bec48d203a6ae3bb62eb
```
