# The Type-I rank-ten anti-Gram residue is killed by its forced glue root

**Status:** exact obstruction with a replayable lattice certificate.

## Result

The unique negative-semidefinite rank-ten residue in the exhaustive Type-I
downstairs anti-Gram census cannot be a good Enriques lift package.

Its abstract generated quotient initially looks viable:

\[
 A_0\cong E_8(-2)\oplus\langle-4\rangle^2.
\]

It is negative definite, has determinant \(2^{12}\), discriminant group

\[
 A_{A_0}\cong (\mathbf Z/2)^8\oplus(\mathbf Z/4)^2,
\]

and every norm is divisible by four, so `A_0` itself is rootless.  The
obstruction is not its Gram matrix or its ordinary primitive-embedding
genus.  It is the simultaneous half-sum glue.

The final two Type-I downstairs roots obey

\[
 r_{11}-r_{12}\in2L_{10},
\]

whereas the residue has

\[
 \delta_{11}^2=\delta_{12}^2=-4,
 \qquad \delta_{11}\cdot\delta_{12}=0.
\]

For any choices of sheets, absorb the signs into the anti-vectors and write

\[
 b_i=\frac{u_i+\delta_i}{2}\in\Lambda_{K3},
 \qquad u_i=(0,r_i,r_i).
\]

Then

\[
 \alpha=
 b_{11}-b_{12}-\frac{u_{11}-u_{12}}2
 =\frac{\delta_{11}-\delta_{12}}2
 \in\Lambda_{K3}^{-\tau}=N.
\]

It satisfies

\[
 \alpha^2=\frac{-4-4}{4}=-2.
\]

Changing either sheet replaces the difference by the corresponding signed
sum and gives the same square.  Thus **every** compatible half-sum gluing
forces an anti-invariant `(-2)` vector.

Consequently the saturated anti-lattice contains the index-two overclosure

\[
 A_1=
 E_8(-2)\oplus
 \left(\langle-4\rangle^2+
       \mathbf Z\frac{\delta_{11}-\delta_{12}}2\right)
 \cong E_8(-2)\oplus A_1(-1)^2.
\]

It is rootful.  The `L2` period-rootlessness gate fails, every invariant
class is orthogonal to these forced roots, and hence the `L3` invariant
root-regular ample gate also fails.  `L4` is not reached and cannot repair
the failure.

Together with the exhaustive 18-orbit census and the already certified
rank-nine positive package, this eliminates the only competing
negative-semidefinite orbit.  Therefore, among all Type-I anti-Grams allowed
by the downstairs two-sheet intersection alphabet,

```text
good Enriques lift package
    => the unique rank-nine E8(-2)+<-4> orbit,
```

up to independent sheet switches.  This closes the anti-Gram projection
stage for Type I; it does not claim a new surface or a new abstract Coxeter
polytope.

## Raw quotient lattice

The twelve-vector residue has inertia `(0,10,2)`.  A saturated integral
basis of its radical is given by the columns

\[
\begin{pmatrix}
-2&-2\\
 3& 2\\
-4&-2\\
 3& 1\\
-2& 0\\
-1& 1\\
 0&-2\\
 1& 2\\
 2& 1\\
 0& 1\\
 0& 0\\
 0& 0
\end{pmatrix}.
\]

Its maximal minors have gcd one.  Eliminating generators 2 and 9 gives an
integral basis for the quotient.  The final two basis vectors are orthogonal
`<-4>` summands.  Dividing minus one half of the first rank-eight block
produces a positive-definite even unimodular lattice of rank eight, hence
`E8`.  The quotient Smith diagonal is

```text
2^8, 4, 4.
```

Thus

\[
q_{A_0}=q_{E_8(-2)}\oplus\langle-1/4\rangle^2.
\]

Adjoining the forced root changes the final two summands to two orthogonal
`(-2)` roots.  The resulting Smith diagonal is `2^10`, determinant
`2^10`, and discriminant form

\[
q_{A_1}=q_{E_8(-2)}\oplus\langle-1/2\rangle^2.
\]

The certificate also gives a primitive embedding of the last
`A1(-1)^2` into `U+U(2)`.  In the standard basis with Gram
`U+U(2)`, its two columns are

\[
(1,-1,0,0)^t,
\qquad (1,1,1,-1)^t.
\]

They are orthogonal roots and have a unit maximal minor.  Adding the
`E8(-2)` identity factor embeds `A1` primitively in

\[
N=U\oplus U(2)\oplus E_8(-2).
\]

Thus primitive embeddability is not the obstruction; rootlessness is.

## The rank-twenty K3 edge case

Adjoin all twelve half-sums to

\[
M\oplus A_0,
\qquad M=L_{10}(2).
\]

Their combined mod-two rank is ten, while their projection to the invariant
discriminant group has rank nine.  The one-dimensional kernel is exactly the
forced anti-saturation above.  The resulting even Picard lattice `P` has

```text
signature:       (1,19)
determinant:     -4
Smith diagonal:  1^18,2,2
nonzero q-values modulo 2Z:  3/2,3/2,1.
```

Hence

\[
q_P\cong\langle-1/2\rangle^2
       =-q_{A_1\oplus A_1}.
\]

Gluing to the positive lattice `A1+A1` gives the K3 lattice.  In other
words, an ordinary singular-K3 period is not obstructed.  But the entire
rank-two positive period plane is orthogonal to the anti-algebraic lattice,
including its forced roots.  It therefore lies on the forbidden Enriques
root hyperplanes.  The corresponding involution cannot be fixed-point-free.

This illustrates why rank-ten inertia alone was insufficient: the raw
generated lattice was rootless, but the **compatible primitive saturation**
was not.

## Replay

- [exact certifier](../scripts/certify_enriques_type_i_rank_ten_residue_obstruction.py)
- [JSON certificate](../data/enumerations/enriques_type_i_rank_ten_residue_obstruction.json)
- [focused tests](../tests/test_enriques_type_i_rank_ten_residue_obstruction.py)
- [upstream exhaustive anti-Gram census](enriques_type_i_downstairs_anti_gram_census.md)
- [positive rank-nine package](enriques_type_i_explicit_lift_package.md)

Run

```bash
PYTHONPATH=src:. python3 \
  scripts/certify_enriques_type_i_rank_ten_residue_obstruction.py --check-archive

PYTHONPATH=src:. pytest -q \
  tests/test_enriques_type_i_rank_ten_residue_obstruction.py
```
