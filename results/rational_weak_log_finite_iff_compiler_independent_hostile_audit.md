# Independent hostile audit of `RLogComp`

**Audit status (2026-08-11): PASS after one necessary endpoint-finiteness
repair, now present in the owner report.**  The audited theorem is a finite
iff/`EMPTY` criterion for the exact weak-log rational species it defines.  It
is not a universal Property-4 existence theorem, and the generic elimination
engine is not implemented.

The audit pins these owner bytes:

| file | SHA-256 |
|---|---|
| `results/rational_weak_log_finite_iff_compiler.md` | `48cf57f67dea15ed7a307047506e1d23ffd6841af6dcb29798d99e9eb5e34ffe` |
| `scripts/certify_rational_weak_log_finite_iff_compiler.py` | `45de835b6a75f2abc0ff9d1bdd02314279ceee09db90f575ff4655a43569f5ee` |
| owner archive | `87f391c8c82f6ef9739c5b034b53c14f24dc70fba397aa8a52011c05dc31906d` |
| owner test | `f8f0111c2a955ec58f14ce5786a848bfcabba25c4b6c4ee8f65379395f49720c` |

The independent replay is
[`audit_rational_weak_log_finite_iff_compiler.py`](../scripts/audit_rational_weak_log_finite_iff_compiler.py),
with its frozen
[`archive`](../data/enumerations/rational_weak_log_finite_iff_compiler_independent_audit.json)
and
[`tests`](../tests/test_rational_weak_log_finite_iff_compiler_independent_audit.py).
It imports no code from the owner certifier.

## 1. Verdict

The implication audited is exactly

```text
fixed maximal arithmetic P
  -> finite HOST/MMP/PRIME/NEF constructible presentation
  -> nonempty iff P has a literal full-Picard weak-log rational return.
```

The following pressure points pass:

1. the physical-square and Picard-overlattice host list is finite and remains
   in the same rational quadratic space;
2. all rational minimalizations are captured by facet-supported contractions,
   including a boundary curve which becomes a `(-1)`-pivot only later;
3. the ruled endpoints form a finite list after the repaired negative-section
   argument;
4. simultaneous prime divisors and the fixed global-generation condition are
   constructible on finite-type cluster charts;
5. the null branch has an independent effective-nef certificate;
6. Fujino's theorem really gives the necessary fixed system `|192D|` in the
   big branch; and
7. reflection maximality is invoked only after `PRIME` and actual `NEF`, so
   the exact-return step is not circular.

`EMPTY` has the narrow, correct meaning: no carrier in this declared
weak-log rational species.  It says nothing by itself about projected K3,
larger-Picard, Enriques, general-type, or future carrier categories.

## 2. `HOST`: finiteness and the rational quadratic form

For every facet line there are only three physical squares,

\[
                            C_i^2\in\{-1,-2,-4\}.
\]

Thus there are `3^f` physical Grams.  An actual full-Picard return has the
integral full-rank wall lattice $R=\sum\mathbf ZC_i$ inside the unimodular
rational-surface lattice $L$.  Integrality gives

\[
                              R\subseteq L\subseteq R^*.
\]

The discriminant group $R^*/R$ is finite, so its integral unimodular
overlattices are a finite exhaustive list.  Facet span, and hence uniqueness
of the adjunction solution for $K$, follows from finite hyperbolic volume.
The characteristic and $K^2=10-\rho$ tests are necessary.  False positive
numerical lattices cause no problem: their later geometric loci are empty.

This construction is performed inside each candidate physical rational
space.  Passing to an abstract real Lorentz lattice with the same inertia
would be invalid; determinant square class, rational isotropy, and local
Hasse invariants cannot be discarded.

The `(2,4,6)` check is robust under rational root rescaling.  Any rational
symmetrized Gram for the connected diagram can be written

\[
 H=\begin{pmatrix}-a&x&0\\x&-b&y\\0&y&-c\end{pmatrix},
 \qquad ab=2x^2,\qquad 3bc=4y^2.
\]

The rational diagonal change

\[
                 \operatorname {diag}(b/x,1,b/y)
\]

takes it, up to the common factor $b$, to

\[
 \begin{pmatrix}-2&1&0\\1&-1&1\\0&1&-4/3\end{pmatrix}.
\]

Hence the pinned rational form is not an artefact of one normalization.
The mod-three descent in the owner report proves it anisotropic.  The
independent replay also exhausts all `27` physical square assignments and
finds no assignment having integral crossings at both the `4`- and
`6`-edges.  Either calculation rejects a rank-three rational Picard host.

For a literal rational full-nef return all facet lines are rational divisor
lines, so all cyclic Cartan products are rational.  Therefore a non-rational
minimal Vinberg field is a valid immediate rejection gate.

## 3. `MMP`: facet support and the repaired endpoint bound

An irreducible negative curve $C$ on a smooth surface exposes a genuine nef
facet.  If $H$ is ample, then

\[
        N=H+\frac{H.C}{-C^2}C
\]

satisfies $N.C=0$, while $N.\Gamma>0$ for every other irreducible curve
$\Gamma$, because $C.\Gamma\ge0$.  Also $N^2>0$.  Thus on a literal
finite polyhedral full-nef return every negative curve is among the displayed
facets.

It follows that a curve contracted at a later MMP stage already has a
negative strict transform on the original surface and is in the finite facet
list.  The numerical pullback

\[
        e_C=C-E(E^{\mathsf t}GE)^{-1}E^{\mathsf t}GC
\]

is the class orthogonal to the previously contracted exceptional
configuration.  The conditions $e_C^2=K.e_C=-1$ are exactly the numerical
conditions for its current integral image to have arithmetic genus zero and
square `-1`; hence that image is a smooth exceptional curve.  This applies
equally to an original `(-2)`- or boundary `(-4)`-facet whose square has risen
after earlier contractions.  Invalid numerical orders simply have empty
reverse cluster loci.

The original report initially omitted a real finiteness issue.  The marked
pair $(\operatorname {NS}(\mathbb F_e),K_{\mathbb F_e})$ remembers only the
parity of $e$, so it cannot by itself make the endpoint list finite.  The
repaired argument is decisive.  For $e\ge2$, the negative section has on
the original surface an irreducible strict transform with

\[
                     \widetilde\Sigma^2=-e-t,
\]

where $t\ge0$ counts blowups on its successive transforms.  It is a
displayed negative facet, so $e+t\in\{1,2,4\}$.  Since
$\mathbb F_1$ is not minimal, the only endpoints are

\[
                         \mathbb P^2,
   \qquad \mathbb F_e\quad(e=0,2,3,4).
\]

The possible pairs with $e\ge2$ are exactly

\[
                    (e,t)=(2,0),(2,2),(3,1),(4,0).
\]

For each fixed endpoint and fixed finite contraction length, the iterated
universal point blowup is of finite type.  Proper/free/satellite proximity
conditions are locally closed, with only finitely many Enriques diagrams.
No quotient by the generally positive-dimensional automorphism group of
$\mathbb F_e$ is needed: the compiler asks for a finite-type presentation,
not a finite set of isomorphism classes.

## 4. `PRIME`: simultaneous constructibility

On every ordered cluster chart the divisor class of each facet defines a
relative line bundle.  After the standard finite flattening stratification,
its effective divisors are represented by a relative projective section
space, equivalently by the corresponding divisor part of the relative
Hilbert scheme.  The geometrically integral fibre locus is constructible.
On explicit blowup charts the same condition is an exact-multiplicity open
condition together with finitely bounded factorization avoidance.

Taking the fibre product over the one cluster base is essential and is done
in the owner construction: all facet curves occur on the same surface.
Their adjunction rows give arithmetic genus zero, so a geometrically integral
facet is automatically smooth rational.  Distinct prime divisors have
nonnegative local intersection multiplicities; prescribed intersection zero
therefore makes all boundary components disjoint.  The coefficient-one-half
boundary pair is klt.

## 5. `NEF`: the null branch

Let $0\ne D$ be null and lie in the closed finite-volume wall cone.  Its
projective ray is an ideal vertex.  The active normals span $D^\perp$.
Choose a linearly minimal spanning subset.  Its Gram is negative
semidefinite of corank one; after sign change it is a symmetric singular
M-matrix.  If it is reducible, exactly one irreducible block is singular.
That block has a strictly positive kernel vector, and its corresponding
normal sum spans the radical $\mathbf RD$.

The sign is not arbitrary.  A chamber-interior ample class pairs positively
with every effective facet and with the chosen future $D$, so the positive
kernel sum is a positive multiple of $D$.  Rationality and denominator
clearing give

\[
                         F=mD=\sum b_jC_j
\]

with $b_j>0$.  For a support component, $F.C_j=mD.C_j=0$; every other
irreducible curve meets all prime support components nonnegatively.  Thus
$F$, and hence $D$, is nef.  This is independent of the desired full-nef
equality and does not require semiampleness.

The replay independently verifies the affine-$E_8$ calibration: nine
coordinate-independent roots span $D^\perp$, their Gram has rank eight,
and the positive coefficients

\[
                        (2,4,6,5,4,3,2,1,3)
\]

reconstruct $D=(3,-1^9)$.

## 6. `NEF`: Fujino's fixed big-branch gate

The application of Fujino is correct.  On an actual return set

\[
       \Delta=B/2,\qquad L_D=2D=-2K-B.
\]

The divisor $L_D$ is nef Cartier.  Since $D$ is nef with $D^2>0$, it
is big.  The pair $(S,\Delta)$ is klt, hence

\[
                 0\cdot L_D-(K+\Delta)=D
\]

is nef and log big.  Theorem 1.1 and Remark 1.2 of Fujino's
[effective base-point-free theorem](https://www.math.nagoya-u.ac.jp/~fujino/effective4.pdf)
allow $a=0$ and give

\[
 m(2,0)=2^3\,3!\,(0+2)=96.
\]

Therefore every actual return has

\[
                         |96L_D|=|192D|
\]

base-point-free.  Conversely, global generation of this Cartier multiple
makes $D$ nef.  The gate is therefore necessary and sufficient; `192` is
not asserted to be optimal.

After stratifying the cluster base so that $H^0$ commutes with base change,
the evaluation cokernel

\[
 \pi^*\pi_*\mathcal O(192D)\longrightarrow\mathcal O(192D)
\]

has closed support.  Its image in the stratum is closed because the surface
family is proper.  The complement is the fibrewise globally generated locus,
so the big `NEF` gate is constructible.

## 7. Decidability and exact return

The endpoint list, cluster charts, relative section/Hilbert spaces,
integrality conditions, factorization strata, and evaluation maps are all
finite-presentation algebraic data over a computable number field.  After a
finite constructible stratification, nonemptiness is a first-order question
over an algebraically closed field.  Gröbner elimination and the
Rabinowitsch trick terminate.  This proves an abstract decision procedure;
the owner correctly states that the enormous generic degree-`192D`
elimination is not implemented by the current PB-Atlas code.

Once a point passes `HOST`, `PRIME`, and an independent `NEF` gate, the
previously audited weak-log exact-return theorem applies.  Adjunction and
nefness classify every hidden negative curve as a boundary `(-4)`, a
boundary-orthogonal `(-2)`, or a `(-1)`; every corresponding reflection is
integral on the adapted host.  Ambient arithmetic reflection maximality
saturates those mirrors, and an ample class fixes their orientation.  This
gives the literal ideal equality

\[
                        P=\mathbf P(\operatorname {Nef}S).
\]

Maximality is not used to manufacture sections or to prove $D$ nef.

## 8. Exact limits

The audit does not certify any of the following stronger statements:

- every raw maximal arithmetic chamber has a weak-log rational carrier;
- `EMPTY` excludes a carrier in another geometric species;
- the current three-branch Property-4 disjunction is exhaustive;
- the current repository implements the generic elimination; or
- real inertia can replace the rational HOST tests.

Consequently the dimension bound $n\le25$, the field conclusion
$k(P)=\mathbf Q$, and carrier-relative finiteness apply only after
`RLogComp` returns a nonempty weak-log carrier.  The negative triangle
examples illustrate possible `EMPTY` outputs; they are not used to prove the
iff compiler itself.

## Replay

```bash
python3 scripts/audit_rational_weak_log_finite_iff_compiler.py \
  --check-archive

python3 -m pytest -q \
  tests/test_rational_weak_log_finite_iff_compiler_independent_audit.py
```
