# Finite congruence monodromy for Coxeter--Borcherds dessins

**Status (2026-08-09): proved; exact E10 calibration and independent hostile
audit passed.** This theorem makes the finite-dessin layer rigorous. It
preserves all marked finite Coxeter labels and separates any prescribed
finite list of matrix nonrelations. A separate bounded-height lemma makes one
finite quadratic dessin recover the exact Gram matrix. No finite layer
detects an infinite, parallel, or ultraparallel pair by itself.

## 1. Exact finite-order level

Let \(A=\mathcal O_k\), let \(L=A^r\) carry an integral quadratic form, and
let

\[
 s_1,\ldots,s_f\in O(L)
\]

be distinct nonidentity integral reflection matrices. Suppose the exact
characteristic-zero order of \(s_i s_j\) is \(m_{ij}<\infty\) for every pair
in a marked finite set \(E_{\mathrm{fin}}\).

Form the following finite collection of nonzero matrices:

\[
 s_i-I,\qquad s_i-s_j,\qquad
 (s_i s_j)^{m_{ij}/q}-I
 \quad(q\text{ prime},\ q\mid m_{ij}).                 \tag{1}
\]

Choose one nonzero entry from each matrix in (1), multiply them, and choose a
prime ideal \(\mathfrak p\) which contains none of those entries. This is a
terminating computation: exact matrix arithmetic proves that the entries are
nonzero, and one enumerates prime ideals away from their finite prime
support. One may simultaneously avoid \(2\), the lattice discriminant, and
all \(m_{ij}\), giving a good nondegenerate tame level.

> **Finite-order theorem.** Reduction modulo \(\mathfrak p\) keeps the
> \(s_i\) as distinct nonidentity involutions and gives
>
> \[
>   \operatorname{ord}(\bar s_i\bar s_j)=m_{ij}
>   \quad ((i,j)\in E_{\mathrm{fin}}).                 \tag{2}
> \]

Every reduced generator still squares to one, while its chosen nonidentity
witness survives. The order \(h\) of a reduced product divides
\(m=m_{ij}\). If \(h<m\), choose a prime \(q\mid m/h\). Then
\(h\mid m/q\), contradicting the surviving last witness in (1). Two
distinct involutions whose product has order \(m\) generate the dihedral
group of order \(2m\). Thus the finite quotient retains the complete marked
finite dihedral data, not just upper bounds on the orders.

This proof does not need the good/tame exclusions; they are added for the
clean orthogonal-cover interpretation.

## 2. Congruence towers separate finite nonrelation lists

Fix any nonzero prime ideal \(\mathfrak p\). For a nonidentity matrix word
\(w\), let

\[
 \nu_{\mathfrak p}(w-I)=
 \min_{\substack{a\ne0\\a\text{ an entry of }w-I}}
 v_{\mathfrak p}(a).
\]

For a prescribed finite set \(F\) of nonidentity matrix words, the level

\[
 a(F)=1+\max_{w\in F}\nu_{\mathfrak p}(w-I)             \tag{3}
\]

keeps every \(w\in F\) nonidentity modulo
\(\mathfrak p^{a(F)}\). The same argument separates finitely many prescribed
inequalities \(w\ne w'\) by applying it to \(w(w')^{-1}\). Since
\(\bigcap_a\mathfrak p^a=0\), the complete tower is faithful on the given
matrix group.

The matrix qualifier matters: this does not prove that a separately
presented abstract group acts faithfully. In the geometric reflection
setting the group is already the matrix subgroup of \(O(L)\).

## 3. The height-to-one-level reconstruction lemma

The preceding level preserves finite orders but need not recover
coefficients. The following is a separate, stronger use of congruence.

> **Bounded-height injection lemma.** Let \([k:\mathbb Q]=d\), and let
>
> \[
> A_B=\{a\in\mathcal O_k:|\sigma(a)|\le B
>       \text{ for every embedding }\sigma\}.
> \]
>
> If a prime ideal satisfies
>
> \[
>       N\mathfrak p>(2B)^d,                            \tag{4}
> \]
>
> then reduction \(A_B\longrightarrow\mathcal O_k/\mathfrak p\) is
> injective.

For distinct \(a,b\in A_B\),

\[
 1\le |N_{k/\mathbb Q}(a-b)|\le(2B)^d.                 \tag{5}
\]

If \(a-b\in\mathfrak p\), then \(N\mathfrak p\) divides the norm of the
principal ideal \((a-b)\), contradicting (4). The strict inequality is
necessary: for \(k=\mathbb Q\), \(B=1\), and \(\mathfrak p=(2)\), the
entries \(-1\) and \(1\) collide at equality.

Consequently, fix the field, order, prime, and a class-uniform bound \(B\).
A marked finite quadratic dessin

\[
 (L/\mathfrak pL,\bar q,\{\bar r_i\},\{\bar s_i\},
   \text{inward-root marking})                          \tag{6}
\]

recovers the exact Gram matrix: enumerate the finite Minkowski alphabet
\(A_B\), compute each reduced pairing, and take its unique lift. For a fixed
known denominator, first clear it and choose \(\mathfrak p\) away from that
denominator. Reflection matrices without the quadratic/root marking are
insufficient because \(r\) and \(-r\) give the same reflection. An unmarked
finite group also does not recover an unknown number field.

This is the genuine classification leverage: a class-uniform coefficient
height bound converts the entire congruence tower into one exact finite
dessin. The missing theorem is precisely such a bound for a broad
surface-free Coxeter class.

## 4. Algebraic orthogonal-cover interpretation

Let \(\mathcal D\) be the determinant-line type-IV domain and let
\(\Gamma\subset O(\widehat L)\) contain the marked reflections, extended by
the identity on the auxiliary line. Put

\[
 \Gamma(\mathfrak n)=
 \ker\!\left(
   \Gamma\longrightarrow O(\widehat L/\mathfrak n\widehat L)
 \right).
\]

The principal-level map

\[
 [\mathcal D/\Gamma(\mathfrak n)]
       \longrightarrow[\mathcal D/\Gamma]               \tag{7}
\]

is a finite étale map of Deligne--Mumford stacks in characteristic zero. On
the original stack, a mirror has reflection inertia \(C_2\), and a finite
mirror crossing has dihedral inertia of order \(2m_{ij}\). At a level from
Section 1, the images of the reflections in the finite deck group retain
exactly those dihedral groups.

After passage to a sufficiently deep neat principal sublevel, the source of
(7) is an algebraic variety and has no stack inertia. Therefore one must not
claim that inertia lives on the neat cover: the labels survive there as
marked deck monodromy together with the retained root/quadratic local system.
A neat cover has no inertia; that information belongs to the original stack.
The induced map of coarse spaces may ramify along the mirror images even
though (7) is étale as a stack map.

This congruence level is distinct from:

- scaling the lattice to separate Heegner divisor sheets;
- refining level to make a subgroup neat; and
- choosing a large-norm prime for the height lemma.

They can be combined by passing to a common deeper level, but none implies
the others automatically. A Borcherds product pulls back to such a cover;
its existence and divisor support still require the independent
vector-valued cusp-obstruction calculation.

## 5. Exact E10 calibration

The frozen \(E_{10}=U+E_8(-1)\) chamber has dimension \(9\), ten facets, and
a unimodular square-minus-two wall Gram matrix. Its ten integral reflection
matrices have 36 pair products of order two and nine of order three.

At \(\mathfrak p=(5)\):

- all ten reduced involutions are nonidentity and distinct;
- all 45 finite pair orders remain exact;
- the first 24 prescribed powers of \(s_0s_1\cdots s_9\) remain
  nonidentity;
- every obstruction-matrix content is one; and
- \(B=2\), \(d=1\), so \(5>(2B)^d=4\), and the exact Gram is the unique
  lift of its marked mod-five quadratic dessin.

The independent hostile replay tests five traps: the involution \(-1\)
collapsing modulo two; a word collapsing modulo five but reappearing modulo
25; an infinite unipotent acquiring order five modulo five; failure at
equality in (4); and the \(r/-r\) ambiguity. All tests pass with the stated
scope.

The unipotent example is decisive: every element of a finite quotient has
finite order. Hence no single congruence dessin certifies \(m=\infty\) or
distinguishes a parallel pair from an ultraparallel pair. Those relations
remain encoded by the exact characteristic-zero Gram and real chamber. The
full tower separates any finite list of their nonpower relations, but a
single layer never replaces that real data.

## 6. Frozen replay

- [Primary certificate](../scripts/certify_coxeter_congruence_dessin.py)
- [Primary archive](../data/enumerations/coxeter_congruence_dessin_e10_h9.json)
- [Independent hostile audit](../scripts/independent_coxeter_congruence_dessin_hostile_audit.py)
- [Hostile archive](../data/enumerations/coxeter_congruence_dessin_hostile_audit.json)

~~~bash
python3 scripts/run_with_rss_limit.py --memory-gib 0.5 --timeout 120 -- \
  python3 scripts/certify_coxeter_congruence_dessin.py \
  --check-output data/enumerations/coxeter_congruence_dessin_e10_h9.json
python3 scripts/run_with_rss_limit.py --memory-gib 0.5 --timeout 120 -- \
  python3 scripts/independent_coxeter_congruence_dessin_hostile_audit.py \
  --check-output data/enumerations/coxeter_congruence_dessin_hostile_audit.json
python3 -m pytest -q tests/test_coxeter_congruence_dessin.py
~~~

Peak RSS in both exact replays is below 0.025 GiB.
