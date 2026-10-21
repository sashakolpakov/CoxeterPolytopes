# Constant audit and effective label reduction for the AGCat Coxeter theorem

**Status (2026-08-12): proved numerical corollaries, with exact certificate.**
This note audits the four conspicuous constants in
[`agcat_global_maximal_arithmetic_coxeter_finiteness.md`](agcat_global_maximal_arithmetic_coxeter_finiteness.md).
It improves the two label bounds using the existing norm arguments, without
changing the geometric or arithmetic hypotheses.  It does not claim that a
pair surviving a norm inequality is realizable.

## 1. Verdict at a glance

| constant | role | provenance | status for computation |
|---|---|---|---|
| \(n\le995\) | makes the rank set finite | imported Khovanskii--Prokhorov finite-volume Coxeter geometry (with Vinberg in the compact branch) | very crude for arithmetic inputs; not improved by the AG field argument |
| \(A=62\) | bounds a connected full-rank wall frame | imported Nikulin geometric appendix lemma | an unoptimized existence constant; sufficient for finiteness but far too large for naive frame enumeration |
| minimality \(14\) | bounds the closed star of a narrow facet | Nikulin-style narrow-facet geometry, reproved by the cone-plane inequality | sharp for that local inequality; it cannot be lowered by tightening the same estimate |
| fork labels \(m,n<10^6\) | makes the compact-\(H^3\) fork bases finite | deliberately crude project norm envelope | **replaced by \(\max(m,n)\le420\)**; only 525 unordered pairs survive the necessary norm filter |
| tail label \(d<20000\) | makes the transition-tail bases finite | deliberately crude project norm envelope | **replaced by \(d\le126\)**; only 66 labels survive the necessary norm filter |

The last two reductions are internal improvements.  They use the same
integral arrangement coordinates and cyclotomic norms as the existing proof;
no ABSW, Nikulin field-finiteness conclusion, or Coxeter classification is
inserted.

## 2. Exact reduction of the general-fork labels

Retain the notation of
[`compact_h3_general_fork_arrangement_field_finiteness.md`](compact_h3_general_fork_arrangement_field_finiteness.md):

\[
 s_r=4\sin^2(\pi/r),\qquad
 E=\mathbf Q(s_m,s_n),\qquad D=[E:\mathbf Q],
\]

and

\[
 \gamma(r)=
 \begin{cases}
 p,&r=p^a,\\
 1,&\text{otherwise}.
 \end{cases}
\]

The already proved norm inequality is

\[
 1< {1024\over s_{m,0}s_{n,0}}
 \left(
 {\gamma(m)^{D/d_m}\gamma(n)^{D/d_n}\over4^D}
 \right)^e,
 \qquad d_r={\varphi(r)\over2},\quad e=[k:E].       \tag{2.1}
\]

Set

\[
 R_{m,n}=
 {\gamma(m)^{D/d_m}\gamma(n)^{D/d_n}\over4^D}.          \tag{2.1a}
\]

This ratio need not be below one for small pairs (for example
\((m,n)=(3,3)\)). If \(R_{m,n}<1\), then \(e\ge1\) gives
\(R_{m,n}^e\le R_{m,n}\); combining (2.1) with
\(s_{r,0}\ge16/r^2\) gives

\[
 4m^2n^2\,
 \gamma(m)^{D/d_m}\gamma(n)^{D/d_n}>4^D.             \tag{2.2}
\]

If \(R_{m,n}\ge1\), inequality (2.2) is automatic because
\(4m^2n^2R_{m,n}>1\). Thus (2.2) is a necessary filter in both cases:
it follows from the norm bound in the contracting branch and holds
tautologically in the noncontracting branch.

### 2.1 First finite reduction: every label is at most 840

Let \(M=\max(m,n)\). Write
\(h(r)=\gamma(r)^{2/\varphi(r)}\). Since \(h(r)\le3\),
\(D\ge\varphi(M)/2\), and \(3h(M)/4<1\) throughout the range used below,
the actual ratio satisfies
\[
 R_{m,n}\le\left({3h(M)\over4}\right)^D<1.
\]
Consequently the contracting branch of the preceding case split applies.
Using \(D\ge\varphi(M)/2\) gives the weaker necessary inequality

\[
 4^{\varphi(M)/2}
 <4M^4\gamma(M)3^{\varphi(M)/2}.                       \tag{2.3}
\]

The exact integer certificate checks that the reverse strict inequality
holds for every \(841\le M<10^6\). It normally uses the elementary exact
lower bound

\[
 (4/3)^{10}>2^4
\]

and falls back to direct integer powers if needed.  The frozen analytic
argument in the compact-\(H^3\) theorem excludes every \(M\ge10^6\).
Consequently

\[
                         M\le840.                       \tag{2.4}
\]

### 2.2 The exact compositum degree

Put \(\ell=\operatorname {lcm}(m,n)\). Because
\(\mathbf Q(s_r)=\mathbf Q(\zeta_r+\zeta_r^{-1})\), both fields lie in the
maximal real subfield of \(\mathbf Q(\zeta_\ell)\). An automorphism
\(a\in(\mathbf Z/\ell)^\times\) fixes both real subfields exactly when

\[
                a\equiv\epsilon_m\pmod m,qquad
                a\equiv\epsilon_n\pmod n,qquad
                \epsilon_m,\epsilon_n\in\{1,-1\}.       \tag{2.5}
\]

The two signs must agree when \(\gcd(m,n)>2\), giving a kernel of order two.
When \(\gcd(m,n)\le2\), all four sign pairs are compatible and distinct.
Therefore

\[
 D=
 \begin{cases}
 \varphi(\ell)/2,&\gcd(m,n)>2,\\
 \varphi(\ell)/4,&\gcd(m,n)\le2.
 \end{cases}                                             \tag{2.6}
\]

The divisibilities \(d_m\mid D\) and \(d_n\mid D\) are checked for the
whole finite box.

### 2.3 Exact pair sieve

Evaluate (2.2), with (2.6), for every unordered pair
\(3\le m\le n\le840\). All comparisons are integer comparisons; no
floating sine or logarithm is used.  Exactly 525 pairs survive.  Their
largest label is 420, and the only endpoint survivor is

\[
                         (m,n)=(3,420),\qquad D=48.       \tag{2.7}
\]

This proves the rigorous replacement

\[
                    \boxed{3\le m,n\le420}.              \tag{2.8}
\]

It does **not** prove that any of the 525 pairs supports an integral
\(\beta\), an admissible field, or a Coxeter polyhedron. It is a complete
first-stage candidate list on which Fekete/minimal-polynomial enumeration can
now operate.

## 3. Exact reduction of the tail label

For the tail coordinate, the existing proof gives

\[
 1<196^2{4\over\ell_{d,0}^2}
       \left({\gamma(d)^2\over4^{\varphi(d)/2}}\right)^e.
                                                               \tag{3.1}
\]

Put
\[
 T_d={\gamma(d)^2\over4^{\varphi(d)/2}}.
\]
If \(T_d<1\), then \(T_d^e\le T_d\), and
\(\ell_{d,0}\ge16/d^2\) turns (3.1) into the necessary inequality

\[
             601d^4\gamma(d)^2>4^{\varphi(d)/2}.        \tag{3.2}
\]

If \(T_d\ge1\), the same inequality is automatic because
\(601d^4T_d>1\). (The only labels with \(T_d>1\) are \(d=3,5\), and
\(T_4=1\).) Hence (3.2) is again a valid necessary filter for
every \(d\), without a false global contraction claim.

The exact integer sieve over \(3\le d<20000\) leaves 66 values, the largest
being \(126\). The analytic endpoint and monotonicity argument already in
the minimality-\(14\) field theorem excludes \(d\ge20000\). Hence

\[
                         \boxed{d\le126}.                \tag{3.3}
\]

Again this is a necessary-filter result, not an existence assertion.

## 4. Why 14, 62, and 995 behave differently

### 4.1 Minimality 14 is locally sharp

For two closed facets meeting a minimum-aperture facet, the normalized
pairing has the form \(4X-2\), with each of the two sine ratios defining
\(X\) strictly below two. Thus the pairing is below 14. Setting the three
cone angles equal to \(\varepsilon\) gives

\[
        16\cos^2(\varepsilon/2)-2\longrightarrow14.
\]

So 14 is the sharp supremum of this local cone estimate.  A smaller
universal constant would require extra Coxeter/arithmetic restrictions or a
different global choice of frame, not a more careful version of the same
inequality. Computationally, the useful optimization is therefore to prune
the integral points inside the minimality-\(14\) charts, as Sections 2--3 do.

### 4.2 The frame constant 62 is imported and nonoptimized

Nikulin's geometric appendix lemma says that one *can take* \(A=62\) for a
connected full-rank wall frame.  The native-host proof uses only that an
absolute finite \(A\) exists. Nothing in its cycle-packet, Burnside, or
conductor-sandwich stages improves 62.  Reducing it would dramatically shrink
the selected-place Minkowski boxes \(62^r\), but it requires a new geometric
frame lemma.  Replacing 62 by 14 merely because a narrow facet exists is not
valid: the closed star need not already contain a connected Lorentz basis in
the noncompact/general frame setting.

### 4.3 The 995 ceiling is a geometric termination device

The number 995 is imported wholesale from the Khovanskii--Prokhorov
finite-volume Coxeter-polytope theorem.  The field funnel cannot lower it:
spherical components of types \(A_r,B_r,D_r\) have unbounded rank while
their labels and ground field stay fixed.

For the **arithmetic** scope there are two optional classical shortcuts,
deliberately absent from the AG-mediated independence theorem:

1. importing Vinberg's crystallographic nonexistence theorem replaces the
   global \(n\le995\) by \(n\le29\);
2. in the noncompact branch the cusp already forces \(k=\mathbf Q\), and
   Esselmann's reflective integral-lattice theorem restricts that branch to
   \(n\le19\) or \(n=21\).

Thus a computational version should normally use ranks at most 30 globally
and at most 22 in the rational noncompact branch.  This is a major practical
improvement, but it is a **classical Vinberg/Esselmann import**, not an AG
achievement.  If independence from arithmetic dimension theorems is part of
the theorem statement, 995 remains the honest current ceiling.

## 5. Computable classification route

The constants now suggest a staged algorithm rather than a scan of every
graph up to rank 996.

1. Enumerate the 525 fork pairs and 66 tail labels, then enumerate integral
   determinant coordinates in their contracting real boxes using exact
   relative minimal polynomials, norm pruning, and Fekete degree bounds.
2. Make the triangle/quadrilateral/orbicurve branch effective; the current
   use of Long--Maclachlan--Reid proves finiteness but does not itself output
   the list.
3. Work one surviving \((k,m)\) at a time. Enumerate connected Cartan supports
   under Coxeter, inertia, cusp, and nonselected-definiteness constraints
   before enumerating cycle values.  A raw search through every
   \(62^m\)-sized Minkowski box is not realistic.
4. Reconstruct representations by spanning-tree gauge, enumerate stable
   projective lattices via the full order and conductor quotient, and run a
   certified Vinberg algorithm to recover all mirrors and test finite volume.
5. Quotient by field, lattice, diagram, and chamber isomorphisms, and test
   reflection maximality by full mirror-host equality.

For discovering examples, the best immediate lanes are the now finite
compact-\(H^3\) packet list and fixed small \((k,m)\) slices. For lowering a
dimension ceiling rather than finding examples, the genuine open geometric
targets are a better connected-frame lemma and rational-cusp link/diagram
obstructions.  Better label sieves alone do not imply a better dimension
bound.

## 6. Reproduction

```text
python3 scripts/certify_coxeter_finiteness_constant_improvements.py
python3 -m pytest -q tests/test_coxeter_finiteness_constant_improvements.py
```

The certificate contains all 525 fork-pair survivors and all 66 tail-label
survivors.  Its comparisons are exact integers.  The geometric provenance of
14, 62, and 995 is audited separately in
[`geometry_vs_ag_provenance_ledger.md`](geometry_vs_ag_provenance_ledger.md).
