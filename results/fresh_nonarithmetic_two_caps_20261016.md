# The proposed two-cap family fails at its first cap

All **128 ordered recipes** in the frozen [nonarithmetic target](fresh_nonarithmetic_escape_target_20261016.md) are excluded at their first stage. The eight first-cap recipes have exactly sixteen isolated real branches and no degeneracies. Nine branches have a positive pairing with a retained wall. The other seven have a finite angle strictly between consecutive Coxeter angles. There are **no accepted or unresolved first caps**, so no second-cap systems were attempted.

This is a complete obstruction to that bounded recipe family, not to all caps perpendicular to the support or to a finite-volume completion of the retained chamber. In particular, it does not impose a finite-angle alphabet on unspecified pairings: each rejection is justified by acuteness or by a gap valid for every integer Coxeter label.

The [exact checker](../scripts/fresh_nonarithmetic_two_caps_20261016.py) and [saved certificate](fresh_nonarithmetic_two_caps_20261016.json) pass generation and saved replay. They reuse the frozen independent octic arithmetic of the [33-wall bridge review](fresh_h9_joint_corner_20261015_bridge_review.md), without changing the old target or any chamber artifact. This producer is not itself an independent implementation of that arithmetic.

## Complete first-stage solution

Retain all33 old walls and let

\[
T=[0,1,2,5,6,10,30].
\]

For \(j\in\{3,7\}\) and \(m\in\{2,3,4,5\}\), the first cap \(u\) is required to satisfy

\[
(u,n_i)=0\ (i\in T),\qquad
(u,n_{31})=0,\qquad
(u,n_j)=-2\cos(\pi/m),\qquad u^2=2.
\]

For each anchor the nine linear covectors are independent. The checker constructs their one-dimensional kernel \(k\), and the unique particular vector \(p\) satisfying

\[
(p,n_i)=0\ (i\in T\cup\{31\}),\quad
(p,n_j)=-1,\quad (p,k)=0.
\]

It checks \(k^2>0\) and \(p^2<0\). Set

\[
s_m=4\cos^2(\pi/m)\in\{0,1,2,\phi^2\},\qquad
d_m=\frac{2-s_m p^2}{k^2}>0.
\]

The complete real solution set is therefore

\[
u_\pm=\sqrt{s_m}\,p\ \pm\sqrt{d_m}\,k.
\]

Both square roots have their nonnegative real meanings. Since \(d_m>0\) and \(k\ne0\), the branches are distinct. The independent linear rank and nonzero kernel square exclude continuous solution sets and singular conics. These formulas also treat \(m=2\) exactly, with \(s_2=0\).

The existing support \(T\cup\{31\}\) has positive-definite rank-eight Gram matrix. The required new root is perpendicular to all eight normals and has norm2, so its proposed rank-nine endpoint Gram is automatically positive definite. This necessary finite-corner condition does not repair the incompatibility with other retained walls found below; no new actual endpoint is claimed.

## All sixteen branches

The signs below use the kernel normalization serialized in the certificate. An entry “gap \((a,b)\) at wall i” means that the Cartan cycle lies strictly between the values for consecutive integer labels \(a\) and \(b\).

| Anchor j | Prescribed label m | Negative branch | Positive branch |
|---:|---:|---|---|
| 3 | 2 | obtuse to7 | obtuse to4 |
| 3 | 3 | gap(2,3) at7 | obtuse to29 |
| 3 | 4 | gap(4,5) at7 | obtuse to29 |
| 3 | 5 | gap(7,8) at7 | obtuse to29 |
| 7 | 2 | gap(2,3) at3 | obtuse to3 |
| 7 | 3 | gap(3,4) at3 | obtuse to3 |
| 7 | 4 | gap(3,4) at3 | obtuse to4 |
| 7 | 5 | gap(4,5) at3 | obtuse to8 |

All528 old-wall pairing signs are recomputed and saved. For an obtuse branch, the displayed wall is one exact positive blocker; the certificate also lists every other positive pairing. Every angle-gap branch is acute to all33 old walls and negative at the saved strict timelike center. Its failure is specifically the non-Coxeter angle at the other anchor.

For a retained normal \(n\), write

\[
(u,n)=\sqrt{s}\,a+\sqrt d\,b.
\]

If the two terms have opposite signs, their sum has the sign of \(a\) times the sign of \(s a^2-d b^2\). The zero-term and same-sign cases are immediate. This reduces every comparison to the old exact field; no primitive element for the new square roots is needed. The saved certificate records the coefficients, signs and squared-magnitude comparison for each displayed blocker.

For the acute branches, the squared normalized pairing is

\[
c=\frac{4(u,n)^2}{u^2n^2}
 =A+B\sqrt{sd},\qquad
A=\frac{2(sa^2+db^2)}{n^2},\quad B=\frac{4ab}{n^2}.
\]

Each comparison with a field element uses the same exact sign rule. Six branches lie strictly in one of

\[
(0,1),\qquad (1,2),\qquad(2,\phi^2),
\]

which are precisely the gaps between labels2/3,3/4 and4/5. Their negative nonsingular pairings therefore cannot be Coxeter relations, and they are not parallel or dotted.

## The label7–label8 comparison

The remaining branch, anchor3 with prescribed label5, has

\[
\frac{33}{10}<c<\frac{10}{3}.
\]

These are exact radical-sign comparisons, not numerical thresholds. Let \(c_7=4\cos^2(\pi/7)=2+2\cos(2\pi/7)\). The seventh-root-of-unity identity gives

\[
P(c_7)=0,\qquad P(X)=X^3-5X^2+6X-1.
\]

Also \(c_7>c_6=3\). On \([3,\infty)\), the derivative is positive: \(P'(3)=3\) and \(P''(X)=6X-10\ge8\). Since

\[
P(33/10)=287/1000>0,
\]

we have \(c_7<33/10\). Meanwhile

\[
c_8=2+\sqrt2>10/3
\]

because \(2-(4/3)^2=2/9>0\). Thus \(c_7<c<c_8\). The function \(4\cos^2(\pi/m)\) strictly increases for real \(m\ge2\), so no integer finite label is possible. This argument makes no assumption about which cosine fields are contained in the new radical extension.

## Scope and next distinction

The first stage has eight recipes, two distinct branches per recipe and no rank-deficient case. Every branch is rejected, with no unresolved algebraic equality. Consequently all128 proposed ordered two-cap recipes are excluded regardless of the second endpoint and mutual labels. Solving their second stages cannot repair an old-wall violation of the first cap.

The negative branches with anchor3/label2 and anchor7/label2 do pair positively with the saved escape, but are respectively obtuse to wall7 and in the forbidden label2/3 gap at wall3. No compatible escape-cutting normal was obtained. The retained33-wall chamber and its strict escape are unchanged.

Other prescribed anchor labels, other anchors, nonright contacts to T, different retained supports, changed old walls, larger cap clusters and arbitrary nonarithmetic completions remain open. No finite-volume or general nonexistence claim follows from this bounded obstruction.

```sh
python3 scripts/fresh_nonarithmetic_two_caps_20261016.py
python3 scripts/fresh_nonarithmetic_two_caps_20261016.py --verify
```
