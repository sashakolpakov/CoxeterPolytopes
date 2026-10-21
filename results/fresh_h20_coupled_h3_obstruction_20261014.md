# Coupled H3 sections: exact irrationality and compact-face obstruction

No new finite-volume polytope is claimed. This note excludes the coupled
two- and three-contact H3 configurations described below, for every real
section normal, in the saved 210-wall Borcherds H21 chamber. It also adds
the one-contact B3/B4 endpoint cases to the preceding
[A3/A4 obstruction](fresh_h20_h4_face_closure_20261014.md).
Together these cases prove that **no finite-volume Coxeter H20 section of
this source can contain an actual H3 face**. This is not a classification
of all H20 sections: a section might have no H3 face.

The [exact producer](../scripts/fresh_h20_coupled_h3_obstruction_20261014.py)
and [algebra certificate](fresh_h20_coupled_h3_obstruction_20261014.json)
use rational and Q(sqrt(5)) arithmetic and irreducible polynomials over Q.
The separate numerical H-type scout is not a premise of this proof.

## Geometric reduction

Write a section normal as e with e²=1, and let a_i=(e,r_i). The projected
wall Gram matrix is H_ij=G_ij-a_i a_j. An actual spherical H3 subdiagram
bounds a compact H17 face F of a hypothetical finite-volume H20 section:
H3 has no affine extension. Its only proper connected finite extension is
H4, by a label3 edge at the end opposite the label5 edge. These standard
finite/affine diagram and Coxeter-face facts are the same premises used in
the preceding A3/A4 proof. In particular this H3 face is itself Coxeter.

If every boundary wall of F is an unchanged integral source root, then
compactness makes these roots span F's Lorentzian linear space U. Hence U
is rational, and its full lattice U intersect Z^22 is preserved by every
boundary reflection. The compact reflection lattice has finite index in
the corresponding integral orthogonal group. Meyer isotropy in rank at
least5 makes that orthogonal group nonuniform, a contradiction. See the
[previous proof](fresh_h20_h4_face_closure_20261014.md) for this lattice
step. Rational coefficients alone would not suffice: lattice preservation
and the spanning argument are essential. The producer rechecks all210
source reflection matrices for integrality.

Thus it suffices to show that an orthogonal boundary of F has (e,r)=0,
and that a nonorthogonal H4 extension is impossible.

Every old Gram block supporting a new spherical block is positive
definite, since G=H+aa^T. The source finite diagrams are crystallographic
forests. Its exact pair pools are short/short {-1,0}, short/long {-2,0},
and long/long {-8,-6,-4,-2}, with norms2 and4. The five oriented old
three-node paths (center c, leaves f,l) therefore have
(d_c,d_f,d_l,g_cf,g_cl) equal to

```
(2,2,2,-1,-1), (2,4,2,-2,-1), (2,2,4,-1,-2),
(4,2,4,-2,-2), (4,4,2,-2,-2).
```

## Joining A1 to a rank-two component

Let the old isolated root i have norm d0 and the old adjacent roots a,b
have norms da,db and squared cosine c=1/4 or1/2. Take contacts (u,v,0)
on (i,a,b), producing the path b--5--a--3--i. Put C=(3+sqrt(5))/8 and

```
Y=v²/da=1-c/C,       X=u²/d0=(1-Y)/(1+3Y).
```

The field norm of X/Y in K=Q(sqrt(5)) is1/11 or1/5. Multiplication by
the positive rational d0/da cannot make it a square in K. Consequently
r=u/v is outside K. An orthogonal boundary wall has integer pairings
g_i=u t and g_a=v t, where t=(e,r_wall); their irrational ratio forces
t=g_i=g_a=0.

An H4 extension must be orthogonal to a,b and attach at i. Orthogonality
to a gives t=g_a/v; the label3 condition is

```
4(g_i-r g_a)² = d0(1-X)(d_wall-g_a²/v²).
```

All terms apart from the linear r term lie in K, so g_i g_a=0. If g_a=0,
the equation equates a rational number to the irrational d0(1-X)d_wall.
If g_i=0, it becomes g_a²=d_wall*da*Y², also irrational. No H4 extension
exists; the rational compact-face contradiction applies.

## Two contacts on an old connected path

Let a_l=0 and a_c,a_f be nonzero. The unchanged leaf l must lie at the
label5 end: a nonzero a_c strictly increases its old squared cosine,
which was1/4 or1/2. Write A=a_c² and R=a_f/a_c. The two target angles give

```
A = d_c-g_cl²/(d_l*C),
A(3A+d_c)R² - 8g_cf*A*R + 4g_cf² - (d_c-A)d_f = 0.
```

For each of the five old paths, the norm of the monic quadratic in R
from K to Q is irreducible of degree4. Eliminating R likewise gives an
irreducible degree4 polynomial for a_f²=A R². Therefore R and a_f² are
irrational. Orthogonal boundary walls have t=0 by the irrational ratio R.

For an H4 extension m, orthogonality to c,l gives g_cm=a_c t, g_lm=0;
m attaches at f. If t=0, its label3 condition would make a_f² rational,
which is impossible. If t is nonzero, the old positive-definite
four-node forest has the three edges cf,cl,cm. It must be D4, with fm
absent. All norms must be2: a label4 would make the star non-positive,
and the all-long D4 is absent from the literal source because distinct
long roots are never orthogonal.

For this sole possibility, set

```
A=2-1/(2C), T=t²=1/A, B=a_f²=2(2-T)/(2+3T),
S=((2-A)(2-B)/4-1-A B)/2.
```

The two label3 conditions require S=a_c a_f and S²=A B. Exact algebra
instead gives S²-A B = 27657/1444 - (3122/361)sqrt(5), which is nonzero.
Thus there is no H4 extension and again F would be compact rational
arithmetic.

## Three nonzero contacts

For a new H3 path, three nonzero contacts on an old finite forest have
mixed signs: otherwise a triangle is created. The old path has a center
with sign opposite both leaves, and one old center/leaf edge must cancel.
Choose the canceled edge cf, orient a_f positive, and write

```
a_f=sqrt(X), a_c=g_cf/sqrt(X), a_l=R sqrt(X), R>0.
```

Let C1 be the target squared cosine on fl and C2 on cl, so
{C1,C2}={C,1/4}. The angle equations imply

```
R² = C1*d_l*(d_f-X) / [X*(C1*d_f+(1-C1)*X)],
c0 = g_cl²+g_cf² R²-C2*(d_c-g_cf²/X)*(d_l-X R²),
P(X) = numerator(c0²-4g_cl² g_cf² R²) = 0.
```

The producer treats both orientations for each of the five old paths.
All ten monic P are quartic over K; their norms to Q are irreducible of
degree8. In particular X is irrational. Eliminating X gives degree8
irreducible norm polynomials for R² when C1=C, and for X²R² when C2=C.
All displayed denominators are nonzero in a physical solution (positive
projected norms, X>0); clearing them or squaring only adds candidates,
so the resulting irrationality restrictions remain necessary.

It follows that the ratio of contacts on the label5 endpoints is
irrational. An H4 extension is orthogonal to both those endpoints, so
its contact t must vanish, and its projected norm remains an integer.
But its label3 attachment would then require a rational squared projected
norm at the third H3 root. That norm is d_c-g_cf²/X if C1=C, and d_f-X
otherwise, both irrational. The extension is impossible. The same
irrational contact ratio forces t=0 for orthogonal boundary walls.
This proves the obstruction for all these continuous three-contact cases.

## One-contact B3/B4 endpoint cases

An old label4 edge becomes label5 by modifying its endpoint a of norm
d_a. Here normalize (e,r_a)=1 rather than e²=1; then
q=e²=(2+sqrt(5))/d_a. A boundary root orthogonal to the projected B3/B4
block has (e,r_wall)=q*g_awall and projected norm d_wall-q*g_awall².
Source integrality/parity makes this negative for every nonzero
g_awall. Thus every such boundary is unchanged. An H4 face is compact
and has no proper finite connected extension, giving the rational-face
contradiction immediately. For H3, a nonorthogonal boundary attaches to
its unchanged far end; its label3 equation first forces g_awall=0 and
then gives exactly an old B4 extension, already excluded. This closes
the B3 case too.

## Scope and verification

The cases exhaust the possible contacts on the three walls of an actual
H3 face, without restricting contacts on any other source wall. The old
three-wall Gram is positive definite and its diagram is crystallographic.
Zero contacts cannot produce H3. One contact must be at an endpoint of an
old A3/B3 path: at the center both old nonright angles would increase in
squared cosine, whereas one target edge must still have label3. With two
contacts, a disconnected old diagram must be the isolated-root/rank-two
join. On an old connected path the contacts cannot be at both leaves:
equal signs create a triangle and opposite signs give a positive projected
pairing. Thus they are at its center and one leaf. The three-contact case
was reduced above to a mixed-sign old path with one canceled edge.

Consequently, together with the preceding one-contact A3 proof, **every
actual H3 face is excluded**, without assuming that the hypothetical
section is arithmetic or restricting its other normal pairings. Every
H4 diagram contains H3, so actual H4 faces are excluded as well. A separate
argument would be needed to show that every possible section has an H3
face; no such claim is made here.

The saved exact calculation checks30 irreducibility claims (ten quartic
and twenty degree8 norm polynomials), the two joined-field norm gates,
the D4 nonzero defect, the B-endpoint inequalities, and literal source
integrality/pair pools. These are exact algebra checks, not yet an
independent implementation of the entire proof. The
[separate review](fresh_h20_coupled_h3_obstruction_20261014_review.md)
records independently reconstructed identities and the completeness check.
