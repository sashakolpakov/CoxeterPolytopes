# Isolated side with a long-arm contact: zero and finite opposite products

19 October 2026. **Written local theorem candidate, awaiting independent
review.** Paper/source work only. No scientific process, enumeration,
optimizer, Lean or build ran, and no pinned source was edited. This
source treats x=0 and finite native x only. It does not exclude the
dotted long-contact cases or install a spectator half refund.

## 1. Whole-factor scope and complete local inventory

Use the exact ambient20 eight-field setup and actual quadrilateral
a,b,c,d of
`arithmetic_rank20_isolated_side_no_long_contact_20261019.md`:
the entire old vertex is D_r+T_s, r+s=20,10<=r<=16,
N=r-2>=8,L=N-1, and the complete white support is
A1(e)+A_L+T_s. The old D_r fork is a and its short leaves are b,e.
Assume d contacts only the old leaf b, c.b=c.d=0, and c is
orthogonal to the ENTIRE spectator T_s. Now require c to have a
nonzero contact with the D_r long arm. Write a.c=-x.

**Claim.** The actual quad is impossible if x=0 or x=cos(pi/j)
for any finite native integer j>=3.

The complete finite W+c inventory is:

* one simple contact at k=1,2,N-2,N-1, and no e contact;
* at N=8, one simple contact at k=3 or5, and no e contact;
* one endpoint4 contact, and no e contact;
* one simple endpoint contact, and an e3 or e4 contact.

Two long contacts make a cycle. A4 bond must be terminal in a
finite path in these ranks. Adding e to a next-to-end fork, a
long4 bond, or the exceptional E8 cases is not finite. The
remaining internal simple attachments give nonfinite arms except
for the two E8 extensions explicitly retained in the second row.

## 2. Important small-rank qualification

At L=7 a classical-only inventory would be FALSE:
attaching c simply at l3 or l5 gives arms1,2,4, the finite E8.
At L=8, attaching at l3 or l6 gives arms1,2,5, which is affine;
at larger L there is no exceptional finite one-root extension.
Thus r=10 also has the extra c.l3 or c.l5 simple contacts, with
c.e=0. They MUST be retained in the centered-moment table below.
They cannot also meet e, since the component would then have rank9.

Accordingly the theorem's full inventory includes k=3,5 at N=8.
No actual exceptional contact is suppressed by an asymptotic
classical classification. This is an essential distinction from
the previous whole-D19/20 sources.

## 3. Common arithmetic identities

The spectator is orthogonal to all roots used, so the proper-block
identities apply at the actual full vertex without assuming a
section is a Coxeter polytope. With c.lk=-alpha,c.e=-gamma,

    t=1-gamma^2-alpha^2*2k(N-k)/N>0,
    kappa=gamma/2+alpha*(N-k)/N,
    qc=1-t+2N*(x+kappa)^2.

For a simple long contact alpha=1/2. The actual old-path cycle
puts2x in O_K. When gamma is0 or1/2, t,kappa are rational,
and the independent old-full-vertex-plus-c block gives

    |x'+kappa|<rho,       rho^2=t/(2N).

It also proves K=Q(x), since an embedding fixing x would preserve
qc>1. The cosine moment estimate

    Tr((x+kappa)^2)>=d*(1/2+kappa^2)-1/2

compared to the displayed interval gives d<=3 whenever

    B=1/2+3*kappa^2-2*kappa-3*rho^2>=0.

The weak inequality B>=0 suffices because the conjugate upper
bound is strict. Thus only j=3,4,5,6,7,9 remain. Exact field equality
in the eight fields then leaves only j=6 (Q(sqrt3)) or j=9 (K9).

## 4. Simple long contact, including the actual E8 exceptions

With gamma=0 the classical four rows have

|k|kappa|t|B|
|---|---|---|---|
|1|(N-1)/(2N)|(N+1)/(2N)|(N-5)/(4N)|
|2|(N-2)/(2N)|2/N|(N-4)/(4N)|
|N-2|1/N|2/N|(N-4)/(2N)|
|N-1|1/(2N)|(N+1)/(2N)|(2N-7)/(4N)|

All B are positive for N>=8. The exceptional N=8 rows are

    k=3: kappa=5/16, t=1/16;
    k=5: kappa=3/16, t=1/16.

Their B values are positive as well: respectively
1/2+3*(5/16)^2-2*(5/16)-3/256=5/32,
and1/2+3*(3/16)^2-2*(3/16)-3/256=7/32.
Their rho=1/16. All rows have kappa<1/2 and rho<=3/16.
For j=6 use the conjugate-sqrt3/2; for j=9 use cos(7pi/9),
which is below-sqrt2/2. The shifted absolute value is greater
than sqrt2/2-1/2>3/16, contradicting the required interval.
Thus every simple-long-only finite product is excluded, including
the two finite E8 continuations at r=10.

For an additional simple e bond, only the two endpoint rows occur:

    near: kappa=(3N-2)/(4N), t=(N+2)/(4N),
          B=(11N-26)/(16N);
    far:  kappa=(N+2)/(4N), t=(N+2)/(4N),
          B=(3N-10)/(16N).

Again only j=6,9 remain in the eight fields. In the far row the
negative conjugates just used have shifted absolute value greater
than sqrt2/2-5/16, whereas rho<=sqrt5/16. For the near j=9 row
use cos(5pi/9)>-1/2; its shifted value exceeds3/16>rho.
For near j=6 the exact squared difference is positive:

    16N^2*((sqrt3/2-3/4+1/(2N))^2-rho^2)
      =N*((21-12sqrt3)*N+8sqrt3-14)>0.

At N>=1 the last bracket is at least7-4sqrt3>0. This excludes
all simple endpoint plus e3 finite products.

## 5. Endpoint4, with no e contact

Put w=2*sqrt2*x, an integral actual cycle scalar. Direct inversion
gives

    qc=1-1/N+N*(w/2+kappa0)^2,
    kappa0=(N-k)/N.

Thus K=Q(w). At the far endpoint every nonphysical
-4/N<w'<0. Physical sqrt2<=w<2sqrt2<3 gives degree at most2,
since3*(4/N)^2<=3/4<1. In degree2 its norm is-1 and trace1
or2, giving exact fieldsQ(sqrt5) or Q(sqrt2), outside the eight.
In degree1 w=2 forces Q. The far4 case is excluded.

At the near endpoint set A=w+2. Physically2+sqrt2<=A<5,
and every other0<A'<4/N<=1/2. Degree at least4 has norm
below5/8, so d<=3. If d=2, its positive norm is1 or2 and trace4
or5. The four discriminants give fieldsQ(sqrt3),Q(sqrt21),
Q(sqrt2),Q(sqrt17). Only the first is one of the eight, and it
would give x^2=3/8, not a native finite cosine-square.

If d=3, the field must be K9. Since w^2=8*x^2 belongs to K9,
its native cosine-square field is Q or K9. The possible labels are
j=3,4,6,9,18. Their exact norms in K9 of8*x^2 are respectively

    8,64,216,8,24.

Only64 is a rational square; it gives w=2 and contradicts d=3.
For the others Norm(w^2) cannot be the square of the integer
Norm(w). The two nonrational entries follow from
Norm(4*cos^2(pi/j))=Phi_j(-1), equal to1 or3 for j=9 or18.
This excludes K9 without a hidden integral-basis assumption.
Degree1 again gives w=2 and Q. Thus near4 is excluded too.

## 6. Simple endpoint plus e4: the exact small-N CRT bound

The old a-e-c-long-path cycle puts sqrt2 in K, and the other
actual cycle then puts x in K. Nonphysical positivity requires

    |x'+kappa0+sqrt2'/4|<1/(2N),
    kappa0=(N-k)/(2N)<=1/2-1/(2N).

If4 does not divide j, the cyclotomic automorphism fixing x and
negating sqrt2 makes the left side at least1/2-sqrt2/4>1/8.
If j=4 modulo8, exponent j-1 negates both x and sqrt2 and
contradicts the inequality immediately. If8 divides j, the same
exponent negates x and fixes sqrt2. Now

    cos(pi/8)>1/2+sqrt2/4,

because the difference of their squares is1/8. Therefore
x-kappa0-sqrt2/4>1/(2N), again a contradiction. This exact
endpoint shift is what permits the argument at N=8; the old
coarse3/70 comparison would not suffice there. Each automorphism
of the indicated cyclotomic subfield extends to the totally real K.
All finite products are now covered.

## 7. The zero x cases and the unresolved dotted branch

When x=0, a simple contact at any k<L gives the positive-null D
tree with path a-l1-...-lk and four leaves b,e,c,l_(k+1), so this
includes the exceptional k=3,5. A far simple contact gives positive
D_(r+1)+T_s of rank21. An endpoint4 gives the positive-null B
tree. Adding e3 gives the affine simple cycle; e4 at the far end
gives the five-node affine B witness, and e4 at the near end has
conjugate cost r/4+(r-3)/2-(r-3)/sqrt2>1 for r>=10.
All are principal actual-wall obstructions.

The present theorem thus covers zero and finite native x with EVERY
actual long contact, including the previously easy-to-miss E8
extensions at r=10. The dotted x>=1 branches remain unresolved
here. Spectator contacts are explicitly absent; their separate
finite continuation obligations have not been imported into this
statement. No half refund, forest cap or dimension bound is claimed.
