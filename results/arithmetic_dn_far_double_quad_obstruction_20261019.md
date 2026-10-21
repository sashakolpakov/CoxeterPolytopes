# The far-double quadrilateral continuation is dominated by the far-fork case

19 October 2026. **New written research corollary, awaiting independent
review.** No computation or build was run. Existing reports are unchanged.

Use the notation and actual-root hypotheses of the
[far-fork note](arithmetic_dn_far_fork_quad_obstruction_20261019.md),
except that the new side d now contacts the far long endpoint lL with
label 4, rather than the next-to-end node with label 3. Thus

    d.lL=-sqrt(2)/2,
    d is perpendicular to every other old root except b,
    c.b=c.d=0,
    b.d=-beta,     beta>=1.

The old vertex is whole D19 or D20; N=n-2, L=N-1 and
W={e,l1,...,lL}. The conclusion is that this actual quadrilateral
cannot occur. The proof uses no arithmetic input.

The finite opposite corner W+c+d contains a whole large B factor
on l1,...,lL,d, with its only 4-bond at the d end, and an isolated e.
Any finite extension meeting that B factor is a larger B path.
Consequently c can either have no long-arm contact, with c.e any
finite dihedral bond or zero, or it can have a simple contact to l1.
In the latter case c.e is zero or a simple bond. An additional
4-bond, an internal 4-bond, a fork or a cycle is impossible in a
finite component of these ranks. This is exactly the c inventory
in Section 2 of the far-fork note.

After projecting off W, the square s of a and the possible squares
t of c, as well as X=-a_bar.c_bar, are the same as in that note.
The quantities involving d are instead

    u=1/N,     r=sqrt(2)/(2*N),
    v=0       if c has no long-arm contact,
    v=sqrt(2)/(2*N) if c.l1=-1/2.

Here r=-a_bar.d_bar and v=-c_bar.d_bar. Compare these with the
far-fork values

    u_*=2/N,     r_*=1/N,
    v_*=0 or 1/N,

respectively. We have u<=u_*, r<=r_* and v<=v_*, while t,s,X and
beta are held fixed.

The exact positive-root expression from that note is

    X = v*(beta/2+r)/(beta^2-u)
        + sqrt[(t+v^2/(beta^2-u))
               *(s+(r*beta+r^2+u/4)/(beta^2-u))].

For beta>=1 and 0<=u<=u_*<1, its right side is coordinatewise
nondecreasing in u,r,v: all numerators are nonnegative, increasing
one of these parameters increases every affected numerator, and
increasing u decreases the positive denominator. Thus X is at most
the far-fork upper expression at the same t,s,beta. That expression
is in turn bounded by its value at beta=1, as proved there.

Every numerical inequality in Sections 4--5 of the far-fork note
therefore applies unchanged. If c has no long-arm contact, it forces
either the impossible completely orthogonal c, the affine D4 star,
the affine E8 arms (1,2,5), or the affine B4 witness. If c contacts
l1, it forces a.c=0 and hence the affine D5 or affine A3 witness.
Their roots do not include d, so their null equations and all their
required zero contacts are unchanged by the present replacement.
The positive-null compactness argument excludes all cases.

This closes beta>=1 for the far-end 4-bond continuation. Its finite
beta branch is separately excluded in the
[finite-product note, type B](arithmetic_dn_nonisolated_quad_finite_opposite_obstruction_20261019.md).
The combination still requires independent review. The remaining
unpaid quadrilateral issue is the far-end simple continuation;
simultaneous weight consumers have not been settled by this corollary.
