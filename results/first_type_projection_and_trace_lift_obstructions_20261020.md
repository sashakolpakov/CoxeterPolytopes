# Two exact obstructions to transferring the right-angled method to first type

20 October 2026. Paper algebra only. No matrix program, optimizer, census,
checker, Lean or build was run. The purpose is to identify two precise
missing structural steps, not to claim a new compact example above
dimension eight or a better uniform dimension bound.

## 1. Maximum-height projection is not a contraction at a nonright angle

Let q be a finite interior unit point, and let n be a maximum-height
outward unit normal, with `t=-<q,n>`. If an adjacent unit normal m has
`<m,n>=-c`, where `c=cos(pi/l)`, its intrinsic unit normal in n-perp is

    mbar=(m+c*n)/sqrt(1-c^2).

The actual orthogonal projection of q to n-perp is still
`q'=(q+t*n)/sqrt(1+t^2)`, and all-wall acuteness still proves feasibility.
But the NEW intrinsic height is

    -<q',mbar>=(t_m+c*t)/(sqrt(1+t^2)*sqrt(1-c^2))
      <=cot(pi/(2*l))*t/sqrt(1+t^2).                         (1)

Thus the unit-normal contraction used for right angles acquires a factor
sqrt3 already at label three. Arithmeticity of first type does not remove
this factor, as the following genuine compact arithmetic triangle shows.

Put r=sqrt2 and take the norm-two root Gram matrix

    B=[[2,-1,-1],[-1,2,-r],[-1,-r,2]].

Its determinant is `-2r`; its first two leading principal minors are
2 and 3. Hence it has signature (2,1). The conjugate r->-r has the same
first minors and determinant `2r>0`, so is positive definite. The root
reflections are integral over Z[r], because every root has norm two and
all Gram entries are integral. The hyperbolic triangle has angles
pi/3,pi/3,pi/4, hence is compact. The Gram cycle generates Q(r); these
integrality and conjugate-positivity facts also verify the usual Vinberg
arithmeticity criterion. Its reflection group is therefore of first type
over Q(sqrt2), not merely a formal local Gram configuration.

For G=B/2, solving `G*c=1` gives

    c0=-2-3r,  c1=c2=-3-3r,  sum c_i=-8-9r.

The incenter has all unit-normal heights equal to t, with

    t^2=1/(8+9r).

Choose the facet whose two adjacent labels are three; it is a legitimate
maximum-height facet because all heights tie. Formula (1) is then an
equality at both neighbors, and gives

    (t')^2=3*t^2/(1+t^2)=(r-1)/3 > t^2.

The last comparison is `(10-r)/3>1`. Therefore even this compact
arithmetic example increases the maximum intrinsic height under the
proposed projection. This does not rule out a large-dimensional theorem
using additional structure; it rules out obtaining that theorem merely
by retaining the right-angled recurrence for first-type forms.

## 2. Even an integral Coxeter trace lift need not have finite volume

There is a second, independent obstacle to applying a rational reflective
lattice dimension bound after restriction of scalars. Here is an explicit
compact first-type chamber whose natural split-root trace chamber has an
unsaturated cusp.

Over Q(r), let a1,a2,a4 be a basis with norm-two Gram matrix

    B0=[[2,-1,-1],[-1,2,-2-r],[-1,-2-r,2]],

and put `a3=r*(a2+a4)-a1`. Direct substitution gives

    a_i^2=2,
    <a_i,a_(i+1)>=-1  (cyclic indices),
    <a1,a3>=-2-2r,   <a2,a4>=-2-r,
    a1+a3=r*(a2+a4).                                        (2)

The determinant of B0 is `-12-10r<0`, with leading minors 2 and 3.
Its conjugate determinant `-12+10r>0` makes the conjugate form positive
definite. All four root reflections preserve the Z[r]-lattice. The vector
`v=a1+a3` has square `-4r` and products `-2r,-2,-2r,-2` with the four
walls, so is a strict feasible timelike point. Adjacent pairs give four
actual finite vertices, and opposite pairs are ultraparallel. Thus these
walls enclose a compact right-angled-free quadrilateral with all four
angles pi/3. The cycle entries generate Q(r); again the integral roots
and positive conjugate verify first-type arithmeticity.

Now view this three-dimensional K-space as a six-dimensional Q-space
with trace bilinear form `T=Tr_(K/Q) B0`, of signature (5,1). Introduce
the eight positive roots

    a1,...,a4,    b1,...,b4,   b_i=r*a_i.

Their T-norms are 4 and 8 respectively. For a lattice vector x with
`B0(x,a_i)=A+B*r`, the rational reflection coefficients are

    2*T(x,a_i)/4=A,    2*T(x,b_i)/8=B.

Both are integers, so all eight rational reflections preserve the trace
lattice. Their normalized products are all native Coxeter products:
cyclic a/a and b/b pairs have product -1/2; corresponding a/b pairs and
cyclic cross pairs are orthogonal; the two within-layer opposite pairs
have product -1; cross pairs for (1,3) have product -sqrt2, and those for
(2,4) have product -1/sqrt2. The remaining entries are obtained by
symmetry. A real point `(q,0)` in the two-embedding decomposition, with
q in the compact original chamber, is strict against every lifted wall.
All walls are actual, for projection of that point to any one of them
is strict against the others by their nonpositive products.

Nevertheless the trace chamber has infinite volume. Equation (2) gives
the nonzero trace-null vector

    z=a1+a3=b2+b4,   T(z,z)=0.

Its products with `(a1,a2,a3,a4;b1,b2,b3,b4)` are exactly

    (0,-4,0,-4; -8,0,-8,0).                                (3)

It is therefore an actual future ideal point, with precisely four
incident walls a1,a3,b2,b4. Their Gram matrix is two mutually orthogonal
affine A1 blocks. In a four-dimensional horosphere, those four walls
bound a two-dimensional rectangle times R2, rather than a finite-volume
cusp section. The other four walls are strictly negative at z and impose
no restriction on sufficiently high horospheres. Consequently the
hyperbolic five-dimensional trace chamber has infinite volume.

This is an obstruction even with integral split reflections and an actual
Coxeter wall system. It is not just a failure to prove integrality, nor
an inference from a set of candidate Gram vectors. The rational trace
reflection subgroup is discrete, but its finite covolume does not follow
from compactness of the original K-reflection chamber. Additional mirrors
of the full trace lattice, if any, are a separate problem and are not
silently supplied by restriction of scalars.

## 3. Consequence for the new uniform research lane

These two examples block direct shortcuts to a uniform compact arithmetic
bound near thirteen: nonright projection loses the sharp height descent,
and trace-form reflection splitting loses finite covolume. A successful
argument must add genuinely global contact, field or reflection-completion
information. No such completion theorem, large right-angled face theorem,
new compact high-dimensional example or improved uniform ceiling has been
proved in this note. The current uniform compact arithmetic bound remains
twenty; the user-deferred field-specific nineteen casework is untouched.
