# Bugaenko facets: exact binary reduction, infinite freedom and a uniform E8 collar

20 October 2026. **Paper-only structural audit, not a new dimension bound.**
No worker, root search, matrix program, Lean or TeX job ran. The closed
literal hosts remain closed; no new isolated host is opened here.

The main new quantitative statement is Section5: every mixed integral
mirror in an even golden lattice E8+B stays at least `asinh(sqrt(phi))`
from the whole E8-orthogonal hyperbolic subspace. Section6 links the
binary edge length to the normalized discriminant even when other root
norm classes intervene. Neither statement supplies a uniform volume
upper bound or proves that every possible binary extension is excluded.

## 1. The facet forces its field

Let P be a compact arithmetic first-type Coxeter polytope in H9 with a
metric facet isometric to Bugaenko's compact eleven-facet H8 chamber F.
The reviewed metric-facet rigidity lemma applies: every F node has an
incident intrinsic3- or5-bond, so all eleven adjacent parent walls are
perpendicular to the base wall. Their actual normals therefore have the
intrinsic F Gram, with no projection correction.

Write the admissible quadratic form over its totally real field K. The
5-bonds give `cos²(pi/5)=phi²/4 in K`, hence Q(phi) is a subfield.
The connected3/5 graph makes every pair of normal squared lengths have
square ratio in K: on an edge `(v,w)²=v²*w²*cos²(theta)`, with cos(theta)
in Q(phi). Rescale the eleven normals over K to one common norm2lambda.
Their rank-nine Gram is lambda times the literal Bugaenko Gram. At any
other embedding of K fixing phi that Gram would still be indefinite,
contradicting admissibility away from the physical embedding. A proper
totally real extension has more than one such embedding. Consequently

    K=Q(phi).

Scale the whole quadratic form by1/lambda; its eleven normals now have
norm2. This is the usual arithmetic-face field principle proved here
for the specific connected facet; no field list for higher dimensions
has been used.

## 2. No hidden evenness, saturation or overlattice assumption is needed

Let Gamma be the reflection group of P. For each gamma in Gamma some
power lies in an integral arithmetic group commensurable with Gamma.
Thus its eigenvalues, and hence its trace in K, are algebraic integers.
For any two norm-two vectors v,w in the Gamma-orbit of the facet normals,

    trace(s_v*s_w)=6+(v,w)².

It follows that `(v,w)² in O`; since O is integrally closed, `(v,w) in O`.
The O-span R of all these orbit vectors is therefore even and integral.
It is a finite O-module: arithmetic commensurability and a finite coset
decomposition place the whole orbit in one fixed bounded-denominator
lattice. Its rank is ten. Otherwise Gamma would preserve the Lorentzian
rank-nine span of the facet; distance to its hyperplane would be an
unbounded Gamma-invariant function, impossible on a compact quotient.

The literal eleven facet vectors span `M0=<-2phi>+E8`; its known nine-root
basis has unit determinant in that model. Since E8 is unimodular, its
orthogonal projection of every R vector is integral. Hence it splits:

    R=E8 orthogonal-sum B.

The time vector u0 of square -2phi is primitive in B: if u0=a*v, evenness
would imply a² divides the unit phi. As O=Z[phi] is a PID, extend u0 to
an O-basis of B. In that basis

    B=[[-2phi,a],[a,2c]],     a,c in O.                 (1)

Every wall reflection of P preserves R and has a primitive normal in R.
Thus (1) is a necessary reduction for **all** arithmetic P in this facet
class, not an assumption that a previously searched host is maximal.

There are just two parity/glue forms. Replacing the second basis vector
t by t+k*u0 changes a by -2phi*k. Rescaling t by a power of phi makes
every nonzero residue modulo2 equal1, since O/2O=F4 and phi generates
its multiplicative group. Therefore (1) can be normalized to a=0 or1.
The first case is the direct orthogonal extension. In the second,
`r0=u0+2phi*t` is perpendicular to u0 and spans an index-two extension
with M0. The remaining coefficient c is not bounded by this reduction.

## 3. An infinite genuine module family survives the local conditions

For u=phi^k consider

    B(u)=[[-2phi,phi*u],[phi*u,2]],
    b_k=4+phi*u²,       det B(u)=-phi*b_k.

Admissibility is exactly k>=-1: `b_k'=4-phi^(-2k-1)>0` holds precisely
in that range. If L_n denotes the Lucas trace `phi^n+phi'^n`, then

    N(b_k)=15+4L_(2k+1)=11,19,31,59,131,319,...

strictly increases for k>=-1. Thus these literal integral modules are
pairwise nonisometric, not just different presentations of finitely
many discriminants.

All have the primitive normals `r=(u,2;0)`, `c=(0,-1;0)`, with
`r²=2b_k`, `c²=2`, `(r,c)=-b_k`. Together with the eleven F walls these
give an actual acute thirteen-facet Coxeter seed. Its compact base facet
is exactly F, and it has the actual finite E8 edge from `(1,0;0)` to
`p=(2,-phi*u;0)`. Products of c with the three positive-height F roots
are `-phi^(k+j+1)`, j=0,1,4: they are native3/5 bonds or ultraparallel
for every admissible k. A strict interior and actuality follow from a
small inward perturbation of an interior point of F and acute projection.
These local facet/germ/finite-edge conditions therefore do not give a
finite list of binary modules.

These are not compact counterexamples. For k>=0 put D=phi*b_k and

    X=(1,-phi*u/2;(D/8)*omega0),  omega0²=4.

It is feasible for all thirteen walls: its F-wall evaluations are the
positive multiple D/(4phi) of those at the known F vertex
`(1,0;(phi/2)*omega0)`; c is tight and r is strict. But
`X²=D(D-8)/16>0`. Perturbation towards the strict center proves that
these finite seeds have infinite volume. The unresolved issue is their
possible full-arrangement completion, not a failure to notice this escape.

Existing exclusions are not a classification. k=-1 is the old literal
discriminant11 form; k=0 is the now-closed H19 form. k=2 is L59 via
`U=[[1,-1],[0,-1]]` from the L59 binary basis. For k=3 the explicit map

    U=[[1,-phi²],[0,phi]],
    U^t*[[-2phi,-1],[-1,6]]*U=[[-2phi,phi^4],[phi^4,2]]

identifies the already-closed norm131 module and its oriented base germ.
This uses the integral basis change, not equality of determinant norms.
The saved norm31 exclusion is over Q(sqrt2), not this k=1 golden form.
That is a bounded repository distinction, not a literature-wide novelty
claim or authorization to start a new host-by-host search.

## 4. What extends uniformly, and what does not

The intrinsic Bugaenko reversal S satisfies `S*u0-u0 in 2M0`. In (1) its
extension fixing the perpendicular base normal has

    A|M0=S,       A*t=t-a*(S*u0-u0)/(2phi).

It is an integral isometric involution for every a,c. This resolves the
module/glue question for this particular extension. It is **not** by
itself a chamber stabilizer theorem. Extra-norm R mirrors can subdivide
F; then A may permute the smaller full chambers. Prime or suitable
squarefree discriminant isolation supplies the missing full-facet premise
in the closed examples, but arbitrary c needs its own norm/isolation proof.

The simplest second-involution mechanism is also not uniform. For a
positive-height F wall g_j of time coordinate phi^j, j=0,1,4, the point
inversion at c+g_j has denominator `phi*u*phi^j-2`. Its spatial covector
has a unit entry and c+g_j has a unit y coordinate, so this inversion is
integral exactly when that denominator is a unit. It is a hyperbolic point
inversion only when n=k+j+1>=2. But

    N(phi^n-2)=(-1)^n-2L_n+4

is +/-1 in that range only for n=2. Thus the immediate known-wall
mechanisms cover k=0,j=1 and k=1,j=0, and do not extend to all large k.
Indeed L2=3, L3=4 and the positive Lucas sequence then increases; for
odd n the norm is 3-2L_n, and for even n it is 5-2L_n.
The special L59 period is genuinely additional geometry, not a uniform
consequence of having the base-facet involution.

## 5. A uniform quantitative collar for every even golden E8 splitting

This lemma works more generally for `R=E8+B`, with B admissible of any
Lorentzian rank at least two. Let m be a primitive integral reflector,
`m²=2d`, d in O totally positive. Reflection integrality and unimodularity
give `m_E8=d*z`, z in E8(O). Therefore

    T=m_E8²/m²=d*z²/2 in O,      0<=T'<=1.            (2)

The upper bound follows from positive definiteness of the conjugate B
block. If m has nonzero E8 part then T,T'>0. The integer norm forces
T>=1. If T=1, its B part vanishes by conjugate positive definiteness.
Otherwise write T=A+Bphi. Since T>1>=T'>0, B is a positive integer and
A>=1. Consequently

    T>=1+phi=phi²,       m_B²<=-phi*m².              (3)

Thus every root is either pure E8, pure B, or has timelike B projection
satisfying (3). Pure E8 mirrors give the ordinary E8 zero arrangement;
pure B mirrors restrict to exactly the integral reflection arrangement
of B. No mixed mirror meets the hyperbolic subspace H(B).

Orient a mixed wall towards a chamber touching H(B). Its B projection
is future timelike, so reverse Lorentz Cauchy gives, for every future
unit x in H(B), `-(x,m)>=sqrt(phi*m²)`. Hence the entire mixed mirror
stays at distance at least

    rho0=asinh(sqrt(phi))

from H(B), uniformly in the determinant and root norm class.

The chamber contains the radius-rho0 normal tube sector over its actual
B chamber, with angular sector one E8 Weyl chamber. In H9, if its binary
interval has length ell, the Fermi-coordinate volume element is
`cosh(rho)*sinh(rho)^7 d rho dt d Omega`. Therefore

    Vol(C)>=ell * [phi^4/8] * Vol(S7)/|W(E8)|.        (4)

This is a genuine volume lower bound, not an inference from a bounded
root prefix. It supplies a geometric collar where arbitrary nonright
projection heuristics had no such control.

## 6. Discriminant growth cannot be hidden by an intervening norm class

In the normalized odd case a=1 put `D=1+4phi*c`. The base root is
`r0=(1,2phi)`, with norm2phi*D and covector `(0,D)`. Let v=(x,y) be the
next actual binary wall, of norm2d, and let ell be their axis distance.
Set `q=4cosh²(ell)`. Exact two-dimensional Gram identities give

    q=D*y²/(phi*d),   q-4=(y-2phi*x)²/(phi*d).

Reflection integrality says `(y-2phi*x)/d in O`. Thus

    q(q-4)=D*[y(y-2phi*x)/(phi*d)]².                 (5)

The binary reflection matrices are integral, so q is integral (also
their product trace is q-2). Conjugate positive definiteness gives
`0<q'<4`. Taking the norm of (5), using the nonzero integral square
factor and `q'(4-q')<=4`, proves

    |N(D)|<=4q(q-4).                                (6)

No squarefree, prime norm or selected first-wall norm is needed. In the
direct case `B=diag(-2phi,2c)` the identical calculation gives
`q(q-4)=phi*c*(4xy/d)²`; since `2x/d in O`, it yields the uniform weaker
bound `|N(det B)|<=64q(q-4)`. Thus a large normalized determinant forces
a long **actual** E8 edge even if c from Section3 is not the first wall.
Together (4)--(6) force chamber volume to tend to infinity with determinant.

This still is not an asymptotic nonreflectivity theorem: a uniform upper
bound for these full integral reflection chambers has not been proved or
imported. A bound stated only for maximal arithmetic reflection groups
cannot be applied silently to arbitrary finite-index reflection subgroups.
Nor does a growing volume lower bound alone prohibit compactness.

## Inputs and remaining gap

Metric-facet rigidity is the reviewed source
`coxeter_metric_facet_angle_rigidity_20261020.md`, SHA256
`f688e520f3ff7ac43225aaacd89413c03f5dd3d2d008d7e4dcdcc48326ac83fb`.
The literal full Bugaenko chamber, reversal, native coordinates and current
H19 closure inputs are pinned in
`compact_h9_discriminant19_d7_edge_period_20261020.md`, SHA256
`ecc61a2dcff6df1678cc94533aa201c70cd0b262e26deab5ef85f5d64ebaf824`.

The audit reduces the entire metric-Bugaenko-facet class to two explicit
binary glue types and supplies a quantitative collar/discriminant bound.
It does not force such a facet in an arbitrary compact arithmetic polytope,
bound the remaining c values by compactness, control every additional norm
class on F, or create a second chamber symmetry in all cases. Those are
the precise missing global inputs; the finite list of closed hosts does
not fill them. No further host census is proposed here.
