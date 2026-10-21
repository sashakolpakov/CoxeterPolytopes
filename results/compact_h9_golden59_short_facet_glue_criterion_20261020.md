# Discriminant59: automatic extension at the short c facet

20 October 2026. New paper lemma, proposed for independent review. No
scientific computation, worker, census, optimizer, Lean or build ran.
This does not settle reflectivity of the literal host or construct H9.

Use the host and actual short wall c from
`compact_h9_golden_discriminant59_edge_route_20261020.md`:

    O=Z[phi], phi^2=phi+1, delta=5+12phi,
    L=E8(O) + [[-2phi,-1],[-1,4+2phi]],
    c=(1,1;0), q(c)=2,
    h=(3+2phi,1+2phi;0), q(h)=-2delta, (h,c)=0.

The proof below uses these literal identities, not a similarity of
quadratic spaces or a claim that every short facet is Coxeter.

## 1. Exact orthogonal lattice and half-glue

The c covector on the binary coordinates is

    (-(1+2phi), 3+2phi)=(-phi^3,3+2phi).

Since phi^3 is a unit, all integral vectors perpendicular to c have
binary part an O-multiple of h. Thus

    Mc=L intersect c-perp = E8(O) orthogonal-sum O*h.

The binary basis (c,h) has determinant2. With p=(1,0;0),

    p=(h-phi^3*c)/2,
    L=(Mc orthogonal-sum O*c) + O*(h-phi^3*c)/2.       (1)

Indeed c and p already form a unimodular binary basis. Equation(1)
therefore reconstructs the full lattice, including its index-two-ideal
glue, rather than only its rational quadratic space.

## 2. Every Mc isometry fixes the glue

**Lemma.** Every integral isometry T of Mc extends uniquely over K to
an integral isometry of L fixing c.

**Proof.** The prime2 is inert in O, so O/(2)=F4. The E8 Gram is
unimodular and its reduction modulo2 is nondegenerate. The h block
reduces to zero. Consequently the radical of the bilinear form on
Mc/2Mc is precisely the one-dimensional h-line. Any integral isometry
sends

    T(h)=a*h modulo 2Mc,

with a in F4 nonzero. The even lattice Mc also has the integral
quadratic refinement Q(x)=(x,x)/2, preserved by T. Here

    Q(h)=-delta=1 modulo2.

On the radical line Q(a*h)=a^2 Q(h). Hence a^2=1 in F4, which forces
a=1. Thus T(h)-h belongs to2Mc. Extending T by c->c preserves the
half-glue vector in(1) modulo Mc, so sends L into L. The same argument
for T inverse gives equality. QED.

Equivalently, restriction identifies the integral stabilizer of the
vector c with O(Mc). There is no additional mod2 obstruction for a
putative integral E8-corner map. This says nothing about nonorthogonal
ambient mirrors cutting c; those remain essential to its actual metric
facet.

## 3. From an actual E8+c endpoint to a genuine integral map

Suppose beta_1,...,beta_8 are outward norm-two integral roots in Mc,
forming the SAME ordered E8 Cartan Gram as the initial eight roots
alpha_i. Suppose they and c are the full inward germ of an actual
second finite vertex of the full L reflection chamber C.

The module E_beta they generate is not merely assumed primitive. Its
Gram determinant is1, so the orthogonal projector onto it is integral:

    P(x)=sum_ij beta_i*(C_E8^-1)_ij*(x,beta_j).

It follows that

    Mc=E_beta orthogonal-sum N,

where N is a rank-one free O-module. Here freeness uses the golden PID.
Choose a generator n of N and let u be the unit determinant of the
change from the original basis (alpha_1,...,alpha_8,h) to
(beta_1,...,beta_8,n). Taking Gram determinants gives

    q(n)=-2delta*u^2.

Thus h_beta=n/u is an integral primitive generator with EXACT norm
-2delta. Choose its sign to be future-pointing. The normalization is
important: a stored endpoint vector may be a nonprimitive K-multiple
of h_beta, so its raw scalar or norm must not silently be used in the
isometry.

The basis map alpha_i->beta_i, h->h_beta is therefore an integral
isometry of Mc, with integral inverse. Section2 extends it to L fixing
c. It preserves time orientation and maps the initial full inward germ
to the stipulated second actual germ. Since an integral isometry
permutes the FULL integral-mirror arrangement, it stabilizes C.

The actual/full-germ premise is indispensable. A norm equation or a
prefix-feasible E8 frame alone does not establish it. Nor does this
argument prove that the resulting chamber stabilizer has infinite order.

## 4. A same-attachment E8 endpoint is never the first wall

There is a useful consequence for a released-node edge. Suppose the
retained seven E8 roots form D7, and a proposed norm-two wall w completes
them back to E8 with EXACTLY the same attachment as the released old
root alpha. Both alpha and w are perpendicular to c. If the proposed
second E8+c vertex is different, the unique germ isometry fixing D7+c
pointwise is

    T=s_q,   q=alpha-w.

Indeed q is perpendicular to the retained roots and c, and
q^2=4-2(alpha,w). At the old endpoint,
(q,x)=-(w,x)>0; at the proposed new one, (q,x)=(alpha,x)<0. Thus q-perp
meets the relative interior of the finite edge, making the nonzero q
spacelike. Reflection in q swaps alpha and w and preserves time
orientation. Section3 therefore makes T integral on Mc, and Section2
extends it integrally to L. Consequently this is itself an integral
reflection of L, not a new external chamber symmetry. Its interior
intersection contradicts w being the first endpoint wall of the FULL
chamber. Equivalently a nonidentity element of the full reflection group
cannot stabilize its own chamber.

Hence a same-attachment E8 completion may be a useful finite upper
reflection, but cannot win a complete first-wall census. An intervening
mirror must occur. The OTHER D7 fork attachment can involve the leaf
swap and is not covered by this single-reflection obstruction; its germ
map and possible period must be considered separately.

This eliminates one tempting false period route. It is not a
nonreflectivity theorem for L59, a classification of its remaining
edges, or an exclusion of its quadratic-form commensurability class.
