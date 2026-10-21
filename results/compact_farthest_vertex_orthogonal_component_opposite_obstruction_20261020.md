# Opposite walls at a farthest vertex respect orthogonal components

20 October 2026. Paper proof only; review pending. No computation,
enumeration or formal verification. Root proposed the actual right-
triangle argument, resource_cleanup independently checked it, and
batch_guard supplied the shorter direct Schur proof used below.

**Theorem.** Let P be a compact hyperbolic Coxeter polytope and z an
interior point. Choose an actual vertex p maximizing d(z,p) over the
vertices. Two outgoing edges at p with orthogonal unit directions
cannot have the same actual opposite wall. In particular, the partition
of outgoing directions by their opposite wall refines the irreducible-
component partition of the finite diagram at p.

There is no arithmetic or one-contact hypothesis. No dimension ceiling
is asserted.

## Proof

Use Lorentz signature(n,1), future unit points of square-1 and outward
unit facet normals. The incident n normals n_i at the finite simple
vertex p have positive-definite Gram B. Its inverse C is entrywise
nonnegative. For their dual vectors omega_i put

    d_i=C_ii,       u_i=omega_i/sqrt(d_i).

Thus the inward direction of edge i is -u_i. If its length is L_i>0,
its other actual endpoint is

    p_i=cosh(L_i)p-sinh(L_i)u_i.

Write z=E(p-x), where E=-z.p>0 and x is tangent at p. The distance
comparison -z.p_i<=-z.p gives the farthest-vertex inequalities

    x_i:=x.u_i>=tau_i:=tanh(L_i/2)>0.                 (1)

Suppose edges i,j share an actual opposite wall c, and u_i.u_j=0.
Let W be all other incident walls. The standard all-wall acuteness
property gives beta_s=-c.n_s>=0 for s in W. The positive-definite
Gram B_W has nonpositive off-diagonal entries, hence B_W^-1>=0.
Put

    alpha=B_W^-1 beta_W>=0,
    c_F=c+sum_(s in W) alpha_s n_s.

Then c_F is orthogonal to every n_s in W. Moreover, because z is
STRICT interior and these are actual walls,

    z.c_F=z.c+sum alpha_s(z.n_s)<0.                    (2)

This sign uses a positive combination of actual wall inequalities;
c_F need not itself be an actual wall normal.

The orthogonal complement of the retained white normals is spanned
by p,u_i,u_j, with Gram diag(-1,1,1). Let A=-p.c>0 and
r_i=tanh L_i, r_j=tanh L_j. The actual endpoints p_i,p_j lie on c
and on all retained white walls. Consequently

    c_F=A p-(A/r_i)u_i-(A/r_j)u_j.

Substituting z=E(p-x) in(2) yields

    x_i/r_i+x_j/r_j<1.                                (3)

But(1) and r_i=2tau_i/(1+tau_i^2) imply

    x_i/r_i+x_j/r_j
      >=tau_i/r_i+tau_j/r_j
       =1+(tau_i^2+tau_j^2)/2>1,

contradicting(3). This proves the theorem for any orthogonal dual pair.
The empty W case is interpreted with no correction terms.

For distinct whole finite components, B and C are block diagonal, so
u_i.u_j=0. Within an irreducible component C is strictly positive,
which identifies these orthogonal pairs exactly. This proves the
component-partition consequence.

## Geometric interpretation and scope

Retaining W gives the actual two-face containing the two outgoing
edges. The two old side walls and the common opposite c bound a right
triangle with vertices p,p_i,p_j. Convexity verifies that this is the
whole face: its three actual vertices give the contained triangle,
while its three wall inequalities give the reverse containment. The
direct proof above establishes the farthest-point contradiction
without any additional projection or recognition premise.

An opposite wall may still close several directions WITHIN one finite
component. The theorem by itself imposes no number or length bound
there. It excludes proposed cross-component continuations once actual
finite endpoints show that one wall closes directions in both blocks,
including arbitrary multiple contacts and bridges inside the blocks.
It does not assume that every vertex can be farthest from an interior
point, or replace an actual endpoint by a formal Gram configuration.
