# Terminal triples at whole B and D vertices

20 October 2026. Paper reasoning only; no scientific process was run.
This is a local necessary condition, not a dimension ceiling.

## Statement and conventions

Let p be an actual vertex of a compact arithmetic hyperbolic Coxeter
polytope whose WHOLE finite diagram is B_n or D_n, n >= 14. No actual
outside wall can close all three directions specified below. In fact the
argument does not require p to be farthest from an interior point.

For B_n, use the path 1,...,n with the sole 4 bond between n-1 and n.
The triple is {n-2,n-1,n}. For D_n, use the simple path 1,...,n-2
and the two whiskers b=n-1 and e=n attached to a=n-2. The triple is
{a,b,e}. Normals are unit outward normals.

The inputs are the finite spherical Coxeter classification, the exclusion
of an actual affine subdiagram in a compact Coxeter chamber, and
Vinberg's conjugate-positive Gram criterion for arithmetic reflection
groups. No list of number fields is used. An affine subset is excluded
even with outside contacts: its positive null combination is a nonzero
future null vector satisfying every actual wall inequality, which would
give an ideal point of the compact chamber.

If c closes direction i, replacing old wall i by c gives a positive
definite rank-n Gram matrix at the actual adjacent vertex. Conversely,
the actual projection lemma identifies precisely those replacements
with finite edge intersections. Below we use only the necessary
positive-definite replacements and the actual all-wall acuteness.

## B_n: complete terminal-triple inventory

Assume c closes n-2,n-1,n. Delete n. The retained old diagram is
A_(n-1), of rank at least 13. Its finite one-root continuations allow
only these contacts with that old path:

* no contact;
* one endpoint, label 3 or 4;
* one node next to an endpoint, label 3.

There cannot be two contacts, since they create a nonpositive cycle.
If there is a contact the connected finite continuation has rank n,
so it is A_n, B_n or D_n; exceptional finite diagrams have smaller rank.
Write beta=-c.n_n, which is the only remaining old contact.

Delete n-2. The retained old tail {n-1,n} has label 4. Unless c also
contacts n-1, its three-node path with c forces beta to be either 0
or cos(pi/3): beta=cos(pi/4) is affine and larger beta is indefinite.
If c contacts n-1 and beta>0, the retained three-node cycle is already
nonpositive. Thus every case is covered as follows.

1. No path contact: beta=0 gives a finite full old-plus-c diagram;
   beta=1/2 gives the affine F4 path on
   {n-3,n-2,n-1,n,c}, with bonds 3,3,4,3.
2. Contact at old endpoint 1: with beta=0, label 3 gives finite
   B_(n+1), while label 4 gives an affine path with 4 at both ends.
   With beta=1/2 the same affine F4 subset as in case 1 remains exact.
3. Contact at old node 2, label 3: beta=0 gives the affine fork with
   two simple leaves at node 2 and a 4 at the other end. A positive
   kernel has coefficients 1 on old node 1 and c, 2 on old nodes
   2,...,n-1, and sqrt(2) on old node n. The case beta=1/2 again
   contains the affine F4 subset of case 1.
4. Contact at old endpoint n-1: beta=0 and label 3 give the affine
   four-node star of labels 3,3,4 centered at n-1; label 4 gives the
   affine two-4 path {c,n-1,n}. Positive beta was already excluded by
   the retained cycle.
5. Contact at old node n-2, label 3: beta=0 gives the affine five-node
   fork on {n-3,n-2,n-1,n,c}, with two simple leaves at n-2 and
   terminal 4 on its other arm. Beta=1/2 is the sole remaining case.

Consequently the ONLY candidate has

    -c.n_(n-2) = -c.n_n = 1/2,
    c.n_j = 0 for all other old j.                    (B)

The old inverse Gram has C_ij=2min(i,j) for i,j<n,
C_in=sqrt(2)i and C_nn=n. Its projection cost is

    q_B = (3n-4+2sqrt(2)(n-2))/4.

The 4-cycle in (B) puts sqrt(2) in the cyclic field. Choose an embedding
of the arithmetic ground field that sends sqrt(2) to -sqrt(2). The old
B_n Gram remains positive definite, by changing the sign of its last
root. Conjugate positivity of the full old-plus-c Gram requires

    q_B' = (3n-4-2sqrt(2)(n-2))/4 < 1.

But q_B'>1 exactly when n>8+4sqrt(2), and n>=14 satisfies this.
Thus (B), and hence the terminal triple, is impossible arithmetically.

## D_n: complete fork-triple inventory

Assume c closes a,b,e. Delete b. The old diagram is A_(n-1), with
path endpoints 1,e and next-to-end nodes 2,a. The same finite
continuation classification applies. Put beta=-c.n_b.

If c has its path contact at 1 or 2, positive beta makes a cycle in
the finite replacement obtained by deleting e. For beta=0, the contact
at 1 with label 3 gives finite D_(n+1), the contact at 1 with label 4
gives an affine B fork, and the contact at 2 with label 3 gives the
affine two-fork D diagram.

If c contacts a, positive beta makes a nonpositive retained triangle
after deleting e. For beta=0, the four simple leaves at a give the
affine D4 star. If c has no path contact, the finite replacement after
deleting e permits beta only 0,1/2,1/sqrt(2). Beta=0 is a finite full
extension. Beta=1/2 extends one whisker and contains affine E8 with
arms 1,2,5; beta=1/sqrt(2) contains the affine five-node fork with
terminal 4 described above.

Finally suppose c contacts e, with label 3 or 4. Deleting e forces
beta to be 0,1/2 or 1/sqrt(2). Beta=0 gives the same affine E8 or
terminal-4 fork. Deleting a forces the two remaining contact squares
to have sum less than 1. Thus two 4 bonds are excluded; two 3 bonds
give the affine four-cycle. The ONLY remaining candidate, up to
interchanging b,e, is

    -c.n_b=1/2,  -c.n_e=1/sqrt(2),
    c.n_j=0 for all other old j.                     (D)

Here C_bb=C_ee=n/2 and C_be=(n-2)/2, so

    q_D = (3n+2sqrt(2)(n-2))/8,
    q_D' = (3n-2sqrt(2)(n-2))/8.

The cycle again puts sqrt(2) in the cyclic field. Conjugate positivity
requires q_D'<1, while q_D'>1 exactly when n>8+4sqrt(2). This excludes
(D) for n>=14 and proves the D_n statement.

## Exact limitation: two terminal B directions can survive locally

Do not replace either triple conclusion by an assertion about every
terminal pair. In B_n take c to contact n-1 and n, both with label 3,
and no other old wall. Then

    q=(3n-2+2sqrt(2)(n-1))/4,
    q'=(3n-2-2sqrt(2)(n-1))/4.

For 14<=n<=18, q>1 and 0<q'<1. With every root normalized to square
2, all simple Gram entries are -1, the old 4-bond entry is -sqrt(2),
and all entries lie in Z[sqrt(2)]. They give integral reflections and
an admissible full-rank Lorentz star over Q(sqrt(2)). Deleting n or
n-1 gives respectively A_n or A_(n-2)+A2. Every other deletion retains
the nonpositive 4,3,3 triangle, so these are exactly its two finite
replacement directions.

This is a genuine integral local Gram configuration; it is NOT a
compact-polytope construction, a finite-volume certificate, or a proof
that such a pair occurs at a farthest vertex of a compact chamber.
The new triple exclusion alone therefore does not close whole B_n or
D_n farthest energy, and does not establish a new uniform ceiling.

## Addendum: the isolated dotted arithmetic envelope

For a whole A_n, B_n or D_n old vertex, let d_i be its UNIT-normal
inverse-Gram diagonal. If an actual opposite c contacts only old node
i, with dotted product -beta, put X=4beta^2. Cyclic integrality gives
X in the ground-field integer ring. At every nonphysical embedding,
old-Gram conjugate positivity gives

    0 < X_sigma < 4/d_i.

The d_i are rational for these three finite types. If D is the ground
field degree, the positive nonzero integer norm consequently gives

    X_0 > (d_i/4)^(D-1),
    cosh L_i = sqrt(d_i X_0/4) > (d_i/4)^(D/2).

Here D>=2: the rational arithmetic compact case in these dimensions
is excluded by the usual rational-form isotropy theorem (Meyer).
Thus, for d_i>=4,

    cosh L_i > d_i/4,
    tanh^2(L_i/2) > (d_i-4)/(d_i+4).

The dotted condition also independently gives cosh L_i>sqrt(d_i).
When d_i>=4 and the degree is at least three, the displayed norm bound
strengthens to cosh L_i>(d_i/4)^(3/2). These are individual lower envelopes ONLY for
isolated dotted actual opposites. They do not cover finite releases,
white attachments, or complete a whole-component energy contradiction.
