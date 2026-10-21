# Independent review: the complete D7 edge and the discriminant-59 period

20 October 2026. **Full paper and combined checker source review PASS.**
No scientific source was imported or executed; no generation, replay,
search, Lean or TeX job ran in this review lane. Generation and separate
replay of the final combined checker remain root-scheduled work.

Reviewed paper sources:

* `results/compact_h9_golden59_d7_edge_interval_20261020.md`, SHA256
  `e4715a78f22628bf7a27051ae6588adae7803f191a2a08f476d0a02fa2223f52`;
* `results/compact_h9_golden59_long_facet_period_20261020.md`, SHA256
  `699e9922fc0cc3ca6f16d5619e02ab345fbcb39ed25ddb355d32a4120ef14f9d`;
* the inherited route, SHA256
  `725d8775cebec51d8c18b4b0a0d3be0255e8b1472621adb4d3cefa656adb429e`,
  and the literal calibrated Bugaenko root/vertex data.
* `scripts/check_compact_h9_golden59_d7_first_edge_20261020.py`, SHA256
  `5022037cc17357e1cbe745a0fc9be40b445372ad266d00bdf7bc2236dc0d70a9`;
* the initial, execution-pending edition of
  `results/compact_h9_golden59_nonreflectivity_conclusion_20261020.md`,
  SHA256 `0d54c15665f0911f426e201cf7e04e6f6ff793664a84e3bff5565301c603d0da`.
  Its mathematical composition is reviewed here; a later execution-status
  update is not being anticipated as already verified.

## Complete edge inventory and the genuine upper

Removing native E8 node 0 leaves D7 with arms 1,1,4. Its one-root finite
extensions are D8 at node 7 and the two E8 placements at nodes 1 and 2.
Joining the retained isolated c as well allows only D9, through node 7
and two simple bonds. Otherwise a short new root is D7-orthogonal and
has with c label 2,3 or5. A long new root meeting the finite corner must
be perpendicular to every short white by the previously reviewed
whole-arrangement isolation theorem. This gives exactly the eight rows,
not a selected list of convenient attachments.

I independently reconstructed the E8 inverse entries: omega0 has square
4, its products with omega7,omega1,omega2 are 2,5,7, and the three
residual squares are 1,7/4,7/4. They give precisely
R=8,6,6-2phi,8b,4,1,1,2 in Z^2-2delta X^2=R. The binary difference
X=x-y is integral. For the four surviving c-orthogonal rows, the unit
h_y forces X in 2O. E8 dual integrality gives the stated Z residues, and
the long reflection condition gives Z/(4b) in O. The four nonsquares
8,6,56,2 modulo59 eliminate the other rows at every height.

I independently found and checked the smaller upper as well. With
alpha=38+61phi and beta=133+213phi=b(9+14phi), the vector
v=alpha h+beta omega0 has norm2b. Its covector divided by b is integral;
the coefficient norms 41 and649 are coprime, confirming primitivity.
The equation Zv^2-2delta Xv^2=8b puts its collision strictly before the
light endpoint. The upper does not assume compatibility with an old cap
list, nor assume itself first. The larger Pell upper and its exact
factorization agree with this smaller bisector; the larger reflected
E8 root cannot replace a first-contact proof.

## Independent finite-list arithmetic

The collision bound and positive conjugate equation give the displayed
physical and conjugate intervals. For D8, X=2Y and Z=2W reduce the
equation to W^2=1+2delta Y^2. For E8 the same equation holds with Y=X
and W=Z. Thus enlarging both E8 domains to the same D8 Y domain is safe.
In that common domain 0<Y<25 and |Y'|<1/2, so the phi coefficient d is
an integer between0 and12. The four integer constants around the exact
floor m=(isqrt(5d^2)-d)//2 cover every possibility, giving52 tests.

I independently obtained the ten Y values in the report and all their
norms. One convenient independent expansion, for Y=p+qphi, is

    N(1+2delta Y^2)
      =1+44(p^2+q^2)+82(2pq+q^2)-236(p^2+pq-q^2)^2.

It gives exactly 99,641,2059,5771,5401,5571,15489,1531,20249,24131.
None is an integer square. This excludes D8 and both E8 placements
without assuming any unproved relative-unit or coefficient-height bound.

For the long row the integral parameter t=Z/(4b) satisfies
A^2=(-4+14phi)t^2+1-phi. The physical and conjugate bounds give
0<t<=9+14phi and |t'|<9/25. Hence0<=d<=16 and the same four constants
give68 tests. I independently reconstructed its ten t values. For
t=p+qphi, the radicand norm is

    -236(p^2+pq-q^2)^2
      +38(p^2+q^2)+44(2pq+q^2)-1.

It gives 305,1181,3475,1895,1189,9481,11071,12659,25205,1681. Only
1681 is square. The last radicand is exactly (38+61phi)^2; positive
contact chooses its positive square root and reconstructs v. I also
checked a separate, tighter enumeration in A=X/2 using its modulo59
condition; it yields the same unique root and is not needed for the
author's complete ten-t argument.

Consequently the full bounded interval has just v. There is no
simultaneous-contact dominance choice. A first contact exists by the
genuine interior upper and local finiteness of the integral reflection
arrangement. Its finite spherical inventory is exhaustive, so it is the
actual first simple wall. The endpoint has the full positive-definite
germ D7+c+v. No finite-volume assumption on the full chamber enters this
argument.

## Full integral period and global implication

I rederived every displayed G column in both bases. The basis p,c,e_i
is a full integral basis, so the integrality claim hides no discriminant
gluing. A short independent Gram check uses the orthogonal frames in
the three-space perpendicular to D7:

    (r,a,P) -> (v,c,H),       Gram=diag(2b,2,-4/phi).

Here P=2p+phi omega0, H=(10+18phi)h+alpha omega0, and
delta(10+18phi)=2beta. These identities verify the frame images,
H^2=-4/phi and orthogonality. Since D7 is fixed pointwise, the full
Gram is preserved. Integral columns and det(G)^2=1 make the inverse
integral; the future apices give the correct time component. Its trace
is exactly40+55phi, greater than10, so it has infinite order.

The old finite germ is genuinely the archived vertex
[0,2,3,4,5,6,7,9] on the actual r facet, with r added. Its nine outward
normals are exactly D7+a+r. The complete new-edge result gives the new
nine-wall germ D7+c+v. G sends the old outward normals and future apex
to the new ones. An integral lattice isometry permutes the entire
integral-mirror arrangement; the two full inward germs therefore show
G(C)=C, not merely equality of two prefix corner Grams.

Finally the full reflection group W has trivial intersection with its
chamber stabilizer. The infinitely many powers of G thus give distinct
cosets modulo W. This proves infinite index in O^+(L), conditional only
on the stated inherited full-arrangement inputs, and excludes the
literal displayed integral L59 host. It does not exclude every golden
binary extension or its entire commensurability class, and gives neither
a new compact example nor a uniform dimension bound.

## Combined checker source

I read the complete 240-line worker and checked its imported pure exact
field and matrix definitions. Its field encoding, integrality and sign
tests agree with O=Z[phi]. It imports none of the earlier producer's
scientific functions. Instead it pins the earlier actual certificate and
producer source, reads the former's sixteen roots and strict center, and
uses the small calibrated Cartan data to reconstruct this batch.

The eight profile rows are rebuilt from their full normals. Their six
quadratic norm coefficients are checked directly, not by assuming an
entered R table equals a saved answer. All nonsquares modulo59 are
verified against the complete finite residue-square set. The small upper,
reflection covector and primitive coefficient norms are checked exactly.
The two integer-strip loops enumerate exactly52 and68 pairs from the
proved boxes, apply physical and conjugate inequalities by exact signs,
and compare their resulting lists and all twenty computed norms with the
paper. The E8 even subdomain is checked separately but is unnecessary
for the larger-domain nonsquare exclusion. The only square long norm is
not accepted by itself: its whole radicand is checked equal to alpha^2,
and the resulting full vector, norm, integrality and reflection matrix
are checked again. Its positive root is unique in the field.

All sixteen old-root relations are tested only after that arithmetic
selection. The endpoint has exactly the listed nine vanishing seed
walls, is feasible, is future timelike and has positive-definite corner
Gram. A new strict interior point is constructed and checked against all
seventeen roots. These checks corroborate the full-germ interface; their
finite cap list is never used to discard an arithmetic first candidate.

The fixed period is reconstructed in the original full integral basis.
The code verifies its full Gram, computes its inverse by the Gram-adjoint
formula and checks both its inverse identity and integral entries. It
checks the old archived vertex, exact frame images, fixed D7 roots,
future endpoint and exact trace greater than10. The earlier actual
certificate's pin is an inherited result, not a claim that this second
worker independently replays its first-edge proof.

Every used source, including the worker itself, is streamed into a
before/after hash map. The main routine takes no inner lock, runs no
background work, refuses certificate replacement, has a 45-second alarm,
and compares a freshly recomputed canonical JSON in replay mode. Root's
outer gate and lock remain required. Its result explicitly leaves the
universal norm inventory, finite spherical inventory, actual-edge
geometry and chamber/coset implication as written inputs. No formal
verification of those geometric arguments is claimed.

I also read the complete root closure note. It faithfully composes the
actual forty-one-vertex facet, complete first transverse edge, complete
D7 census and integral period. It correctly distinguishes an independent
mathematical source review from a fresh-process replay of the same
implementation, and marks final execution pending in the reviewed
edition. No missing geometric interface or stronger unsupported scope
was found. The written literal-host nonreflectivity argument passes;
the final worker has not been run by this reviewer.
