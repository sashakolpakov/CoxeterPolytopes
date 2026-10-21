# Independent review of the quadratic13 far edge and chamber period

20 October 2026. **PASS: no mathematical flaw found.** This is a full
independent written-proof review. No scientific program, import, exact
worker, Lean compilation or TeX build was run. Metadata hashing and source
reading do not constitute an execution of the proposed checker.

Reviewed source:
`compact_h9_sqrt13_a9_far_end_and_period_20261020.md`, SHA256
`86ffab56f4696298fac4ae34b8bc8772dd1e7a6daf2e95e48d56291e810bbcf2`.
I also read the complete inherited two-edge proof
`compact_h9_sqrt13_b9_two_edge_route_20261020.md`, SHA256
`b922647902d32beccd98fbe5a2781ddcbb9cafaeb10a12902eb3ebbce809547f`,
and A8-edge proof `compact_h9_sqrt13_a8_edge_interval_20261020.md`,
SHA256 `e3e62a6ec4d834e2e0c7d4abc4d28961839f1a467ac06966ef08854757814228`.
Their actual primitive alphabet, chamber germs and first-contact arguments
are inherited inputs. This review does not claim independent execution
of their certificates or of a future far-edge certificate.

## The plane and finite upper

I independently rederived the retained plane from the eight literal
equations for `a,alpha_1,...,alpha_7`. Its integral `(x,z)` basis P0,R0
has determinant eta, and the displayed norms are respectively -eta and9,
with zero cross product. Thus the plane is integrally split as stated;
there is no discriminant-gluing assumption hidden in this step. The
vector P=P0/eta is integral and has unit norm -1/eta.

The longest A8 Weyl element is integral and fixes this plane. The point
inversion J_P is integral because its rank-one coefficient is -2eta,
and is future-preserving because it fixes a future timelike vector.
It commutes with the longest element. Their product H is therefore an
integral future-preserving involution before any chamber claim is made.
I obtain the stated b0 and b=H alpha_8, including b's norm2 and its
products -1 with a, zero with alpha_1,...,alpha_7, -2eta-1 with alpha_8,
-eta squared with alpha_9, -15eta-5 with v and -40eta-11 with d.
These give the finite upper; they are not filters in the census.

The two expressions Q=(eta+3)P0-eta R0 and
Q'=(eta+3)P0+eta R0 both have norm -10eta and positive time coordinate.
They are on opposite sides of the timelike P0 ray on the same released
edge line. The construction does not require P or Q' to be in C before
the first-contact proof.

## Complete profiles, integrality and orientation

Finite rank-nine classification containing A8 leaves exactly the eight
literal placements printed in the source: isolated norms1,2; short B9
at either end; long D9 next to either end; and long A9 at either end.
There is no exceptional rank-nine continuation. This covers the actual
new simple wall, including any simultaneous finite-root enlargement;
it does not presume that every mirror through the endpoint is simple.

For an attached placement I independently obtain

    A=2omega*x-(eta+4)z-k-eta+5,
    Z=9[(eta+2)z-eta*x-1]+(eta+3)k.

The two-variable coefficient determinant for A and the bracketed linear
form is 1/eta, again a unit. Consequently A integral and
Z congruent to (eta+3)k modulo9 are both necessary and sufficient for
the displayed native coordinates to be integral. The source's
reconstruction formula agrees with this direct calculation. Orthogonal
projection gives norm equation Z squared-9eta A squared=R, with
R=9q-k(9-k); no linear or constant correction was dropped.

The isolated case is separate and correctly uses its stronger native
equation 9B squared-eta A squared=q. At the two indicated split primes
above3, q=1 and q=2 respectively force a square equal to2 in F3.
Both isolated norms are therefore excluded at all heights.

The sign step uses BOTH endpoints. Strict feasibility at Q and contact
by Q' imply

    Z>0,   -(Z/(eta+3)) < A <= Z/(eta+3).

Thus Z>=(eta+3)|A|, with equality on the starting side excluded, and
10A squared<=R follows. Conjugating the binary equation gives
A' squared<=R eta/9. As R<=10, these imply |A|<=1 and |A'|<2;
the omega coefficient is an integer of absolute value less than
3/sqrt13<1. Therefore A=0 for R=1,4 and A=-1,0,1 for R=10.
This is a complete interval, without a height prefix or a cap test.

The congruences exclude R=1,4 by their nonzero eta coefficient modulo9.
For R=10,A=0, sqrt10 is not in Q(sqrt13). For A=+/-1, the positive
square root is eta+3 and only k=1 reconstructs integrally. I checked
that A=-1 reconstructs exactly -b0; its product with Q is zero, so it
is excluded by the strict outgoing condition. The surviving A=1 root
is b. This distinction matters: the former is a nonsimple starting
mirror, not an earlier outgoing wall.

Local finiteness, the known actual inward A9 germ and the finite upper
then prove that b is the first simple wall of the full arrangement.
There is no unresolved simultaneous candidate. The new full germ has
chain beta_8,...,beta_1,b, and H maps the old germ to precisely this
germ. Arrangement-cell uniqueness consequently gives H(C)=C.

## Independent trace and global implication

I rederived the four pairwise orthogonal norm-two gamma roots in the
source. Their reflection product is w0, with
w0=I-sum gamma_j tensor gamma_j flat. Directly from the native G columns
I obtain their four products

    -63eta-19,  1-2eta,  -eta,  2.

Their sum is -66eta-16. Independently P dot GP=-51eta-14 and trace G=6.
Thus

    trace(HG)=-6+(-66eta-16)-2eta*(-51eta-14)
             =268eta+80>10.

This proves infinite order in a ten-dimensional real representation;
merely composing the two involutions would not prove it. Both preserve
the same full chamber, so every power preserves C. The standard free
action of a reflection group on its open chambers makes these powers
distinct cosets modulo the full integral reflection group W.

The final inference uses the admissible arithmetic lattice O^+(L): an
infinite-index subgroup of this lattice cannot itself be a lattice.
Equivalently C has infinite volume. A compact Coxeter reflection group
whose reflections all preserve this literal L would lie in W; its
discrete overgroup W would then be a lattice, a contradiction.

Accordingly the written proof excludes this literal integral host. It
does not exclude its similarity or commensurability class, change a
uniform dimension bound, or supply a new compact example. The finite
spherical classification, admissible-arrangement local finiteness,
actual-germ continuation and arithmetic-lattice/index facts remain
explicit mathematical inputs. A scalar/matrix worker may check the
displayed identities but does not by itself establish those inputs.
