# D components: a common opposite on the internal long arm

20 October2026. Complete written candidate, awaiting outside-author
review. No arithmeticity, scientific execution, optimization, TeX or
Lean. Spectator components of the actual finite vertex are arbitrary.

Use unit D_r normals on the chain1,...,r-2, with fork leaves r-1,r
at r-2. All old bonds are3. Suppose

    r>=8, 4<=i<=r-4,

and one actual outside wall c is opposite to BOTH outgoing directions
i,i+1 of an actual compact Coxeter vertex. Then its complete contacts
with the WHOLE old vertex are

    c.n_i=-1/sqrt2, c.n_(i+1)=-1/2,
    every other old product0.                         (1)

The actual edge length in direction i+1 consequently satisfies

    cosh L_(i+1)=i+(i+1)/sqrt2.                       (2)

In particular its half-tangent square is the strictly increasing
function

    tanh^2(L_(i+1)/2)
      =[i+(i+1)/sqrt2-1]/[i+(i+1)/sqrt2+1].         (3)

The D8 proof already handles the first case r8,i4. This note proves
the uniform classification, not an assumed harmless-spectator transfer.
It does not give a component energy bound for every D_r or a new
ambient dimension ceiling.

## Actual replacements and the complete right-hand inventory

Delete old i. The right component is D_s on i+1,...,r, with s=r-i>=4;
the left component is A_(i-1). The actual endpoint missing i has a
positive-definite old-plus-c diagram, including every spectator.

A finite one-node extension of D_s allows no contact to it, or one
simple3 at a leaf. A3 at an internal vertex creates a second fork
or a degree-four affine star; no finite connected Coxeter diagram
has that pattern. A4 at a leaf already has Schur cost>=1, since
the inverse-D leaf entries are2 at the long leaf and s/2>=2 at the
short leaves. Larger bonds or added contacts only worsen a positive
nonpositive witness. Two contacts create a nonpositive cycle. Thus
the right contact is none,3 at long leaf i+1, or3 at one fork leaf.
The fork option is finite only in the appropriate D/E small ranks,
but the subsequent exclusion does not need that further restriction.

Delete old i+1 instead. The complete left A_i on1,...,i remains.
Positive definiteness permits at most one c contact to this chain.
Its single native contact is3, or4 at a path endpoint; i>=4 rules
out the H/F smaller-rank exceptions. All later nonpositive retained
diagrams are genuinely subdiagrams of one of these actual finite
replacements. Exact affine diagrams in the full star are forbidden
because the original polytope is compact.

## Right contact at the long leaf

Suppose c contacts i+1 by3. If it has no contact to i, select the
full old path from i to the fork center and add c and the two old
fork leaves. Its two forks have leaves(i,c) and(r-1,r); this is an
exact affine D diagram. Additional left/spectator contacts lie outside
this principal subdiagram. Hence c must contact i.

At the endpoint of the retained A_i that contact can be3 or4. A3
makes the exact affine triangle c,i,i+1, so it must be4. No second
left contact is possible. For any spectator contact x, the retained
path x,c,i,i-1,i-2 has bonds at least3,4,3,3, hence is nonpositive
by the affine F4 witness. Therefore no spectator contact occurs.
The right-hand inventory already excludes further right contacts.
This gives precisely(1).

## Right contact at a fork leaf

Suppose c has a3 contact at b=r-1; the other fork is symmetric.
A positive left contact must be3. An internal4 or larger gives the
retained4,3,3 star. A4 at left endpoint i gives the retained F4 path
b,c,i,i-1,i-2; at endpoint1 use old3,old2,old1,c,b. Larger bonds only
worsen the same positive quadratic witness.

A left3 at j creates an exact all3 cycle through c,b,r-2 and the old
long path to j, forbidden in a compact polytope. With no left contact,
the full old D_r plus c contains the exact affine E8 tree with arms
1,2,5: use the other fork, b-c, and five old long-arm nodes. The latter
exist since r>=8. Spectator contacts cannot alter that principal
diagram. Thus this entire right-fork case is impossible.

## No right contact

An internal left3 at j>=2 creates an exact two-fork affine D diagram:
its first fork has leaves c,j-1, and its second has the old short
leaves, joined by the old path to r-2. An end4 at i gives the full
4,3,3 star, and an end4 at1 gives the exact terminal-4 affine D
extension along the old long arm. These are forbidden independently
of spectators. The remaining possibility is a3 at1, possibly with
spectator contacts.

Choose a minimal nonpositive principal subdiagram of the old vertex
and c. Such a diagram exists because all n+1 normals cannot be
positive definite in signature(n,1). It contains c and i+1, since
the actual replacement deleting i+1 is positive definite. It is
connected and is Lanner, as a compact polytope has no affine principal
subdiagram. Its unique old path from c to i+1 is c,1,...,i+1, with
i+2>=6 vertices. Spectators attach to this D component only through c,
so do not shorten that path. This contradicts the explicitly imported
classical Lanner order bound of5.

If there is no left contact either, the D-component derivative of c
is zero, so it cannot close direction i+1. This exhausts the cases.

## Exact actual length

For i+1<=r-3 the unit inverse-D entries are

    C_ii=2i, C_(i+1,i+1)=2(i+1), C_(i,i+1)=2i.

At the actual endpoint missing i+1, the left A_i contact4 has white
Schur cost i/(i+1), leaving residual delta=1/(i+1). The outgoing
derivative is

    [sqrt2*i+(i+1)]/sqrt(2(i+1)).

Dividing by sqrt(delta) gives(2). Formula(3) is the actual hyperbolic
half-tangent identity, not the full tangent used in the distinct-wall
shadow law. It increases strictly with i. For instance i4 recovers
4+5/sqrt2>7 and half-tangent square>3/4.

## Inputs and intended use

Finite spherical classification and Lanner order<=5 are explicit
classical inputs. Actual endpoint/covering conventions are those of
`compact_universal_e6_e7_e8_component_energy_20261020.md`, SHA256
`acf261d97f6eb1a11a64502e17e079ef7f418b9746fbf85dbdf6915750fe9a18`.
The earlier center/neighbor theorem
`compact_dn_center_neighbor_common_opposite_energy_20261020.md`, SHA256
`e371611c4903ad362654ee0d4656457140cdc5e4374b4e19a2697718c2d679b1`,
has different center/fork cases and is not silently replaced by(1).

If a proposed component energy upper lies below(3), the two actual
opposites at i,i+1 must be distinct. This supplies an exact consecutive
full-tangent constraint for multi-direction arguments. Center/fork
exceptions and larger-rank energy estimates still require separate
proofs. No rank/14 rate, field-wide lattice exclusion or formal
verification is claimed solely from this lemma.
