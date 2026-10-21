# The nine old sealed E8 quad consumers force both E8 triangle bundles

18 October 2026. **A concrete inactive-state exclusion, with D5 contacts
initially allowed.** Each of the nine E8 quadrilateral consumers in the
old component-weight barrier forces the F0 and F6 third walls at its
E8 corner. At a D5+2E8 vertex this gives both designated E8 triangle
bundles, so none of these nine consumers can occur in the state with
no designated E8 triangle. This is not yet completeness among all
possible E8 long-pair quadrilateral consumers.

The nine classes are 940,958,982,984,993,995,1107,1143,1145 in the
frozen union queue, as recorded in the
[old barrier](fresh_sqrt5_sealed_weights_20261015.md). Their actual
compact local quadrilateral geometry was checked separately in the
[geometry certificate](fresh_sqrt5_sealed_weights_20261015_geometry.json).
The new [small checker](../scripts/fresh_compact21_sealed_e8_quad_triangle_activation_20261018.py)
streams precisely their nine union objects from the 5.7-MiB queue,
with at most one 64-KiB object in memory, and writes a
[canonical mapping certificate](fresh_compact21_sealed_e8_quad_triangle_activation_20261018.json).

## Literal mapping

At the unique connected marked corner each old component is E8.
Map its branch to canonical root 2; its arms of lengths 1,2,4 map to
7, then 1,0, then 3,4,5,6, respectively. The checker reconstructs
this unique mapping directly from the literal labelled graph.
For all nine records the two outside walls are mutually orthogonal.
Each has exactly one contact with the old E8, labelled 5, and the
two contact positions are exactly canonical roots 0 and 6.

The canonical marked pairs are exactly

  {0,1,2} x {3,4,5,6} minus {(2,3),(1,3),(2,4)}.

These are the nine cross pairs with graph distance greater than two.
Thus each literal consumer has the old E8 contact pattern of F0 and
F6. It remains to check that hidden D5 spectator contacts cannot
change these walls; they are not assumed absent.

## Excluding all extra D5 contacts

The second whole E8 factor forces every new side to be orthogonal
to it by the spherical classification. Normalize the two new walls
to square two with old E8 contacts -phi, where phi=(1+sqrt5)/2.
This can be done over the ground field: if an initial field root
has square a and old-root contact b, label 5 gives a=2b^2/phi^2,
and multiplication by -phi/b gives the asserted norm and contact.

First consider the root with canonical contact 0. Its conjugate
old-E8 projection has square 4(phi')^2. Only

  2-4(phi')^2=2sqrt5-4,

strictly between zero and one, remains for any D5 projection.
At a finite corner the whole D5 plus this new root is spherical.
Its only nonorthogonal extensions are D6 and E6, attaching by one
label-3 edge at a D5 leaf. A second contact would create a cycle;
any other label would require a finite irreducible rank-six type
containing a D5 with a 4 or 5 bond, and no such type exists. The
inverse D5 diagonal at an allowed attachment is at least one
(one at its long endpoint, 5/4 at a short leaf). Hence every nonzero
D5 projection has conjugate square at least one, exceeding the
available amount. The F0-pattern wall is therefore orthogonal to
the whole D5, without any finite-height search.

Now consider the F6-pattern wall, still allowing arbitrary D5
contacts compatible with a spherical corner. With no such contacts,
its timelike squared height is H6=1+sqrt5. The F0 wall has height
H0=4+2sqrt5, and the old-E8 mutual projection is p=3+sqrt5.
These satisfy exactly H0 H6=p^2. An extra nonzero D5 projection of
the F6-pattern wall has positive square delta and changes its
height to H6+delta. Since the F0 wall is orthogonal to D5, the
mutual positive projection remains p. The outward timelike
components have the same sign at the common old vertex, so the
new walls' product is p-sqrt(H0(H6+delta)). Their literal mutual
product is zero. Thus H0 delta=0. Since H0>0 and D5 is positive
definite, delta=0 and the F6-pattern wall is orthogonal to D5 too.

This excludes every D5 spectator contact initially, rather than
just checking a guessed orthogonal lift. No particular label or
size of a free opposite product is sampled in this argument.

## Consequence and remaining scope

Both new walls are now precisely the F0 and F6 roots of the
[complete E8 triangle-domain theorem](review_compact21_d5_two_e8_triangle_coupling_20261018.md).
That theorem proves the actual triangle closures, not just spherical
corner compatibility: F0 supplies the simplex on {0,1,2} and F6
the simplex on {3,4,5,6}. Both together give the actual product
seven-face Delta^3 x Delta^4. The nine marked quadrilaterals are
its long cross-pair faces. In particular at least one, in fact both,
designated nonadjacent E8 triangle types are present at this corner.

Accordingly these nine old consumer records have no realization at
a D5+2E8 vertex in the E8-triangle-inactive state. Their relevance
to the old dual is concrete, but limited: deleting their inactive
state copies is legitimate, whereas deleting every E8 long-pair
bonus is not yet justified. Other finite-corner patterns, other
minimal cores and their spectator-contact domains still need full
coverage. No successor linear program or dimension bound is claimed.

Generation passed in chunk `a9bc2c`; separate replay passed in
`4be8bf`. Peak RSS was 28,180,480 bytes for both. The jobs held
the shared global nonblocking proof lock, used a 59-second alarm
and checked RSS below one GiB. Large input was streamed, all
source hashes were streamed, and no queue was rebuilt. No temporary
files, bytecode cache, optimizer or Lean process was used.
Independent review is pending.

```sh
python3 -B scripts/fresh_compact21_sealed_e8_quad_triangle_activation_20261018.py --verify
```
