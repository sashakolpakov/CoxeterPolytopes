# Independent review of the complementary135th H20 wall

**PASS:** the new wall coherently extends the134-wall chamber to135 actual
facets and cuts the weak escape E0 used in the chain-contact obstructions.
It has no new finite vertices. A different strict spacelike escape survives,
so the chamber still has infinite volume.

Artifacts: [independent source](../scripts/fresh_h20_complementary_cap_20261016_review.py),
[certificate](fresh_h20_complementary_cap_20261016_review.json), and
[primary construction](fresh_h20_complementary_cap_20261016.md).
The reviewer imports no producer arithmetic. It uses the independent
sparse-radical field and interval engines from the earlier physical reviews,
reconstructs all107 old normalized covectors from literal coordinates, and
recomputes all26 chain normals with the independent recurrence. The old
physical chamber and parameter isolation remain reviewed dependencies.

## Exact defining equations and all contacts

The cap is `r=y-omega5`, with `omega5^2=30/11`. Its frame coordinates are
`y=(t,0,z)`. The right contact to wall95 is exact because the second frame
axis is its normalized root and omega5 is orthogonal to it. Let p and q
be the time and third-axis pairings with wall0. The exact pi/6 contact is

`p*t+q*z=-6/11-sqrt(3)`.

Write `t=alpha+beta*z`, with alpha and beta obtained from this equation.
The norm condition is `t^2-z^2=8/11`; the reviewer symbolically derives its
quadratic and uses the positive root. Rational interval isolation proves
the discriminant positive and identifies the selected branch. The norm2
and finite-angle identities follow from these defining equations, not from
an interval merely containing zero.

Independent strict checks give the full134-contact census:

| Contact | Count |
|---|---:|
| Right |18|
| pi/3, at support root5 |1|
| pi/6, at wall0 |1|
| Dotted |114|

All26 chain contacts are dotted. The preceding arc cap is also dotted to
the new wall, with pairing approximately-4.68959. Its physical frame is
reconstructed from the independently reviewed recipe before this mutual
pairing is evaluated. Every contact agrees with the primary certificate.

## Coherence, facets and the complete new finite-corner exclusion

The inherited strict center has strictly negative new-wall pairing. The
review reconstructs its negative norm and the new pairing independently.
Projecting this center to each norm2 wall gives a future timelike point
on that wall; acuteness keeps every other wall inequality strict. Thus all
135 listed walls are actual facets. The reviewed acute-extension preservation
argument retains the old finite vertices, giving at least170 in total.

There are exactly20 old finite-contact neighbors of the new wall: the18
support roots, wall0 and wall95. A new finite vertex requires19 of them.
The reviewer explicitly enumerates all20 subsets. Eighteen retain the
dotted pair0/87. Omitting0 instead retains the support plus the new wall,
whose Schur complement is `2-30/11=-8/11`. Omitting87 retains the affine
triple consisting of the cap,0 and1. Its exact Gram matrix has bonds6,3,2
and positive null vector `(sqrt(3),2,1)`. None of these subsets is positive
definite, proving there are no new finite vertices. The cap's dotted contact
with wall133 prevents that wall from supplying another finite neighbor.

The same computation checks strictly negative pairings with the three
inherited proper cusps. They survive. It does not classify all new ideal
faces or provide a complete ambient cone census.

## The old escape is cut; a new escape remains

The reviewer independently obtains `r.E0` approximately1.77808, strictly
positive. The previous epsilon-perturbed strict escape is cut as well.
This verifies a useful construction step even though no new finite vertex
is created.

For the replacement escape, use frame coordinates `(1,-1/100,-49/50)` and
support-dual coefficients-1 at root17 and-1/100000 at the other17 roots.
Its exact squared norm is

`1842000049/4000000000 > 0`.

The independent physical calculation checks all135 wall pairings strictly
negative and the pairing with the common future center negative. This gives
an open feasible patch on the ideal boundary and proves infinite volume.
The existence of this escape prevents a finite-volume or new-dimension claim.

Generation and a fresh saved replay pass, with recursive dependency hashes
checked. The scope is this coherent135-wall extension, its complete new
finite-corner exclusion, and its remaining escape. It neither exhausts
all possible further caps nor proves a general dimension20 obstruction.

```sh
python3 scripts/fresh_h20_complementary_cap_20261016_review.py --verify
```
