# Independent review of the quadratic13 two-edge route

20 October2026. **PASS: no mathematical flaw found in the full written
route.** This is an independent paper/source review, not an executed
census, a finite-volume proof or an external review.

Reviewed source:
`results/compact_h9_sqrt13_b9_two_edge_route_20261020.md`, SHA256
`b922647902d32beccd98fbe5a2781ddcbb9cafaeb10a12902eb3ebbce809547f`.
The source was read in full. All calculations below were rederived on
paper; no scientific script, import, numerical solve or proof job was run.

## Lattice and zero germ

For omega=(1+sqrt13)/2 and eta=omega+1, the norm-Euclidean argument is
valid: coefficient rounding gives a residual norm in [-13/16,1/2],
strictly inside (-1,1). The unit interval argument excludes every unit
between1 and eta, so positive units are eta^Z and totally positive units
are unit squares. The inertness of2 and the primitive coordinate/covector
Bezout argument force a primitive reflective norm to divide2. Both real
norms are positive, giving exactly the normalized classes1 and2.

The field and integral-product restrictions leave fundamental finite
bonds3 between long roots and4 between short/long roots. These statements
are about simple finite angles; they do not impose an incorrect restriction
on arbitrary mirror angles.

At p, the two-embedding bound on a spatial coordinate gives
|b|<=2sqrt2/sqrt13<1. Thus the full zero-root arrangement is the ordinary
B9 arrangement, not just a chosen nine-root subgroup. Admissibility gives
the required discreteness and local finiteness. Standard finite spherical
Coxeter classification and the all-wall acute property are the geometric
inputs used by the edge profile arguments.

## Both complete intervals

The first white B8 permits only the isolated short/long profiles and the
long-end B9 continuation. The upper v has norm2, and its collision bound
gives x^2<=R/2. Combining the conjugate bound with
sqrt(2eta)<sqrt13-1 forces x to be an integer. Only x=1,R=2 survives;
the full first simple wall is v. The endpoint q has norm -2eta and full
B8+A1 germ. The factor2 obstruction to an isometry between the p and q
rays is valid over Q(sqrt13).

For the second edge, the six profiles exhaust finite extensions of
B7+A1(v): isolated1/2; B2/A2 on v; B8; and the B9 bridge. In the source's
coordinates I obtain

    4(r,r)=Z^2-2eta X^2+2ell^2+4t^2.

This gives R=4,8,2,6,4,2 in the displayed inventory. Binary integrality
is exactly omega X=ell modulo2: X is even for ell0, while both its
coefficients are odd for ell1. The finite upper u yields
X^2<=R eta^2/2; conjugate positivity gives X'^2<=R eta/2.

I checked each strip, including the potentially tight R4 exclusion
omega'^2>eta/2, equivalent to11>3sqrt13. The full domains are precisely
Y=1 for R4; Y=1,omega,eta for R8; X=eta-2,eta,eta+2 for R6; and X=eta
for R2. The norm obstructions3,12,29,51 and the nonsquareness of2
eliminate exactly the claimed cases. The surviving isolated long root,
short root and bridge are all at the same contact, with no omitted sign
branch having X,Z>0.

## Actual corner and symmetry

The simultaneous-root selection is sound. With A=-e3 in the retained
B7 positive cone,

    u=d+A,                 w=2d+v+2A.

Thus d is the unique new simple root; u and w cannot be promoted to
separate walls. Its complete products give B9, and the omitted old roots
alpha1 and alpha2 are strict at H, with values -eta and -(3eta+1).
I independently obtain H^2=-eta. Local finiteness on the compact segment
to the timelike upper justifies the first-contact passage without any
assumption that the chamber is cofinite.

The displayed integral G maps alpha1 to v and alpha2 to d, fixes the
remaining simple roots, and maps p to H. Its columns preserve the Gram
form. Conversely Gv=alpha1 and Gd=alpha2; hence G^2=I. Its trace is6.
The full actual inward germs and preservation of the entire integral
mirror arrangement imply G(C)=C. This proves a finite chamber symmetry,
and supplies neither an infinite-order stabilizer nor a reflectivity or
nonreflectivity conclusion.

The literature screen is explicitly incomplete and is not used as a
proof input. I have not independently retrieved its three references in
this review. The valid outcome is a distinct literal admissible host,
two complete actual finite edges and an integral involution. Compactness,
finite volume and the full reflection index remain unresolved.
