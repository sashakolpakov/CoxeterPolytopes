# A mirror closing the last known initial H9 escape

19 October2026. Root written delta proof using the executed18-profile
search and its separate replay. The28-wall system's finite volume is
unresolved. A dedicated fixed-system and further edge-walking check is
being prepared, not claimed executed here.

The unchanged script
`scripts/search_h19_omega3_finite_profiles_20261019.py`, SHA256
`da6be4a933eba17b74701b1aac0699123410427cab4440045e4822c285b636e0`,
passed generation and separate replay under the shared-lock wrapper.
The exact certificate is `results/h19_omega3_finite_profiles_20261019.json`.
The two execution ledgers of the same prefix record0.404/0.392s,
under57MB process-tree RSS,187.53GiB disk and6.20/6.29GiB available RAM.
An earlier attempt deferred before child launch at5.7GiB RAM. No gates
were relaxed and no build intermediates were created.

There were19,314 profile/control square tests. Two field squares gave
four signed branches, of which two were integral norm-two roots.
Exactly one passed all27 old products and the finite endpoint test.
The search scope is0<x<=600, |x'|<=2; no all-height enumeration is asserted.

## The new short mirror

Write phi²=phi+1 and D=1+5phi. The surviving root h is

    h=(22+36phi,-18-29phi;
       86+136phi,172+272phi,258+408phi,215+340phi,
       162+255phi,108+170phi,54+85phi,129+204phi).

Its sparse covector is

    Gh=(-101-163phi,0;0,0,0,10+17phi,1,0,0,0).

Direct substitution gives h²=2: the binary contribution is
-8090-13090phi and the spatial contribution8092+13090phi.
All native coordinates and covectors are integral, so its reflection
preserves the lattice. Its27 old pairs are two label3, seven right
and eighteen strictly ultraparallel. Together with the separately
replayed27-wall system this supplies all378 pair checks for28 roots.

The old strict centre remains strict: (c,h)=-(418+695phi)/5<0.
Hence acute projection proves that every one of the28 inequalities
is an actual facet. The dedicated fixed checker will record these
projection witnesses separately.

## Four finite vertices outside the two H8 facets

The new endpoint on the old4 initial edge is

    p_h=p+tau_h*omega3,
    tau_h=(202+326phi)/(215+340phi),
    p_h²=-50D/(215+340phi)²<0.

The executed search checks its time cone, all27 old inequalities and
the new equality. Exactly the eight fixed white walls are active
among the old walls. Their Gram with h is positive definite, of
type2A4+A1. Thus p_h is an actual finite vertex.

The old finite vertices survive. We have(p,h)=-202-326phi<0.
Since p_s=p+Du*s and(s,h)<0 by the exact pair ledger, also(p_s,h)<0.
For p_q=p+((11-2phi)/10)*omega2, direct evaluation gives

    (p_q,h)=(1-52phi)/5<0.

These three are distinct from p_h. The new root is strictly ultraparallel
to both inherited complete compact H8 facet walls r,w, so it preserves
their two disjoint41-vertex sets. None of the four displayed vertices
belongs to either facet. Consequently the28-wall chamber has at least
82+4=86 finite vertices. This count inherits the earlier H8 facet
identification and vertex enumeration; it is not a new full census.

The former old4 spacelike endpoint has positive product with h and is
cut off. The other eight original p-edges already ended finitely in
the27-wall chamber. This closes the known initial escapes, but says
nothing yet about every edge at other vertices. No compactness,
finite-volume example, whole-lattice reflectivity or new dimension
record is asserted.
