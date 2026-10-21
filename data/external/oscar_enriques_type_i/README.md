# OSCAR Type-I Enriques lattice datum

`TauTaubarGenericEnriquesSurfaceNo172.mrdi` is the serialized datum used by
OSCAR 1.4.0 for `generic_enriques_surface(172)`.  OSCAR identifies this as the
very general finite-automorphism Enriques surface of Type I: its K3 cover has
Picard rank 19 and determinant 4, root invariant `A1 + E8`, and the Enriques
nef cone has twelve rational-curve walls.

Upstream source:

```text
https://raw.githubusercontent.com/oscar-system/Oscar.jl/v1.4.0/data/
TauTaubarGenericEnriquesSurfaces/TauTaubarGenericEnriquesSurfaceNo172.mrdi
```

The corresponding OSCAR 1.4.0 loader is

```text
https://raw.githubusercontent.com/oscar-system/Oscar.jl/v1.4.0/src/
AlgebraicGeometry/Surfaces/EnriquesAut.jl
```

and has SHA-256
`48cebfc2e924b914bfdc4d4d4f9d430a06439f8af35b798c8f7a556395989924`.
It fixes the serialization order by the assignment

```julia
SY, SX, L26, w, u = load(...)
```

and constructs `EnriquesBorcherdsCtx(SY, SX, L26, w)`.  Thus the fourth
serialized entry used in the project certificate really is OSCAR's Weyl
vector; the fifth is the unused preprocessing vector.  The same source
documents `SY` as the invariant lattice of the Enriques involution in the
K3-cover numerical lattice `SX`.

The upstream byte stream has SHA-256
`73e8e8612a9223156fb6024c62d51b32dc8a739ec9cd1cb09129a65acb73caa1`.
The archived copy has one terminal newline added by the repository patch
format; after removing that newline its bytes have exactly the upstream hash.

The upstream OSCAR package is copyright the OSCAR Development Team and
licensed under GPL-3.0-or-later.  See
<https://github.com/oscar-system/Oscar.jl/blob/v1.4.0/LICENSE.md>.

Mathematical provenance is Brandhorst--Shimada, *Automorphism Groups of
Certain Enriques Surfaces*, Foundations of Computational Mathematics 22
(2022), 1463--1512, DOI `10.1007/s10208-021-09530-y`.

The serialized object is a lattice/Borcherds context, not equations for a
surface.  Geometric existence is therefore supplied by the cited paper (in
particular Section 7.4) and not by the byte stream alone.  That section
identifies No. 172 with a very general Type-I Enriques surface, identifies
its K3 cover as Picard rank 19, and identifies its nef chamber with the
`12A` chamber having exactly twelve smooth-rational-curve walls.
