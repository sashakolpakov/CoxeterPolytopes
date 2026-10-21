# Complete ambient cone after closing the r8 quadrilateral

The exact60-wall H9 chamber has **188 finite vertices, zero ideal extreme rays and2,680 exterior extreme rays**. It retains all184 old finite vertices. Both the H4+I2(8)-supported compact tetrahedron and the new r8 compact quadrilateral are whole actual faces. An explicit strict spacelike escape still proves infinite volume.

The [clipping certificate](../scripts/fresh_nonarithmetic_pi8_quad_cone_20261017.py) and [saved output](fresh_nonarithmetic_pi8_quad_cone_20261017.json) start from the [complete58-wall cone](fresh_nonarithmetic_pi8_tetrahedron_cone_20261017.md) and use the [exact two-wall construction](fresh_nonarithmetic_pi8_quad_completion_20261017.md). Its physical normals are assembled in literal order: the inherited55, the three tetrahedral additions55–57, then the quadrilateral additions58–59.

| Added wall | Old rays removed | New crossing rays | Finite vertices after clipping | Exterior rays after clipping |
|---|---:|---:|---:|---:|
|58|40|80|186|2,652|
|59|30|60|188|2,680|

Every old ray has a strict sign at each new clipping wall. For each opposite-sign pair with at least eight common active normals, a certified rank-eight minor proves that they span an edge. Their positive linear combination supplies the crossing ray. All genuine crossing edges occur in this enumeration. Reconstructing each new ray from nine independent active equations and one nonzero anchored coordinate tightens its interval representation. Uniqueness proves exact projective equality; interval overlap is only a consistency check.

The final certificate checks172,080 wall inequalities. Nonincident products and every norm sign are strict. Incident equalities follow from the inherited equations and exact clipping identities. The four new finite supports are

    S6 + [7,49,58]
    S6 + [7,55,58]
    S6 + [30,55,59]
    S6 + [55,58,59],

where S6={0,1,2,3,5,6}. Consequently the original H4+A2-supported three-face has18 finite and8 exterior extreme rays and remains open. Its r8 facet, supported by S6+[55], has precisely four finite rays. The distinct H4+I2(8) three-face still has precisely its four tetrahedral vertices.

For the infinite-volume witness, take the saved exterior ray with support

    [0,1,2,3,5,6,7,33,49]

and add2^-34 times the saved strict future center. Its square is approximately7.19446734207e-11, strictly positive by rational bounds. Its center product and every one of its60 wall products are strictly negative. This exposes an open part of the ideal boundary; the absence of ideal extreme rays does not imply compactness.

Generation and a fresh saved replay pass. A separate agent source review checked the complete clipping criterion, the direct ray reconstruction argument, normal order and face supports. Its separate local S6 calculation agrees with18 finite and8 exterior rays. These checks are internal; this package does not independently reconstruct the algebraic root recipes.

```sh
python3 scripts/fresh_nonarithmetic_pi8_quad_cone_20261017.py --verify
```
