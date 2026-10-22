# Norm239 original-base/fifth-germ period: full outside-author review

22 October 2026. **Full mathematical source review PASS; no correction requested.** Reviewed source: `results/compact_h9_prime239_original_base_to_fifth_e8_long_a1_period_candidate_20261022.md`, SHA-256 `dd15a2b4f51d17c0da3d26e828ebabadf06ca1dabe0f7c97fb87971b90dab848`,11,622 bytes/247 lines. This review covers the complete new period argument, at its explicitly named actual-germ and arithmetic inputs. It is separate from frozen publication acceptance `0c067d5d…` and adds no execution or Lean credit.

I did not author this period construction. I authored earlier norm79/finite57 supporting work and reviewed the norm239 released-S domain, implementation, actual execution and fifth-germ composition. Those inherited contributions are disclosed rather than presented as fresh independent reproduction of the edge census. For this review I read the entire new source, the complete original-base and odd-glue inputs, the fifth-edge composition and root acceptance, and the precise general arithmetic/chamber criterion in Section4 of the cited norm79 proof. I independently rederived the new algebra and geometric implications below. No worker, mathematical import, matrix execution, enumeration, test, Lean/TeX compilation, Git operation or network operation ran.

## Input and actual-germ boundary

All six tabled direct pins match current bytes: original-base `879a3feb…`, odd-glue `b107f0de…`, fifth-edge `30241d6a…`, its full other-author review `a3092c7a…`, root acceptance `1caee2d4…`, and actual runtime manifest `fc6b592c…`. The seventh, prose-pinned arithmetic/chamber interface `b9ddebcc…` also matches. The inherited finite-domain executions are not rerun or inferred merely from their hashes.

The actual input is the ORIGINAL module

    L=E8(Z[phi]) perpendicular-sum [[-2phi,1],[1,8]],
    phi²=phi+1, D=1+16phi, b=phi D=16+17phi.

Here Norm(b)=239. The binary determinant is−D, giving physical signature(1,1); its conjugate has positive first diagonal and positive determinant, so the full rank-ten form has physical signature(9,1) and definite conjugate. The named admissible arithmetic-lattice input is applicable to this literal module, independently of whether its reflection subgroup has finite covolume.

The original-base input establishes the FULL actual outward germ `-e0,...,-e7,r0` at future ray X, with `r0=X+2phi Y`, square2b. The fifth-edge input establishes the FULL actual outward germ `-e0,...,-e6,U,T` at its future negative endpoint P3, in the SAME full integral-arrangement chamber C. The E8 factor has the stated native ordering; T has square2b and is orthogonal to it. These are stronger premises than a matching diagram or the displayed known-wall inequalities. The earlier first-binary endpoint has isolated square2 and is correctly not substituted for the base square2b germ.

## Rational factors and all-ten-coordinate integrality

I checked the full calculation in source lines84–196 using phi²=phi+1. For `w=U+e7=aX+yY+hH`, the native identities `H²=2`, `H.e7=1` give

    w.e_i=0 (i<=6), w.e7=h,
    w²=2h, h=135+218phi, Norm(h)=131.

Thus `s_w(v)=v-(v.w/h)w` is the correctly normalized spacelike reflection, future preserving, fixing e0 through e6 and sending−e7 to U. Individual original-module integrality of s_w is neither established nor needed. From `T.U=0` I recover `T.w=hT`, `hT=1446+2343phi`, and

    hT/h=(-336+1077phi)/131,
    s_w(T)=[(4412+6963phi)X-(3470+5353phi)Y]/131=Tbar.

For `Tbar=rX+sY`, the vector `zeta=((r+8s)X+(2phi r-s)Y)/D` is orthogonal to Tbar and has square−2phi. The source matrix F sends X to−zeta and r0 to Tbar. The two source vectors span the binary space over K: their coordinate determinant is2phi, not an assumed O-unit. Therefore their exact norms and orthogonality prove F is a binary isometry. The determinants of these ordered pairs are2phi and−2phi, so det F=−1. Independently, writing the three entries of its numerator as A,B,C, I obtain `A²+BC=131²`. Also

    F(X).X=(-3470-5353phi)/131<0,

which proves the required future orientation for its timelike image. F fixes the original E8 subspace, with no claim that F alone preserves L.

The original covector of w is precisely `(-4663-7545phi,-6012-9725phi;0,...,0,h)`. I expanded both F-column products with this covector: the numerators are−57902−93665phi and−40872−65762phi, yielding

    F(X).w=-442-715phi, F(Y).w=-312-502phi,
    (442+715phi)/h=(156+169phi)/131,
    (312+502phi)/h=(700-246phi)/131.

Substitution into `G=s_w F` gives exactly

    GX=(3821+6184phi)X-(2945+4765phi)Y+(442+715phi)H,
    GY=(2691+4354phi)X-(2074+3355phi)Y+(312+502phi)H.

All six active numerators in the source table are131 times these displayed O-coefficients; I checked every entry, including the negative Y entries. On the remaining eight ORIGINAL basis vectors,

    G(e_i)=e_i (i<=6), G(e7)=e7-w=-U.

This establishes integrality on all ten original vectors. In particular H is expanded with e7 coefficient2; the e7 diagonal is `1-2h=-269-436phi`. There is no hidden replacement of the E7/H frame by an O-basis and no omitted old-e7 glue condition. Since s_w and F are K-isometries with determinant−1 each, det G=1. Its original integral matrix therefore has an O-integral adjugate inverse, proving equality G(L)=L rather than mere sublattice containment. Finally direct expansion gives `G(X)+2phi G(Y)=T`, exactly, so r0 maps to the outward T rather than−T or a rescaled squareclass.

## Future corner, full chamber and infinite order

Lines198–219 use a valid actual-sector argument. The future isometry G sends X to a future vector of square−2phi. Its product with X is also directly `-15313-24775phi<0`. It is orthogonal to every fifth-germ normal by the exact normal images. That rank-nine positive definite Gram has a unique negative orthogonal line; hence GX is a POSITIVE multiple of accepted future P3. This supplies the corner-ray sign, not just the correct hyperplanes.

Because G and its inverse preserve L and the form, conjugation permutes ALL original integral reflections and therefore the entire arrangement. The normal images are exactly

    G(-e_i)=-e_i (i<=6), G(-e7)=U, G(r0)=T.

They preserve every outward inequality at the two full actual germs. Local finiteness and the full-germ premise give open local interior sectors; G(C) and C share such a sector and hence are the same arrangement chamber. This does not require either rational factor separately to act on the arrangement, nor a finite-volume assumption on C.

The original-basis trace computation retains the seven identity diagonals and the corrected e7 diagonal:

    tr G=(3821+6184phi)-(2074+3355phi)+(1-2h)+7
         =1485+2393phi>10.

A finite-order real ten-dimensional matrix has complex eigenvalues of modulus1, so its trace has absolute value at most10. Thus G has infinite order. The argument is exact and does not depend on an approximate spectral radius or a germ-diagram automorphism.

## Arithmetic conclusion and exact scope

The last inference uses the explicitly inherited facts that Gamma=O⁺(L) is an arithmetic lattice of finite hyperbolic covolume and that the group W generated by ALL integral reflections has the full arrangement chambers as fundamental domains, with trivial W chamber stabilizer. These are applicable general facts, not norm79-specific lattice or period data. G belongs to Gamma by the new all-ten-basis calculation, normality of W follows by conjugation, and all powers of G stabilize C. Equal W cosets of two powers would put a nonidentity power in the trivial W chamber stabilizer. Infinite order therefore gives infinite index `[Gamma:W]`. The subgroup-index/covolume formula then gives infinite covolume for W and infinite volume of its full chamber C.

Consequently the source proves nonreflectivity of this exact original norm239 module at the stated arithmetic/full-germ boundary. It does not prove nonreflectivity of an entire field, family, form class or another invariant lattice; it yields no new compact H9 example or changed dimension bound. No fresh census, period checker, execution or formalization is claimed. All stages of the new argument pass without a source change.
