# Hostile scope audit: rank-four K3 infinite-field no-go

**Verdict (2026-08-10): PASS at the stated weak boundary; the maximal
same-chamber field gate remains open.**

This is a self-contained hostile scope replay, not a second provenance claim.
It attacks both the exact arithmetic and the temptation to promote the family
to reflection-maximal Coxeter return.

## Frozen owner inputs

| artifact | SHA-256 |
|---|---|
| results/rank4_k3_correspondence_infinite_field_no_go.md | 3406e7d4f43ac1c96c13bd93dd8f6eae7e90564efcf06667262de4c155dd26c8 |
| scripts/certify_rank4_k3_correspondence_infinite_field_no_go.py | cb70ef2f39da8309cac648a87fd0731187812208a955dd3ba587a1a8c83e82f9 |
| data/enumerations/rank4_k3_correspondence_infinite_field_no_go.json | d9792dd05caa49335b593ae621296f1c48721328e02b56a6823d715a04c5cd88 |
| tests/test_rank4_k3_correspondence_infinite_field_no_go.py | 3d67e694fd4a8192b9da7c90b2499e115141d8b596f0ce30a143d1278f8268ff |

## Attack table

| gate | hostile check | verdict |
|---|---|---|
| field | \(a^2-1\) is nonsquare for \(a\ge2\), so \(x^2-2ax+1\) is irreducible and defines a real quadratic field | PASS |
| full trace | each companion block is the regular \(A_a=\mathbf Z[t_a]\) module; two blocks give \(A_a^2\), and self-adjointness recovers the \(k_a\)-bilinear trace form | PASS |
| self-adjointness | independent multiplication gives \(C_a^tB_0=B_0C_a\) and \(C_a^tB_a=B_aC_a\); the quadratic polynomial annihilates both blocks | PASS |
| signature | \(B_0\) has signature \((1,1)\); \(B_a\) has negative leading minor and determinant \(2(a-1)>0\), hence is negative definite | PASS |
| real-place split | on the \(C_a\)-eigenline \((1,-\lambda)\), the \(B_0\)-square is \(\lambda^2-1\); the large place contributes one positive direction and the conjugate place none | PASS |
| primitive K3 embedding | pullback of \(U^3\oplus\langle-2\rangle\) is exactly \(2B_0\oplus2B_a\); after \(y_1\leftarrow y_1-y_2\), rows \(e_1,e_2,e_3,v\) give a determinant-one minor | PASS |
| exact Picard carrier | the primitive complement has signature \((2,16)\); a generic period avoids all extra rational classes, while the square-two class makes the K3 projective and minimal | PASS |
| correspondence | the divisor-product cycle with intersection-dual first factor acts as \(\tau_a\) on \(\operatorname {NS}_{\mathbf Q}\); the operator is integral although the cycle is only rational | PASS |
| restricted nef equality | the conjugate eigenspace contains no nonzero rational vector, so no root mirror contains the selected eigenspace; chamber intersection, Weyl framing, and closure give equality with one complete restricted-root component | PASS |
| infinitely many fields | a fixed squarefree field contributes only \(O_d(\log X)\) Pell first coordinates up to \(X\); finitely many fields cannot cover all \(a=2,\ldots,X\) | PASS |

## Promotion traps rejected

The complete restricted arrangement need not have a Coxeter component, a
finite number of sides, or finite hyperbolic volume.  Nothing in the matrices
proves that its facet reflections saturate the ambient reflection group.
Likewise, Weyl framing returns the component selected after constructing the
K3 surface; it does not return an independently prescribed chamber.

Therefore the valid implication is exactly:

    rank 4 + minimal K3 + primitive full trace NS package
    + integral self-adjoint NS action from a rational correspondence
    + exact equality with some complete restricted-root nef component
    does not force finitely many fields.

The following stronger statement is **not certified**:

    reflection-saturated, reflection-maximal, finite-volume, prescribed same-P
    eigen-nef packages occur for infinitely many fields.

The standalone hostile tests reconstruct the integer matrices and embedding
without importing the producer, pin every owner byte, replay the archive, and
fail if any of the excluded maximality claims is promoted.
