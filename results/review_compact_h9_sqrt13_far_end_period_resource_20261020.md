# Independent review of the quadratic13 far edge and chamber period

20 October2026. **PASS; no mathematical flaw found.** Full source read
and independent paper derivation. No checker, matrix program, import,
Lean or other scientific job was run in this review.

Reviewed source:
`results/compact_h9_sqrt13_a9_far_end_and_period_20261020.md`, SHA256
`86ffab56f4696298fac4ae34b8bc8772dd1e7a6daf2e95e48d56291e810bbcf2`.
Its inputs are the two-edge proof b9226479 and A8-edge proof e3e62a6e;
my complete earlier reviews20b3381d and ad0325a2 cover those sources.

## Exact plane, reconstruction and whole interval

I independently obtain the stated P0,R0 vectors, norms -eta,9 and
zero mutual product. Their determinant in native(time,first-spatial)
coordinates is eta, a unit. Hence they are an integral basis of the
entire retained-plane lattice, not merely a K-basis. In particular both
isolated profiles reduce to the genuine native equation
9B^2-eta A^2=q. The different split places over3 exclude q1 and q2,
respectively, at all heights.

For an attached profile I rederive

    A=(2eta-2)x-(eta+4)z+5-eta-k,
    Z=9[-eta x+(eta+2)z-1]+(eta+3)k.

The first eight spatial entries and the ninth entry in the source have
exactly the required white products. Their white projection has squared
norm k(9-k)/9. Thus Z^2-9eta A^2=R is correct for all eight placements.
The (x,z) coefficient matrix has unit determinant eta-3=1/eta. Equivalently,
the orthogonal vector with the same(time,first-spatial) coordinates has
P0 coefficient A+k+eta-5 and R0 coefficient
(Z-(eta+3)k)/9+1. This proves both necessity and sufficiency of the stated
glue; no abstract A8 symmetry is silently extended through a nonunit index.

The finite upper b is primitive, norm2 and an integral reflector. Its
postselection native/cap products all check. More importantly the interval
does not assume A>0: the actual start/end inequalities imply Z>0,
Z>=(eta+3)|A|, with strictness on the starting side. Therefore10A^2<=R,
and conjugate positivity gives A'^2<=R eta/9. These force A to be an
integer, leaving only0 at R1/R4 and -1,0,1 at R10. The modulus9 exclusions
and sqrt10 obstruction are valid. At R10 the A=-1 vector is exactly the
negative highest old A9 root -b0 through Q; the strict starting inequality
correctly excludes it. It is not a second outgoing wall. The remaining
A=1,k=1 vector is b.

The complete primitive alphabet, finite rank9 extension classification
and local finiteness therefore give a genuine first simple wall and the
full new A9 inward germ. No old-cap filter or finite-height assumption
enters that passage. I also obtain

    Q'=(53eta+15;(33eta+10)^8,21eta+8),  Q'^2=-10eta.

The removed alpha8 is strict there with value -12eta-2; the other old
nonwhite walls are strict too. Thus there is no hidden feasibility failure
in interpreting the endpoint geometrically.

## Integral symmetry, trace and global conclusion

The longest white Weyl element w0 fixes P0,R0. The point inversion is

    J_P=-I-(2/eta)P0*(P0,.).

Since eta is a unit and P0 is integral, this is an isometry of the full
native lattice, with integral inverse. Its future orientation follows
from fixing the future timelike P0. It commutes with w0. I independently
check w0 alpha8=b0, H alpha8=b, H beta_j=beta_(9-j), and H Q=Q'.
The full actual inward germs, rather than a partial root configuration,
therefore imply H(C)=C. The inherited G has the same global property.

For an independent trace calculation, the four gamma roots are mutually
orthogonal of norm2, so

    H=-I-(2/eta)P0*(P0,.)+sum gamma_j*(gamma_j,.).

I obtain (P0,G P0)=-552eta-167, and
sum(gamma_j,G gamma_j)=-66eta-16, with the individual rows exactly as
printed. Consequently

    tr(HG)=-6-(2/eta)(-552eta-167)-66eta-16
          =268eta+80>10.

This verifies the period independently of an unexecuted numerical matrix
calculation. H and G are individually involutions, but the trace proves
their product has infinite order.

The chamber-germ continuation and free action of the full reflection group
on its chambers are standard explicit geometric inputs. They give infinitely
many distinct cosets of W in O^+(L). The standard arithmetic-lattice theorem
for the admissible integral form, followed by the lattice index/covolume
fact, gives an infinite-volume reflection chamber. The exclusion is of the
literal L=<-eta>+I9 and any compact reflection chamber whose mirrors preserve
it. No exclusion of K-similarity classes, commensurability classes or all
compact9 polytopes follows. Exact arithmetic execution remains a separate
pending validation step, not a replacement for these geometric inputs.
