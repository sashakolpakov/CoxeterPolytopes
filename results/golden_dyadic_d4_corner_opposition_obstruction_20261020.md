# The full second-corner opposition is not integral

20 October 2026. Exact paper computation; no program or search.

Import the actual second corner from
`results/golden_dyadic_d4_extension_second_edge_20261020.md`, SHA256
`a894a30bb3b17e63a434049195931e8d6bd51a160b44d85141a820a08dd2e72e`.
Use its lattice `R_n=<-phi>+I_n+D4(O)` and roots a1,a2,r,a4.
The vertex is the ray of

    p2=(3;0^n;-1,-1,-1,-3), p2^2=3(4-3phi)<0,

with full positive germ `B_n+A2(a1,a2)+I2(5)(r,a4)`.

Let J fix p2 and negate its orthogonal complement. Compose J with
the longest element of the full finite germ. On B_n both factors
are minus identity. Since the longest elements of A2 and I2(5)
are negative diagram reversals, the resulting isometry G fixes B_n
pointwise and interchanges a1,a2 and r,a4 simultaneously. This
specifies G uniquely: p2 and the finite germ span the ambient space.

Put

    d=2phi-1=sqrt(5),
    v=a1-a2=(0;0^n;1,-2,1,0),
    Delta=r-a4=(2;0^n;1,1,1,-2phi-1).

Then `v.Delta=0`, `v^2=6`, and `Delta^2=4(phi+2)`. Both are
orthogonal to p2 and B_n. The two commuting reflections interchange
exactly the required pairs. Therefore the FULL simultaneous map is

    G=s_v s_Delta,
    Gx=x-(x.v)v/3-(x.Delta)Delta/[2(phi+2)].       (1)

In particular, for the actual integral time basis vector e0,

    Ge0=e0+phi*Delta/(phi+2)=e0+Delta/d,           (2)

because `phi+2=phi*d`. Its first D4 coordinate is `1/d`, which is
not in O: `1/d=(2phi-1)/5`. Hence `G(R_n) != R_n`. There is no
cancellation between the two component reversals on this column,
because e0 and Delta are orthogonal to v.

For an explicit check in a genuine integral D4 basis take

    beta1=e1-e2, beta2=e2-e3,
    beta3=e3-e4, beta4=e3+e4.

The full matrix is identity on I_n. Its other five columns, in the
basis `(e0,beta1,beta2,beta3,beta4)`, are equivalently the following
exact vectors:

    Delta=2e0+beta1+2beta2+(phi+2)beta3+(1-phi)beta4,
    v=beta1-beta2,

    Ge0   =e0+Delta/d,
    Gbeta1=beta2,
    Gbeta2=beta1,
    Gbeta3=beta3-v/3-phi*Delta/d,
    Gbeta4=beta4-v/3+Delta/d.                    (3)

These follow directly from (1); in particular they do not infer
nonintegrality from just one isolated factor. The prime above5 is
already visible in the full time column, and the prime3 remains in
the displayed D4 columns.

This rules out this particular integral chamber-period construction.
It proves no nonreflectivity, compactness, or commensurability
obstruction for R_n. A rational Gram isometry is insufficient to
preserve the full integral mirror arrangement.
