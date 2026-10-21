# A quartic order-exchange obstruction in dimension nine

20 October 2026. Own paper proof; no scientific worker, Lean or TeX run.
This concerns the literal integral forms below, not all lattices or
commensurability classes over this field.

Import the arithmetic and full actual first edge from
`compact_h9_biquadratic_unit_diagonal_first_edge_20261020.md`, SHA256
`05b6473eeadb4befa89c82af2f7002c7441ada3d1878c93d4630355a42685efe`,
and the full second edge from
`compact_h9_biquadratic_unit_second_edge_20261020.md`, SHA256
`b3a7c8f4e2afffe13fb49dd80ee2ed03095608003fabea85e2e6b74c9df2a52b`.
In particular the second source retains all three simultaneous primitive
mirrors before proving that only c is simple. No candidate enumeration
or independent execution is silently imported here.

Put

    Kfield=Q(s,phi), s^2=2, phi^2=phi+1, O=Z[s,phi],
    u=phi-1, alpha=1+s, k=s-1, d=2-s, a=u+s,
    N_n=<-a>+I_n.

The physical place has s,phi positive. N_n is positive definite at each
other real place, O is a PID, and all primitive reflector norms normalize
to1,2,d. Write C_n for the FULL integral-arrangement chamber at the
standard actual B_n vertex `(1;0^n)`.

**Theorem.** The full integral reflection group of N9 has infinite
covolume. In particular no compact Coxeter polytope has all its wall
reflections in this literal integral lattice group. The same proof
works for N_n with n>=8, with the inventory threshold explained below.
No lower-dimensional exclusion is asserted.

## 1. Two different actual rank-two faces

In the active three coordinates `(o,e1,e2)`, set

    beta=-e1+e2, r=o+phi*e1,
    c=C*o+U*e1+V*e2,
    C=s*alpha*phi, U=alpha*phi^2, V=phi*(phi+s).

The first two imported edge proofs establish that r,c are actual walls.
Their Gram matrix is

    G_P=[[d,-1],[-1,2]], det(G_P)=k^2,             (1)

so their I2(8) plane P is unimodular over O. The full actual face
perpendicular to P has the literal N_(n-2) arrangement. Its integral
unit time vector and the useful constants are

    T_P=H*o+I*e1+J*e2,
    K=a*alpha*phi=phi^3+s*phi^2,
    H=alpha*K, I=K*(phi+s), J=K*phi,
    T_P^2=-a, det[T_P,r,c]=k.                    (2)

All future time coordinates in these formulas are positive.

There is a second actual corner at `Z=(1;1,1,0^(n-2))`.
Indeed retain beta and the old B_(n-2) at the original B_n vertex.
The already ACTUAL c completes them to a positive-definite B2+B_(n-2)
diagram: beta^2=c^2=2, beta.c=-s. The actual-facet completion argument
can be checked directly here. Subtracting the retained projection from
c uses nonnegative coefficients, since the retained Gram matrix is
positive definite with nonpositive off-diagonal entries. Project the
old vertex toward this normal until its product with c is zero.
The resulting timelike point pairs nonpositively with every actual
wall by acuteness. Nonincident walls are strict; in particular the
released beta2 is strict because beta2.c<0. Thus this is the true
opposite vertex. Direct solution gives Z, whose norm is2-a<0.
This uses actualness of c; it is not a claim that any reflective
finite candidate must be a chamber wall.

The identity

    a*(a-2)=(alpha*phi)^(-2)                     (3)

gives a second unit time vector `T_Q=K*Z`, with `T_Q^2=-a`.
The full positive integral plane Q=T_Q^perp in these three coordinates
has the O-basis

    E=o+a*e2, F=e1-e2,
    G_Q=[[a*(a-1),-a],[-a,2]],
    det(G_Q)=a*(a-2), det[T_Q,E,F]=-1/(alpha*phi).
                                                               (4)

Consequently Q, not merely its root span, is unimodular and is an
orthogonal direct summand of the ambient lattice. Its oriented B2
simples are beta=-F and c=C*E+U*F; the latter identity uses U+V=a*C.
The root span O*beta+O*c has index ideal(s) in Q. It must NOT be
mistaken for an integral basis, or for the positive coordinate I2.

## 2. The B2 plane is fully sealed

Any ambient mirror crossing Q^perp has a spacelike projection there.
Its span with Q is positive definite. The generated integral reflection
subgroup is discrete and physically compact, hence finite. If the new
mirror is nonorthogonal to Q and not contained in it, its essential
rank-three group is irreducible. The finite rank-three classification
leaves only B3 containing a pair of mirrors at angle pi/4: A3 and H3
do not contain such a pair. A reducible rank-three group either acts
orthogonally to Q or has its nonorthogonal mirrors entirely inside Q.
This also handles three-generator subgroups of a larger finite group;
the generated essential rank is still at most three.

In B3 every pi/4 pair is the standard B2 pair. The B3 root orbit then
contains a root conjugate to beta or c, hence with primitive norm
EXACTLY2, having products (-1,0) or (0,-1) with beta,c. This is seen
in the orthonormal model: take the short mirror normal sqrt2*f1 and
long normal f2-f1 for the old pair; the long root f3-f2 has contacts
(0,-1). Interchanging which old mirror is short gives the other case.
All conjugating reflections preserve the ambient O-lattice, so the
new root is an integral vector, even when the old B2 was not a chosen
parabolic subsystem of a larger group.

The Q-projection of that root is respectively

    -beta-c/s  or  -c-beta/s.                    (5)

Neither is integral: c has the unit coordinate U, while beta has
unit spatial coordinates. But the unimodular direct summand (4)
forces the projection of EVERY ambient integral vector onto Q to be
integral. This contradiction excludes all nonorthogonal crossing
mirrors. It is stronger than testing native B3 extensions only.

Conversely every intrinsic reflector of O*T_Q+I_(n-2) extends by
identity on Q. Thus Q^perp, like P^perp, is a FULL literal N_(n-2)
face. The positive full reflection arrangement inside Q is exactly
the actual B2 germ at Z. There is no assumed B2 diagram-swap isometry;
in fact swapping beta,c alone is not integral.

## 3. Exchange the order of the two constructions

Use four of the original spatial coordinates; all other coordinates
are fixed. Construct first P then Q, or first Q then P. The full-face
identities above justify both operations in the actual full chamber.
The two resulting actual vertices have oriented germ

    I2(8) + I2(4) + B_(n-4).

The remaining factor is the full zero-height B_(n-4) of the final
literal N_(n-4). The sealing statements prevent any extra root from
mixing either positive plane with another component.

Here is the COMPLETE explicit active integral map. Work in the
original basis `(o,e1,e2,e3,e4)` with Gram diag(-a,1,1,1,1), and
regard each following ordered list as the columns of a 5-by-5 matrix:

    A0=K*T_P+K*e3+K*e4,
    A1=r,
    A2=c,
    A3=T_P+a*e4,
    A4=e3-e4;

    B0=H*T_Q+I*e3+J*e4,
    B1=T_Q+phi*e3,
    B2=C*T_Q+U*e3+V*e4,
    B3=o+a*e2,
    B4=e1-e2.                                  (6)

Define `D=B*A^(-1)` on this active space and identity on the remaining
spatial coordinates. Formula (6), together with (2), specifies every
native matrix entry without an unspecified choice of root or basis.
The A and B lists have the identical block Gram matrix

    diag(-a,G_P,G_Q).

They are O-bases: each is the composition of the two unit-determinant
orthogonal decompositions (2),(4); their determinants are the same
unit `-k/(alpha*phi)`. Thus D and D^(-1) are integral and D is an
isometry. It preserves the future sheet because it maps the displayed
future timelike A0 to future timelike B0.

The positive P simples A1,A2 map to B1,B2. The positive Q simples
`-A4, C*A3+U*A4` map to `-B4, C*B3+U*B4`. The remaining positive
B_(n-4) simples are fixed. These are FULL actual vertex germs, so
the integral full-arrangement isometry D maps an interior chamber
germ to another interior germ of the SAME C_n. Hence

    D(C_n)=C_n.                                (7)

This chamber identification, rather than a matching positive Gram
matrix alone, is essential to the argument.

## 4. Exact trace and infinite index

For transparency, compute the trace by the Gram-adjoint formula in
the three orthogonal blocks of (6). Their contributions are

    time: K^2*H^2-(K^2*H+K)*(I+J)/a,
    P:    -K*(2+C)/k,
    Q:    (I+J-2H)/(a-2).                       (8)

The identities

    I+J=a*(H-1), I-J=s*K, H=alpha*K,
    C=K-alpha, K^2=a/(a-2)

simplify their sum to

    tr_active(D)=alpha*K^3-alpha^2*K^2+(3+s)*K
                =196phi+121+s*(138phi+86)>5.    (9)

For example `K^2=14phi+9+s*(10phi+6)` and
`K^3=112phi+69+s*(79phi+49)`, so the last equality is an ordinary
four-coordinate ring identity. In N9 the full trace is

    196phi+126+s*(138phi+86)>10.

A finite-order real linear isometry has trace at most its dimension,
since every eigenvalue has absolute value1. Thus D has infinite order.

The full reflection group W acts freely on its open chambers. From
(7), the powers D^j therefore represent distinct cosets modulo W:
an element of W stabilizing its own open fundamental chamber is the
identity. Hence W has infinite index in O^+(N9).

The external arithmetic input is stated explicitly: the usual
Borel--Harish-Chandra arithmetic lattice theorem, together with its
anisotropic compactness criterion, makes O^+(N_n) a uniform lattice.
Here anisotropy follows immediately from any nonphysical positive
definite embedding: a nonzero isotropic K-vector would remain
isotropic there. Therefore infinite index implies infinite covolume
of W. Alternatively, already (7),(9) contradict compactness of its
chamber, since a compact hyperbolic polytope has a finite isometry
group. This completes the literal N9 obstruction.

For the stated uniform n>=8 version, both operations need the first
and second complete inventories in dimensions n and n-2. Those proofs
hold for dimensions at least6: their retained B_(m-2) has rank at
least4, so no exceptional finite extension competes with B_(m-1).
All arithmetic trace bounds, simultaneous-root formulas and full-face
arguments are dimension independent. The remaining zero germ works
in every rank. Thus (6)--(9) apply with n-4 fixed spatial coordinates,
and the full trace is `196phi+n+117+s*(138phi+86)>n+1`.
This does NOT extrapolate an unchecked small-B3 inventory to n=4.

Only the finite Coxeter classification, the elementary number-field
inputs imported above, and the explicitly named arithmetic lattice
theorem are external. The actual edge and order-exchange calculations
are the written arguments here; no numerical search, global novelty
claim, end-to-end formalization or execution certificate is asserted.
