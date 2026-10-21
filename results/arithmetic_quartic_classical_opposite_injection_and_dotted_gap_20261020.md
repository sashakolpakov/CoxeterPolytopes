# Exact quartic field: large classical opposite injection and a sharp dotted gap

20 October 2026. Complete paper candidate awaiting independent review.
No scientific execution, search, lattice census, Lean, TeX or accepted-source
edits. No new ambient dimension bound or literal-host exclusion is claimed.

Let K=Q(sqrt2,sqrt5), s=sqrt2, phi=(1+sqrt5)/2, u=phi-1, and
O=Z[s,phi]. The arithmetic ground field is assumed to be EXACTLY K.
All assertions allow arbitrary finite spectator components.

**Results.**

1. At ANY vertex of a compact arithmetic Coxeter polytope over K, an
   outgoing direction in an A_r,B_r or D_r component, r>=5, shares its
   actual opposite wall with NO other outgoing direction.
2. For any actual dotted pair of unit walls, with product -gamma,

       4gamma^2>=3+s*phi, hence gamma>57/50>27/25.

   The scalar arithmetic gap is sharp; a high-dimensional compact
   realization of its equality parameter is not asserted.
3. At a farthest vertex, every A_r component, r>=10, satisfies

       ||x_A||^2>9/14.

The first result concerns actual finite replacements, not arbitrary
metric faces. The third concerns the ACTUAL component projection of
one farthest vector, not an edge-only covering relaxation.

## 1. Explicit inherited arithmetic and geometric inputs

Pin `results/fresh_sqrt2sqrt5_lanner_catalogue_20261016.md`, SHA256
`549a1152c68a8a70f55f404b9b6c3df64116e172092e196b460763613d85d9c1`,
and its independent review
`results/fresh_sqrt2sqrt5_lanner_catalogue_20261016_review.md`, SHA256
`9a266583591f50204b10050192d74e98eaba8f8294b34e020567a9d0dc35e403`.
Their inherited exact certificate is
`results/fresh_sqrt2sqrt5_lanner_catalogue_20261016.json`, SHA256
`37fd29a4ab6af46f82906ba1c0fd89c2d4fd96d95748dcdc4f66b4192743b278`.
No predecessor checker is imported or rerun here.

The precise inputs used are these:

* Vinberg's cyclic-product arithmeticity criterion: square-two cyclic
  products are algebraic integers in the ground field, and every
  nonphysical conjugate principal Gram is positive semidefinite,
  and positive definite when the principal block is nonsingular.
* A minimal indefinite principal diagram in a compact chamber is a
  Lanner diagram; affine principal diagrams are forbidden. A hyperbolic
  principal subset retains the exact arithmetic ground field.
* For exact K, the ONLY order-three arithmetic Lanner diagrams are
  (2,5,8) and(4,5,5); none has order at least four. Order two is the
  parameter family tau in O, tau>4 and all other conjugates in(0,4).
* Compact simplicity/acute-face correspondence and finite spherical
  Coxeter classification. A finite irreducible containing label5 has
  rank at most4; one containing label8 has rank2.

The catalogue derives the triangle list by complete exact finite-angle
tests; absence of higher orders uses the published arithmetic Lanner
classification, whose ground fields are only Q,Q(sqrt2),Q(sqrt5).
The pinned review records its primary verification of Nikulin,
Sections2 and3.1, https://arxiv.org/pdf/0708.3991 . That classification
is an explicit inherited input here, not a new own classification.

For clarity, the same-field assertion also follows from conjugate
positivity. A nonsingular Lanner determinant is negative physically
and positive at every nonphysical embedding. If it belonged to a
proper subfield, a nonidentity embedding fixing it would preserve
its negative value. Thus its cyclic field is the whole ground field.
This argument does not assert that every normalized Gram entry is in K.

## 2. No opposite sharing for a large classical component

Let S be the full positive definite vertex diagram and c an actual
outside wall. Suppose c closes two distinct directions i,j, one of
which lies in A_r,B_r or D_r, r>=5. Then both

    (S minus i)+c and (S minus j)+c

are the positive definite diagrams of the two ACTUAL endpoints.
The full S+c spans the ambient Lorentz space and is indefinite.
Choose a minimal indefinite principal set T in it. T must contain
c,i,j: omitting c leaves the old finite vertex, while omitting either
i or j leaves a subset of an actual finite endpoint. The catalogue
therefore forces T to be exactly the triangle c,i,j.

Consider its old bond i-j.

**Old bond2.** The triangle is(2,5,8), so c has a5 contact to one
old node and an8 contact to the other. If the8 node belongs to a
nontrivial old component, it has an old neighbor retained upon
deleting the other node (the two marked old nodes are orthogonal).
This contradicts the finite endpoint: an8 bond cannot have any
nonright extension in a spherical diagram. Thus the8 node must be
an old A1. The5 node is then in the large classical component.
Deleting the A1 retains the ENTIRE large component together with its
5-attached c, a connected finite diagram of rank at least6 containing
5. This contradicts the rank-four finite bound.

**Old bond3.** Neither allowed arithmetic triangle contains3.

**Old bond4.** This is possible in a large classical component only
at the terminal4 of B_r. The triangle would be(4,5,5). Deleting
the old short terminal node retains A_(r-1) and the5-attached c.
Their connected positive diagram would have rank r>=5 and contain5,
again impossible. There are no other old native bonds in A/B/D.

These cases exhaust T, proving opposite injection. In particular,
no triangular two-face at this vertex can contain an outgoing edge
belonging to such a large classical component: its third wall would
be the common opposite of its two outgoing edges.

No assertion that an arbitrary spectator-orthogonal metric face is
Coxeter, or that an ambient lattice splits, was used.

## 3. Sharp scalar arithmetic gap for every dotted pair

Let tau=4gamma^2. By Section1,

    tau in O, tau>4, 0<sigma(tau)<4 for all three sigma!=id.

Put T0=3+s*phi. We prove tau>=T0, retaining equality. Suppose
4<tau<=T0 and write uniquely

    tau=alpha+s*v, alpha=A+B*phi, v=C+D*phi,
    A,B,C,D in Z.

Use a prime for the phi conjugation. Since alpha-s*v is in(0,4),
the physical tau>4 forces v>0 and alpha>2. The other two conjugates
give |v'|<s. Consequently

    0<v<(T0-2)/s=phi+1/s,
    Tr_(Q(phi)/Q)(v^2)<(phi+1/s)^2+2
                         =phi*(1+s)+7/2<8.

For the last inequality phi<13/8 and s<10/7 give
phi*(1+s)<221/56<4. The trace is the integer

    2C^2+2CD+3D^2.

Its complete positive-element list at values at most7 is

    v in {1,u,phi,u^2,phi^2}.

Indeed |D|>=2 gives trace>=10; D0 and D=+-1 give exactly this
list and its negatives. The elements u,u^2 violate |v'|<s.
The element phi^2 violates v<phi+1/s. Only v=1 or phi remains.

If v=1, then

    4-s<alpha<=3+s*u,   s<alpha'<4-s.

Hence 0<alpha-alpha'=B*sqrt5<3<2sqrt5, so B=1.
But alpha'=A-u>s, and s+u>2, give A>=3. This implies
alpha>=3+phi>3+s*u, a contradiction.

If v=phi, then

    s*phi<alpha<=3,   s*u<alpha'<4-s*u.

Here -sqrt5<alpha-alpha'<sqrt5: the lower bound follows from
s*phi>2 and alpha'<4, while the upper uses
3-s*u<sqrt5, equivalently 2u<s. Therefore B=0. Since
2<s*phi<alpha<=3, the integer alpha must be3, and tau=T0.

Finally T0 itself has conjugates

    3+s*phi, 3-s*phi, 3-s*u, 3+s*u.

The first is greater than4 and the other three lie in(0,4).
It generates K: (T0-3)^2=2phi^2 recovers phi, then s.
Thus T0 is the SHARP scalar parameter, not a rounded bound.
Since s>7/5 and phi>8/5,

    gamma^2>=T0/4>131/100>(57/50)^2,

which proves the displayed rational gap.

This does NOT by itself rule out every generic far-helper weak
packet in every rank: its large-r asymptotic loss can exceed .61.
Any such use still needs the exact finite index range or additional
profile constraints.

## 4. A short actual A_r projection consequence, r>=10

Now choose a farthest vertex, put N=r+1>=11, and use unit A_r
normals n_i=(e_i-e_(i+1))/sqrt2 and normalized duals u_i. Let
h_i=x_A.u_i. Every relevant actual opposite is distinct by Section2.
The actual reciprocal-shadow/acute-wall inequality is

    F(h_i)F(h_j)>=u_i.u_j, F(t)=2t/(1+t^2).

This is the geometric interface pinned in
`compact_large_an_actual_half_uniform_transfer_20261020.md`, SHA256
`7351b471ba3b0edaf0a524ccf082010e458b7fedd8d063491a75632223bd70e7`.
Because h_i in(0,1) and log F is concave, the outer pairs1/2 and
(r-1)/r have sums>26/25, while the inner3/4 and(r-3)/(r-2)
have sums>5/4. The least correlations occur at N11; exact seals are

    F(13/25)^2=105625/157609<1341/2000<sqrt(9/20),
    F(5/8)^2=6400/7921<81/100<sqrt(21/32).

The first cross-product margin is103669 and its final square margin
is1719; the second margins are1601 and48, respectively.

Set v_i=u_i+u_(N-i), i=1,...,4, and W=v1+v2+v3+v4. In the
standard sum-zero coordinates,

    v_i=sqrt(N/[i(N-i)])*(1_(first i)-1_(last i)),
    v_i.v_j=2N*sqrt(min(i,j)/max(i,j))
                     /sqrt((N-i)(N-j)).

Every displayed Gram entry decreases with N. At N11, the diagonal
sum is13277/1260. The six off-diagonal contributions and exact
positive rational upper bounds are

    44/sqrt180<329/100,  44/sqrt240<57/20,
    22/sqrt70<263/100,   44/sqrt108<106/25,
    44/sqrt126<98/25,    11sqrt(3/14)<51/10.

They sum to less than2203/100. Squaring positive numbers verifies
each, so W^2<102587/3150<163/5 for every N>=11.
The four actual pair laws give

    x_A.W>2*(26/25)+2*(5/4)=229/50.

Cauchy proves ||x_A||^2>9/14, since

    (229/50)^2-(9/14)*(163/5)=337/17500>0.

No spectator direction or energy was counted in W, so the theorem
is an arbitrary-spectator ACTUAL projection statement.

## 5. Module boundary and contribution

These results apply to every compatible arithmetic lattice over
exact K, without choosing one. They do not identify that lattice
with a closed literal diagonal host. For example the admissible
integral binary Gram

    [[2,-T0],[-T0,2T0]]

has integral reflection Cartan entries -T0 and -1 and determinant
T0*(4-T0). Its determinant has rational norm -31: Norm(T0)=31
and Norm(T0-4)=-1. Thus even the sharp dotted parameter need not
give a unimodular binary block. Adding an orthogonal positive block
does not prove compactness or make this a new construction candidate;
it simply shows why a literal unimodular-host cutoff cannot be
silently substituted for an all-lattice argument.

Finite57_research derived opposite injection, the exact four-place
gap and the fixed A-window corollary. Root identified the useful
field-forcing direction and explicitly cautioned against promoting
the scalar gap to an all-r weak-packet exclusion. The inherited
catalogue, geometric tools and their authors retain their stated
contributions. No new ambient ceiling is asserted before separate
composition and source-bound review.
