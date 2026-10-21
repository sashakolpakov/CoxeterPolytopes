# Odd silver diagonal chambers: two exact edges and lower-face identities

Status: paper derivation, 20 October 2026. No worker, root enumeration,
Lean job or geometric compactness census has been run for this note.

Put `s=√2`, `α=1+s`, `d=2+s=sα`, and `O=Z[s]`. Consider

    L_n = <-α> ⊥ I_n,        M_n = <-s> ⊥ I_n.

The real embedding `s>0` is understood. Both forms are positive definite at
the other embedding, so their integral isometry groups are discrete
cocompact arithmetic groups. This does not make their reflection subgroups
cofinite. Write `x` for the time coordinate, and orient positive-time roots
so that their product with `p=(1;0,…,0)` is negative.

**Proved here:** the full `L_9` chamber at its B9 vertex has a uniquely
determined first coordinate edge and a uniquely determined first A8 edge.
The latter ends at the full A9 corner displayed in Section 4. There are also
uniform full-arrangement face identities. **Not proved:** compactness,
reflectivity, or nonreflectivity of `L_9`, or an integral identification with
the previously excluded silver E8 host. The latter issue must be settled
before any scientific execution or claim of a new construction host.

## 1. Complete norm alphabet and the actual initial chamber

The ring O is norm Euclidean, hence a PID. Its units are `±α^j`, and its
totally positive units are squares of units. If a primitive vector of the
unimodular lattice `L_n` defines an integral reflection, its squared norm
divides 2: apply the reflection integrality condition to all coordinates
and use the primitivity Bezout identity. As `(2)=(s)^2`, the complete
totally positive norm alphabet, modulo multiplying the vector by a unit,
is

    q ∈ {1, 2, d}.

Conversely every integral vector with one of these norms is reflective,
because each norm divides 2. A norm-2 vector can still be nonprimitive;
retaining such vectors only overcounts the same mirrors.

At time zero, write the spatial vector as `a+s b`, `a,b∈Z^n`. For norm 1
the constant equation `||a||²+2||b||²=1` gives `±e_i`. For norm 2 it gives
`±e_i±e_j`, together with `±s e_i`, the latter parallel to the norm-1
mirrors. Norm d is impossible: the s coefficient of a sum of integral
squares is even. Thus the *entire* zero-height root arrangement is Bn.
Its actual inward germ has outward roots

    a_i=(0;−e_i+e_(i+1)), 1≤i<n,       a_n=(0;−e_n),

and inequalities `x_1≥…≥x_n≥0`.

Local finiteness supplies the usual actual-edge interpretation below. A
compact hyperbolic set bounds the physical reflection matrices. The other
embedding is orthogonal for a positive definite form and bounds their
conjugate matrices. Bounded elements in both embeddings of O form a finite
set. Hence only finitely many integral mirrors meet a compact hyperbolic
set. In particular a finite proposed upper endpoint guarantees a first
wall on the relevant ray. Its old incident actual walls, together with
the new simple wall, form a finite spherical Coxeter diagram.

## 2. The coordinate edge: a unique long wall

Assume `n≥5`, and release `a_1`, retaining `B_(n−1)` on coordinates 2,…,n.
The explicit integral primitive root

    r=(1;α,0,…,0),          r²=d,

has `r.a_1=−α`, is orthogonal to the retained walls, and meets the edge at

    p_r=(1;1,0,…,0),       p_r²=−s<0.

The pair `a_1,r` has angle π/8.

The complete finite continuations of the retained B diagram are an
isolated root or a Bn extension at its long end. In the isolated case a
root has spatial coordinates `(y,0,…,0)` and

    y²−αx²=q,       q=1,2,d.

In the connected case its norm is 2, its spatial coordinates are
`(y,1,0,…,0)`, and `y²−αx²=1`. This connected case is essential; it must
not be dropped by assuming every new wall is orthogonal to the old B
component. The finite spherical classification gives no other case for
the retained ranks in question.

If such a wall occurs no later than r, then

    0<x,       x²≤q/d,       x'²≤q'α,

where for the connected case q means its binary residual 1. The first
inequality follows by comparing the slope `αx/y` with the upper slope 1;
the second is the conjugate positive-definite inequality. For
`x=a+b s∈O` these strips give:

| residual q | safe bounds | possible x |
|---|---|---|
| 1 | `x<3/5`, `|x'|<8/5` | none |
| 2 | `x<4/5`, `|x'|<9/4` | none |
| d | `x≤1`, `|x'|<6/5` | 1 |

Indeed `b=(x−x')/(2s)`. In the first and third rows b=0. In the second
row b can only be 0 or 1; b=1 would require `x=s−1`, whose conjugate has
absolute value α>9/4. The final row gives `y=α`, hence precisely r.
This proves r is the unique first simple wall of the *full arrangement*,
and `p_r` is an actual `B_(n−1)+A1` vertex.

## 3. Full lower-face identities

### 3.1 The long facet is the full M chamber

The equation of r-perpendicularity is `x_1=x`. In the integral basis
`(1;1,0,…,0),e_2,…,e_n` its lattice is exactly `M_(n−1)`.

First its own complete primitive norm alphabet is again `{1,2,d}`. The
reflection integrality condition allows a norm ideal dividing `(2s)`.
The only extra valuation would be 3 at `(s)`, represented by `2d`. It
forces every spatial coordinate to be divisible by s. Primitivity then
forces the time coordinate to be an s-unit. Its contribution `−s x²`
has valuation 1 and the spatial contribution has valuation at least 2,
contradicting norm valuation 3. This removes the extra class. Therefore
every intrinsic M root extends as an integral ambient L root.

Conversely each ambient mirror meeting the interior of this hyperbolic
facet induces an intrinsic M root mirror. Here one must permit all
angles between arbitrary arrangement mirrors, not just native simple
angles. The finite values of `4cos²θ` in `Q(s)` are

    0,1,2,3,2+s,2−s.

Integrality and the three norm classes leave the following nonorthogonal
cases with r (signed versions are identical). If `v²=2` and
`v.r=−α`, then `s v+r` is perpendicular to r, integral, and has norm
`2−s`, a unit-square associate of d. If `v²=2` and `v.r=−1`, then
`s v+α^(-1)r` is perpendicular and has norm d. If `v²=d` and
`v.r=−α`, then `s v+r` is perpendicular and has norm d. Norm-1
mirrors have no nonorthogonal intersecting case. Orthogonal mirrors are
already intrinsic roots. Thus no induced mirrors are omitted, and the
actual r facet is a chamber of the *full* M arrangement.

### 3.2 Coordinate short facets

A norm-1 coordinate facet in `L_n`, respectively `M_n`, has intrinsic
lattice `L_(n−1)`, respectively `M_(n−1)`, and the complete induced
arrangement equals its full integral arrangement. The only possible
nonorthogonal crossing is a norm-2 ambient root with product ±1; its
orthogonal projection is an integral norm-1 root. Norm-1 and norm-d
ambient roots cross a coordinate short facet only orthogonally. The
converse extension follows from the same norm alphabet.

### 3.3 A unimodular I2(8) face, uniformly in n

The Gram matrix of `a_1,r` is

    [[2,−α],[−α,d]],       determinant 1.

Its positive rank-two lattice is therefore a direct orthogonal integral
summand of L. Its orthogonal complement has basis
`u=(1;1,1,0,…,0),e_3,…,e_n`, where `u²=1−s=−α^(-1)`.
Replacing u by the unit multiple αu identifies the complement with
`L_(n−2)`.

At any interior point of this codimension-two face a crossing mirror
would give a finite rank-three root system containing the two actual
simple I2(8) walls. No irreducible finite rank-three Coxeter system has
an edge labelled 8. Hence such a mirror is orthogonal to the I2(8)
span. Conversely every root of the complementary lattice extends to L,
by the alphabet in Section 1. Thus this *actual* face has exactly the
full lower L arrangement. For example compactness of L9 would imply
compactness of its actual L7, L5, L3, L1 faces. The implications alone
give no obstruction: low-rank chambers can be compact.

## 4. The L9 A8 edge: complete tiny census

Now release the short wall `a_9`, retaining `a_1,…,a_8`, an A8 chain.
For an added connected root write its spatial coordinates as

    ((y+δ)^k,y^(9−k)),      Z=9y+kδ.

The full finite rank-nine inventory is:

| type | root norm | δ, k | equation `Z²−9αx²=R` |
|---|---|---|---|
| isolated | 1,2,d | δ=0 | R=9,18,9d |
| A9 | 2 | δ=1, k=1 or 8 | R=10 |
| D9 | 2 | δ=1, k=2 or 7 | R=4 |
| B9, short normalization | 1 | δ=1, k=1 or 8 | R=1 |
| B9, equal-norm normalization | 2 | δ=s, k=1 or 8 | R=2 |

Both B9 rows are necessary over O. In the last row the new norm-2 root
has product `−s` with a norm-2 white endpoint, producing angle π/4.
The field/ring conditions exclude all other scalings of the finite
spherical A9/B9/D9 diagrams. The isolated cases retain all three
primitive norm classes. For a forward wall `x>0,Z>0`.

An explicit finite upper wall is

    a=(2α;α+1,α,α,α,α,α,α,α,α),       a²=2,
    Z_a=9α+1=10+9s,                    x_a=2α.

It has `a.a_1=−1`, is perpendicular to `a_2,…,a_8`, and
`a.a_9=−α` is ultraparallel. Its product with the already actual r is
`−d`, also ultraparallel. It meets the A8 edge at a timelike point.

Comparison with this *independently supplied* upper wall and conjugate
positivity give

    x²≤(2R/5)α²,       x'²≤R'α/9.

No cap-product rejection is used in the census. For R=1,2,4,10, the
complete coefficient strips are:

| R | coarse enclosing bounds | complete x list |
|---|---|---|
| 1 | `x<8/5`, `|x'|<3/5` | none |
| 2 | `x<11/5`, `|x'|<3/4` | none |
| 4 | `x<31/10`, `|x'|<21/20` | 1, α |
| 10 | `x≤2α<5`, `|x'|<5/3` | 1, s, α, 2+s, 3+s, 2α |

These are literal two-coefficient checks, not a height cutoff. For
`x=a+b s`, the bounds on `(x−x')/(2s)` give respectively
`b=0`, `b∈{0,1}`, `b∈{0,1}`, and `b∈{0,1,2}`. Intersecting the
remaining integer a intervals gives the displayed lists. The very
slight coarse overcount, where any, is harmless; every exact candidate
is tested below.

The field norms of `9αx²+R` are:

| R | x order as above | field norms |
|---|---|---|
| 4 | 1, α | 7, 439 |
| 10 | 1, s, α, 2+s, 3+s, 2α | 199,136,1279,2296,271,3844 |

A field square has square rational norm. All entries fail except
`3844=62²`; there the full identity, rather than norm alone, is

    9α(2α)²+10=262+180s=(10+9s)².

The positive Z is `10+9s`, congruent to 1 modulo 9O, not to 8. Thus
only the A9 attachment k=1 is integral, and it gives the displayed a.

It remains to exclude the isolated rows. For norm q=1 or 2 use the
native equation `9y²−αx²=q` modulo 3, not the vacuous multiplied
R equation. The quotient O/(3) is F9; `−q/α` has F3 norm 2 and is
therefore not a square in F9. These two rows are empty at every height.

For q=d the bounds give

    x²≤(18/5)dα²=36+(126/5)s <(17/2)²,     |x'|²≤s <(6/5)².

They enclose exactly

    1, α, 2+s, 2α, α², 4+2s, 4+3s.

Indeed b can only be 0,1,2,3; its allowed a values are respectively
`1`; `1,2`; `2,3,4`; `4`. The native equation modulo 3 requires
`x²=−s`, leaving only α and 2α. Their exact tests are

    x=α:       (3y)²=α(α²+s)=3α²,
    x=2α:      (3y)²=α(4α²+s)=30+21s,   norm 18.

The first would put √3 in Q(s), and the second has nonsquare norm.
Both fail. This completes every isolated and connected continuation.

Consequently a is the unique first simple wall of the *full* L9 chamber
on the A8 edge. Its actual full A9 vertex is

    Q=(10+9s;(6+4s)^9),        Q²=−10α.

The retained eight actual walls and a form the full local chamber at Q.
No hidden zero roots or simultaneous alternatives remain: any such
simple wall would have appeared in the complete inventory. This is an
actual edge result, not a claim that a compatible finite prefix is
compact.

## 5. Literature and duplicate-host boundary

The primary low-rank classification
[Bogachev, *From geometry to arithmeticity of compact hyperbolic Coxeter
polytopes*, Theorem C and Table 4](https://nvbogachev.netlify.app/Bogachev-GeomArithmCoxeter-2020.pdf)
includes both `L_3=<-1−√2>+I3` and `M_3=<-√2>+I3` as stably
reflective lattices. The former has a five-wall chamber and the latter
a six-wall chamber. This confirms that low-rank descent is not already
a nonreflectivity obstruction; it supplies no rank-ten classification.
The associated
[2019 thesis summary, Theorem 8 and Table 5](https://www.hse.ru/data/2018/12/29/1142187645/Bogachev-Summary.pdf)
has the same rank-four classification and additional low-dimensional
examples. No exclusion of literal L9 was found in this bounded primary
screen. This is not an exhaustive literature claim.

The earlier excluded literal host is

    S=<-α> ⊥ E8(O) ⊥ <1>.

Its full-arrangement nonreflectivity, with exact infinite-order chamber
symmetry, is in
[the preserved silver E8 proof](fresh_silver_e8_h9_continuation_20261016.md)
and [its scope audit](fresh_silver_e8_chamber_scope_20261016.md).
The positive definite summands `I9` and `E8+<1>` differ integrally, but
this does **not** separate the full indefinite lattices: an integral
isometry may mix time and space. Equal rational forms, determinant
units, and the observed B9 versus E8+A1 initial germs likewise neither
prove nor disprove such an isometry. The dyadic integral identification
is explicitly open. Mere commensurability would not transfer the
literal nonreflectivity conclusion either.

Before any new worker or extension census, settle this duplicate issue
or identify an independent structural payoff from the full-face
identities. In particular Q has norm ratio 10 from p; neither a formal
A9 reversal nor a diagram resemblance alone produces an integral
chamber period. No new compact H9 example or uniform dimension bound
is claimed.
