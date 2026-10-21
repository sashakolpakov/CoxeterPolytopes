# Compact Coxeter search map in dimensions 22--29

## Answer

No new compact Coxeter polytope was found in the exact searches archived
here.  The search was concentrated on

\[
\mathbf Q(\sqrt2)\quad\hbox{and}\quad
\mathbf Q(2\cos(2\pi/7)),
\]

with \(\mathbf Q(\sqrt5)\) initially deprioritized and now restored as an
explicit AGCat decision-manifest target.  Most completed searches are still
bounded prefixes: each has an exact hyperideal ray or too few walls for finite
volume.  There is now one genuine global exception.  The silver
\(\Lambda_{22}\) host below is proved nonreflective: its full reflection
subgroup has infinite index and its full reflection chamber has infinite
volume.

The best nontrivial lane is the silver-field host

\[
  \langle-(\sqrt2-1)\rangle\oplus\Lambda_{22}
  \quad\text{in }\mathbf H^{22}.
\]

Its earlier computation supplied a complete first packet, four explicitly
accepted later walls, and strict replacement rays.  More decisively, an
explicit infinite-order integral isometry preserves a chamber of the full
mirror arrangement.  A complete separator bound reduces that global claim to
the already exact \(H\le150\) frontier.  The resulting infinite cyclic chamber
stabilizer proves infinite chamber volume, so this exact host cannot contain a
finite-covolume integral reflection subgroup.

This document maps the **host families actually searched**.  It is not the
complete finite classification promised by the cycle-packet \(\to\)
reflection-order \(\to\) conductor-lattice architecture: that global packet
and conductor enumeration has not yet been materialized in ranks 23--30.
Consequently, failure of the hosts below does not empty a dimension or a
field.

## Eight-rank coverage matrix

| hyperbolic dimension | \(\mathbf Q(\sqrt2)\) | \(\mathbf Q(\sqrt5)\) | \(\mathbf Q(2\cos(2\pi/7))\) |
|---:|---|---|---|
| 22 | global nonreflectivity theorem for \(\langle-(\sqrt2-1)\rangle\oplus\Lambda_{22}\); separate prime-seven coefficient first packet and bounded non-orthogonal hosts | early \(D_{22}\) norm-two screen; full-rank \(c=2+2\varphi\) path host complete through \(H=10\), with one new wall and 41 hyperideal rays | diagonal two-seed screen; \(\Lambda_{22}\) through \(H=100\); \(E_8^2E_6\) screen |
| 23 | determinant-three odd-Leech section first packets for two coefficients | determinant-three odd-Leech section first packets for \(\varphi,2\varphi\) | diagonal two-seed screen; determinant-three odd-Leech section first packet |
| 24 | Leech, \(K_{12}^2\), \(A_1^{24}\)/Golay, \(D_{24}^{+}\), \(E_8^3\), \(N(A_2^{12})\), \(N(E_6^4)\), \(N(A_{11}D_7E_6)\), \(N(A_{12}^2)\), \(N(A_{15}D_9)\), \(N(D_{10}E_7^2)\), and \(N(D_{12}^2)\) bounded screens; \(K_{12}^2\) continued through its first later cutter at \(H=393+278\sqrt2\) with a strict replacement | early \(D_{24}/D_{24}^{+}\) norm-two screen only | diagonal two-seed screen; rootless Leech two-packet chronology and all-wall obstruction through \(H=1000\); \(D_{24}^{+}\) through \(H=1000\); \(E_8^3\) screen |
| 25 | no dedicated host completion | early \(D_{25}\) norm-two screen only | diagonal two-seed screen |
| 26 | rootless odd-unimodular \(S_{26}\) host through real \(H=1000\); exact full-shell-polar first packet for the determinant-three \(T_{26}\) host | early \(D_{26}\) norm-two screen only | diagonal two-seed screen |
| 27 | rootless odd-unimodular \(T_{27}\) host: exact 119,574-wall first packet and all-wall obstruction through real \(H=1000\) | early \(D_{27}\) norm-two screen only | diagonal two-seed screen |
| 28 | rootless odd-unimodular \(U_{28}\) host through real \(H=1000\); first accepted cutters and strict replacements on the 3-modular `DualExtremal(28,3)` and 5-modular `J2dim14` hosts | early \(D_{28}\) norm-two screen only | diagonal two-seed screen |
| 29 | rootless odd-unimodular \(\mathrm{dim29odd}\) host: exact 93,786-wall first packet and all-wall obstruction through real \(H=1000\) | early \(D_{29}\) norm-two screen only | diagonal two-seed screen |

Here “screen” always means an exact bounded Vinberg prefix for the displayed
host, not an exhaustive field/rank decision.

## Field-by-dimension map

| field | dimensions tested | exact host/family | certified search | outcome at cutoff |
|---|---:|---|---|---|
| \(\mathbf Q(\sqrt2)\) | 22 | \(\langle-(\sqrt2-1)\rangle\oplus\Lambda_{22}\) | complete first packet and exact \(H\le150\) separator frontier; explicit infinite-order integral isometry stabilizing a full reflection chamber | **global host-level nonreflectivity:** the full chamber has infinite volume and the reflection subgroup has infinite index; hence no integral reflection subgroup in this host has finite covolume |
| \(\mathbf Q(\sqrt2)\) | 22 | non-orthogonal rank-23 arithmetic chains/trees with edge product \(2+\sqrt2\) | complete 66-attachment \(A_{22}/D_{22}/B_{22}\) seed census; chain continued through its first two simultaneous norm-two walls | best seed is 3 finite/0 cusp/20 hyperideal; the 25-facet first-cut cone is 6/0/97, hence not finite volume |
| \(\mathbf Q(\sqrt2)\) | 22 | \(\langle-(2\sqrt2-1)\rangle\oplus\Lambda_{22}\) | complete local alphabet and first packet at \(H=5+4\sqrt2\) | 49,896 walls, but a strict hyperideal witness survives |
| \(\mathbf Q(\sqrt2)\) | 23 | \(\langle-a\rangle\oplus S_{23}\), \(a=1+\sqrt2,\sqrt2\), where \(S_{23}\) is an odd-Leech section | complete local alphabets, height-zero exclusions, and first packets | 2,576 and 3,742,928 primitive walls respectively; strict hyperideal witnesses survive |
| \(\mathbf Q(\sqrt2)\) | 24 | \(\langle-(2\sqrt2-1)\rangle\oplus D_{24}^{+}\) | glue-complete through \(H=500\), 4,374 walls | exact rank-24 hyperideal extreme ray |
| \(\mathbf Q(\sqrt2)\) | 24 | \(\langle-(12\sqrt2-13)\rangle\oplus K_{12}^{2}\), with \(K_{12}\) the rootless Coxeter--Todd lattice | complete 16-class local alphabet; 128 controls at 87 heights through the first later cutting packet \(H=393+278\sqrt2\) | the 613,872-wall first packet and both intermediate necessary levels preserve the old ray; an explicit accepted norm-two cutter then appears, but \((1,(7/4)h,0)\) is a strict hyperideal replacement through the cutoff |
| \(\mathbf Q(\sqrt2)\) | 24 | \(\langle-\sqrt2\rangle\oplus N(A_2^{12})\) and \(\langle-2\sqrt2\rangle\oplus N(A_2^{12})\) | extended ternary-Golay candidate frontier complete through real \(H=200\); eight nonempty necessary levels | strict uniform hyperideal points survive every complete necessary superset |
| \(\mathbf Q(\sqrt2)\) | 24 | \(\langle-(1+\sqrt2)\rangle\oplus N(E_6^4)\) | tetracode-glue complete through \(H=15\); packets of 4 and 44 walls | 24 glue-coupled walls occur, but a strict hyperideal point survives |
| \(\mathbf Q(\sqrt2)\) | 24 | \(\langle-(1+\sqrt2)\rangle\oplus N(A_{15}D_9)\) | cyclic index-eight glue complete through \(H=50\); nonempty packets of 2, 15, and 191 walls | 120 cross-component glue walls occur, but a strict D9-axis hyperideal ray survives |
| \(\mathbf Q(\sqrt2)\) | 24 | \(\langle-(2\sqrt2-1)\rangle\oplus N(A_{15}D_9)\) | four complete local norm classes and cyclic glue through \(H=63\); 1,666 positive-height walls | dense 1,655-wall second packet is acute, but a strict D9-axis hyperideal ray survives |
| \(\mathbf Q(\sqrt2)\) | 24 | \(\langle-(1+\sqrt2)\rangle\oplus N(D_{12}^2)\), plus coefficients \(\sqrt2,2\sqrt2\) | exact vector--spinor graph glue; unit/\(\sqrt2\)/\(2\sqrt2\) prefixes through real heights 30/40/27 with 99/19/27 accepted walls | all three have strict sparse hyperideal rays; the unit branch contains 37 cross-glue walls |
| \(\mathbf Q(\sqrt2)\) | 24 | \(\langle-(1+\sqrt2)\rangle\oplus N(A_{12}^2)\) | exact prime-13 graph glue complete through real \(H=150\), 280 accepted walls | all walls through the cutoff have zero glue residue; a strict node-one hyperideal ray survives |
| \(\mathbf Q(\sqrt2)\) | 24 | \(\langle-(1+\sqrt2)\rangle\oplus N(A_{11}D_7E_6)\) | cyclic index-12 glue complete through the first two nonempty packets; 3 and 35 walls | 22 second-packet walls use glue and 12 couple all three factors, but a strict \(A_{11}\)-axis hyperideal ray survives |
| \(\mathbf Q(\sqrt2)\) | 24 | \(\langle-(1+\sqrt2)\rangle\oplus N(D_{10}E_7^2)\) | binary index-four glue complete through the first two nonempty packets; 3 and 16 walls | seven second-packet walls use glue and three couple all factors, but a strict \(D_{10}\)-axis hyperideal ray survives |
| \(\mathbf Q(\sqrt2)\) | 24 | Leech, \(A_1^{24}\)/Golay, \(D_{24}^{+}\), and \(E_8^3\) spatial hosts | coefficient-dependent exact packet/frontier searches | every completed prefix has an explicit hyperideal witness |
| \(\mathbf Q(\sqrt2)\) | 24 | \(\langle-(-13+14\sqrt2)\rangle\oplus\Lambda_{24}\) (norm-223 coefficient) | complete locally gated necessary census through rational \(H=250000\); 1,129 controls and 2,855 surviving Gram types; limiting norm-two root at \(H=114998+81316\sqrt2\) independently realized and Vinberg-accepted | a strict one-direction hyperideal witness still survives with excess about \(1.68\times10^{-9}\); later completion remains open |
| \(\mathbf Q(\sqrt2)\) | 26 | \(\langle-(2\sqrt2-1)\rangle\oplus S_{26}\), with \(S_{26}\) the rootless Conway--Borcherds odd unimodular lattice | complete six-class local alphabet, exact 102,180-wall first packet, and every possible wall through real \(H=1000\) | a fixed minimal direction with scale \(43/55\) remains strictly hyperideal against all 256 necessary Gram types |
| \(\mathbf Q(\sqrt2)\) | 27--29 | \(\langle-(2\sqrt2-1)\rangle\oplus M\), for the rootless odd-unimodular catalogue hosts \(T_{27},U_{28},\mathrm{dim29odd}\) | primary Grams independently validated; complete six-class alphabet and no height-zero walls; exact first packets of 119,574, 98,280, and 93,786 walls; every possible wall through real \(H=1000\) | the same rank-independent \((1,(43/55)d)\), \(d^2=3\), remains strictly hyperideal against all 256 necessary Gram types |
| \(\mathbf Q(\sqrt2)\) | 26 | \(\langle-(2\sqrt2-1)\rangle\oplus T_{26}\) | complete first packet: all 117,936 minimum-shell walls; exact 66-contact polar vertex | a strict shrink of the polar vertex remains hyperideal; later walls are open |
| \(\mathbf Q(\sqrt2)\) | 28 | \(\langle-(2\sqrt2-1)\rangle\oplus\mathrm{DualExtremal}(28,3)\) | complete eight-class chronology through the first accepted cutter at \(H=362+256\sqrt2\); 103 controls, six nonempty levels, and all 17 earlier cross-Gram bounds independently hostile-audited | the cutter removes \((1,d/2)\), but \((1,(12/25)d)\) is a strict hyperideal replacement through the cutoff |
| \(\mathbf Q(\sqrt2)\) | 28 | \(\langle-(1+\sqrt2)\rangle\oplus J2\mathrm{dim}14\) | complete 374,976-wall first packet and the first accepted norm-two cutter at \(H=137+97\sqrt2\); eight coefficient lanes screened through real \(H=300\) | the cutter removes \((1,(3/2)G^{-1}e_{13})\), but \((1,(5/4)G^{-1}e_{13})\) is a strict hyperideal replacement through the cutoff |
| \(\mathbf Q(2\cos(2\pi/7))\) | 22--29 | \(\operatorname{diag}(-\theta,1^{d-1},t)\), \(t=2,2-\theta\) | complete height-5 geometry in every dimension; root patterns stabilize through height 50 | every cone has hyperideal rays; no closure |
| \(\mathbf Q(2\cos(2\pi/7))\) | 22 | non-orthogonal rank-23 chain with coupling \(c=2+\theta\) | exact coefficient optimality near the Lorentz threshold; selected signature \((22,1)\), both conjugates positive definite, genuinely cubic minimal field | new anisotropic/cocompact arithmetic host, but seed cone is 3 finite/0 cusp/20 hyperideal |
| \(\mathbf Q(2\cos(2\pi/7))\) | 23 | \(\langle-\theta\rangle\oplus S_{23}\), where \(S_{23}\) is an odd-Leech section | complete local alphabet and 2,576-wall first packet | a strict norm-four section witness survives |
| \(\mathbf Q(2\cos(2\pi/7))\) | 22 | \(\langle-\theta\rangle\oplus\Lambda_{22}\) | complete local norm/gate and Gram scan through \(H=100\) | no wall at all through the cutoff; no chamber closure |
| \(\mathbf Q(2\cos(2\pi/7))\) | 22, 24 | \(E_8^2E_6\) and \(E_8^3\) spatial blocks | complete through \(H=10\) | the only nonempty early packet leaves hundreds of hyperideal rays |
| \(\mathbf Q(2\cos(2\pi/7))\) | 24 | \(\langle-2\theta\rangle\oplus\Lambda_{24}\) | complete first two packets; every possible later wall through \(H=1000\) bounded by exact augmented Gram minima | packets of 196,560 and 9,258,762,240 accepted walls; a strict norm-six hyperideal point survives all 18 controls |
| \(\mathbf Q(2\cos(2\pi/7))\) | 24 | Leech coefficients \(\theta,\theta(2-\theta),2\theta(2-\theta)\) | complete mirror norms/gates and necessary Gram prefixes through \(H=300\) | respectively one, zero, and two nonempty necessary levels; strict hyperideal witnesses in every case |
| \(\mathbf Q(2\cos(2\pi/7))\) | 24 | \(\langle-\theta(2-\theta)\rangle\oplus D_{24}^{+}\) | glue-complete through \(H=1000\), 179 walls | old rays are cut, but an exact replacement hyperideal extreme ray survives |
| \(\mathbf Q(\sqrt5)\) | 23 | \(\langle-a\rangle\oplus S_{23}\), \(a=\varphi,2\varphi\), where \(S_{23}\) is an odd-Leech section | complete norm alphabets, height-zero exclusions, and first packets; independently reconstructed from Golay/Construction B | 553,104 and 53,682 walls respectively; both exact prefixes retain strict hyperideal points |
| \(\mathbf Q(\sqrt5)\) | 22 | full-rank terminal-dotted path host with \(c=2+2\varphi\) | complete 16-class primitive norm alphabet and all 383 controls through \(H\leq10\); one new post-seed norm-two wall | the true 24-wall cone has 2 finite vertices, 0 cusps, and 41 hyperideal rays; later completion and reflection index remain open |
| \(\mathbf Q(\sqrt5)\) | 22--29 | universal AGCat frame packets | complete dotted three-node tranche; connected four-node finite/orthogonal, exactly-one-dotted, and exactly-two-dotted-tree tranches; explicit full-rank terminal-dotted path family | 2,633,885 three-node packets, 180,521 at-most-one-dotted four-node packets, and 97,172,716 two-dotted trees survive; 3,590 down to 3,552 path products give cocompact arithmetic host/frame seeds in ranks 23--30, but reflectivity and global cyclic joins remain open |
| \(\mathbf Q(\sqrt5)\) | 22--29 | early \(D_n/D_{24}^{+}\) screen | norm-two controls through \(\varphi^2\) | hyperideal count grows with dimension; incomplete as a full mirror-norm search and not promoted |

The diagonal cubic-seven search is the only single family currently mapped
uniformly across all eight dimensions.  At height five its exact
finite/cusp/hyperideal ray counts are

\[
  (3,0,d-2)\quad(t=2),\qquad
  (5,0,2d-5)\quad(t=2-\theta),
\]

for every \(22\le d\le29\).  At height 50 the complete root patterns at
the endpoints \(d=22,29\) agree by zero padding, which certifies the same
bounded pattern in all intermediate dimensions.  It does not prove that
later walls cannot close a chamber.

## Strongest lane: silver \(\Lambda_{22}\) in \(\mathbf H^{22}\)

The spatial lattice has determinant 12, minimum 4, kissing number 49,896,
and Smith form \(1^{20},2,6\).  The complete primitive mirror-norm alphabet
is

\[
1,\ 2,\ 4+2\sqrt2,\ 4,\ 3,\ 6,\ 12+6\sqrt2,\ 12.
\]

The first complete simultaneous packet occurs at
\(H=14+10\sqrt2\).  It contains

- 467,596,800 norm-two roots from the ordinary \(\Lambda_{22}\) norm-ten
  shell; and
- 1,368,576 norm-four roots from the dyadic radical's norm-twenty shell.

All 468,965,376 roots are mutually Vinberg-acute.  The exact ray

\[
  z=(1,d/3),\qquad d^2=4,
\]

has norm \(13/9-\sqrt2>0\) and is strict against the complete packet.

The first explicit accepted wall cutting this ray is at
\(H=164+116\sqrt2\).  With \(d=e_1\),

\[
e=(-6,31,21,15,10,20,14,7,4,2,2,1,0^{10}),
\]

put \(u=8d-e\), \(v=6d-e\).  Then

\[
r_{164}=(28+20\sqrt2,u+\sqrt2v),\qquad r_{164}^2=4.
\]

It is primitive, reflection-integral, and chronologically accepted.  It
cuts \(d/3\), but the replacement \((1,13d/40)\) remains hyperideal.

The next accepted cutter found is

\[
r_{738}=(42+30\sqrt2,14d+g+\sqrt2(10d+g)),
\quad g=(0,-1,0^{19},-1),
\]

with norm two and height \(738+522\sqrt2\).  An exact chronological Schur
scan checks 195 controls and 517 gate-compatible necessary types.  All
earlier positive threats are arithmetically empty; the only remaining
positive extension is the displayed wall itself.  It cuts \(13d/40\).

Nevertheless,

\[
  z''=(1,(32/99)d),\qquad
  (z'',z'')=13897/9801-\sqrt2>0,
\]

is strict against every possible wall through the second cutter.  A third
accepted norm-two wall occurs at

\[
 H=1598+1130\sqrt2,
 \qquad
 r=(62+44\sqrt2,(20d+e)+\sqrt2(14d+e)),
 \quad e=-e_3-e_{17}.
\]

The exact chronological audit covers 284 controls and 1,134 gated necessary
types.  This wall cuts \(32d/99\), but the new replacement

\[
 (1,(10/31)d),\qquad 1361/961-\sqrt2>0,
\]

is strict against every possible wall through the new cutoff.  The same
replacement was subsequently certified through real height \(5000\): 69 new
controls and 635 gate-compatible necessary types leave seven relaxed threats,
all excluded by the ambient minimum four or the dyadic-radical minimum twelve.

A fourth accepted norm-two wall occurs at

\[
 H=2786+1970\sqrt2,
 \qquad r=(82+58\sqrt2,(26d+e)+\sqrt2(18d+e)),
 \qquad e=e_2+e_{22}.
\]

The exact chronology covers 374 controls and 1,997 gated post-first-packet
types.  Its sole new earlier positive row is empty by a finite 13,824-vector
mod-nine coset exhaustion, so the wall is accepted.  It cuts `10d/31`, but

\[
 z^{(4)}=(1,(28/87)d),\qquad
 (z^{(4)},z^{(4)})=10705/7569-\sqrt2>0,
\]

is strict against every possible wall through the fourth cutter.  Continuing
the complete necessary scan to real height 10000 adds 129 controls and 1,555
gated types.  Its unique relaxed threat has 665 integral PSD extensions left
after the explicit `[-3,3]^3` short-combination gate (every discarded case
already contradicts the ambient minimum four); exhaustive enumeration finds exactly one nonnegative
extension.  That extension forces a nonzero \(\Lambda_{22}\) vector of norm
two, hence is empty.  These bounded computations are now finite lemmas inside
a global argument.  For `v0=(1,0)`, the integral isometry

\[
 g(x,y)=\bigl((5+4\sqrt2)x+(4+3\sqrt2)(d,y),
 y+((2+\sqrt2)x+(1+\sqrt2)(d,y))d\bigr)
\]

has infinite order.  Any mirror separating `v0` from `g*v0` has height
at most `56+40sqrt(2)<150` and pairs positively with `z=(1,d/3)`, which
the complete H14/H150 certificates rule out.  Thus `g` stabilizes the full
reflection chamber containing `v0`.  Infinitely many disjoint translates
inside that chamber force infinite volume and infinite reflection index.
Therefore this exact arithmetic host is globally nonreflective; the scope is
host-level, not all of \(\mathbf Q(\sqrt2)\) or all of \(\mathbf H^{22}\).

The prime-seven silver coefficient gives a useful independent \(H^{22}\)
check on the same spatial lattice.  For

\[
  \langle-(2\sqrt2-1)\rangle\oplus\Lambda_{22},
\]

the complete local mirror alphabet has 16 classes.  All 13 controls before
\(H=5+4\sqrt2\) are empty.  At that height, imprimitive duplicates removed,
the complete packet is exactly the 49,896 walls

\[
 r_d=(2+\sqrt2,(1+\sqrt2)d),\qquad d^2=4.
\]

For the exact dual-lattice vector \(k\) with \(Gk=3e_{22}\) and \(k^2=24\),

\[
 z=(1,(4-\sqrt2)k/8),\qquad z^2=31/4-5\sqrt2>0,
\]

and its maximum pairing with the full packet is
\(-(2+3\sqrt2)/4<0\).  This complete first packet therefore also fails the
finite-volume test.

## Rootless 3-modular silver host

The non-Niemeier spatial lattice (K_{12}^{2}), a direct sum of two
Coxeter--Todd lattices, gives a genuinely different silver (H^{24}) lane:

\[
  \langle-(12\sqrt2-13)\rangle\oplus K_{12}^{2}.
\]

The coefficient has norm (-119=-7\cdot17) and is unit-square-equivalent
to one archived norm-119 Leech direction, but the rootless 3-modular spatial
lattice and its ternary reflection gates are new.  The complete local
alphabet has sixteen classes.  All height-zero classes and all 48 controls
before (H=67+48\sqrt2) are empty.  The first packet has 1,512 norm-one
walls from (K_{12}^{2}(4)) and 612,360 norm-two walls from
(K_{12}^{2}(8)), for 613,872 walls total.

Let (h=G_{K_{12}}^{-1}e_1), so (h^2=4/3) and
(h\cdot K_{12}\subset\mathbf Z).  The point
(z=(1,(9/5)h,0)) has norm (433/25-12\sqrt2>0).  Its maximum pairings
with the two packet components are (9/5-(14/5)\sqrt2) and
(-19/5+(9/5)\sqrt2), both negative.  Hence the complete first-packet
prefix has an open hyperideal patch and infinite volume.

The exact continuation contains 128 controls at 87 heights through
(H=393+278\sqrt2).  The fixed-coordinate identity

\[
 \min\{x^2:x_1=k\}=3k^2/4+\rho(k\bmod4),\qquad
 \rho=(0,13/4,1,13/4),
\]

proves that no earlier possible wall cuts (z).  At the endpoint the explicit
primitive norm-two root with coefficient Gram (380,204,278) and
((h\cdot U,h\cdot V)=(22,16)) is acute to the entire first packet, every
intermediate necessary Gram type, and every simultaneous wall, while its
pairing with (z) is (8/5-\sqrt2/5>0).  The replacement
(z'=(1,(7/4)h,0)) has norm (205/12-12\sqrt2>0) and strict endpoint margin
(1/2-\sqrt2<0). Thus the complete prefix through the first cutter still has
infinite volume. This remains a bounded obstruction, not a global
nonreflectivity result. A clean-room hostile replay independently verifies the
chronology, cutter, same-height lemma, and replacement margins.

## Rootless odd-unimodular silver family in H27--H29

The same prime-seven coefficient gives a rank-independent theorem for any
integral unimodular spatial lattice of minimum three with a norm-three vector.
Applying it to the literal Nebe--Sloane catalogue Grams for
\(T_{27},U_{28}\), and \(\mathrm{dim29odd}\) produces genuinely new
target-rank hosts. Exact fraction-free determinant replay gives determinant
one, and PARI gives minimum three together with norm-three shells of 1,640,
2,240, and 1,856. The first packet at \(H=5+4\sqrt2\) is the complete
norm-four shell, with respectively 119,574, 98,280, and 93,786 walls.

The local proof has six primitive mirror-norm classes and no height-zero
walls. The uniform census through real \(H=1000\) has 134 controls at 104
heights and 256 retained necessary binary Gram types. For any fixed
\(d^2=3\), the point \((1,(43/55)d)\) is strict against all of them. Its
tight row has height \(455+322\sqrt2\), maximum pairing
\(37+26\sqrt2\), margin \((-51+37\sqrt2)/55>0\), and Lorentz excess
\(8572/3025-2\sqrt2>0\). Thus all three bounded prefixes have infinite
volume; the statement is neither Vinberg termination nor global
nonreflectivity.

## Cubic-seven controls

Two complementary cubic-seven screens are useful.

First, \(\langle-\theta\rangle\oplus\Lambda_{22}\) has complete mirror
norms \(1,2,4,3,6,12\).  Thirty exact controls at 24 heights and 89 ambient
Gram types reduce to eleven gate-compatible types, all excluded by the
radical minima.  Hence there is no reflecting wall through \(H=100\).
That is a rootless bounded prefix, not a compactness or nonreflectivity
result.

Second, the nonunimodular \(D_{24}^{+}\) host

\[
  \langle-\theta(2-\theta)\rangle\oplus D_{24}^{+}
\]

has 179 walls through \(H=1000\).  Both new 77-wall packets cut every old
hyperideal ray, but the exact extreme ray

\[
\left((23+5\theta-3\theta^2)/14,1,1,1,0^{21}\right)
\]

has positive Lorentz norm, active rank 24, and pairing census
\(154/25/0\) (negative/zero/positive).  It therefore remains noncompact at
the cutoff.

The rootless-Leech host

\[
\langle-2\theta\rangle\oplus\Lambda_{24}
\]

has now advanced beyond its formerly empty `H=10` shell.  Its first two
complete simultaneous packets occur at heights
`5+14theta+6theta^2` and `16+45theta+20theta^2`, with respectively
196,560 and 9,258,762,240 accepted walls.  The second count is the exact
ordered incidence count for minimal Leech pairs of inner product `-1`, and
an exhaustive four-vector Gram comparison gives same-height maximum pairing
`-1`.  Nevertheless, for an explicit norm-six Leech vector `d`, the point
`(1,2d/3)` is hyperideal and strict against all 47 necessary Gram types in
all 18 controls through `H=1000`.  This rules out finite volume only for the
bounded prefix, not for the eventual chamber.

The genuinely new rootless host

\[
\langle-2\theta\rangle\oplus K_{12}^{2}
\]

has a much smaller but still multidirectional continuation.  Its first two
complete packets contain 1,512 and 241,920 walls, respectively: the minimal
shells in either factor, then ordered minimal K12 pairs of inner product
`-1`.  Exact Gram comparison proves both packets simultaneously acute.  A
19-squareclass support screen through `H=100` leaves only the inert-2,
inert-2/ramified-7, and inert-2/discriminant-3 lanes necessary-nonempty; the
other 16 are rigorously empty in that window.  This does not yet settle the
volume of the selected lane.  The first later wall cutting
`(1,(3/2)G_K12^-1 e1,0)` is now explicit and accepted at
`H=33+94theta+42theta^2`.  A hostile replay exposed an omitted necessary H33
Gram type in the original `36/25` proof; the repaired `7/5` replacement is
strict against both H33 types and through that same bound, so the bounded
prefix still has infinite volume.

## Decision and next target

There is no result here that raises the compact-existence lower bound beyond
8, and no field/dimension slice has been emptied.  There is, however, now a
global nonexistence theorem for the single silver \(\Lambda_{22}\) host: it
cannot support a finite-covolume integral reflection group.

The priority order supported by the data is:

1. retire the globally nonreflective silver \(\Lambda_{22}\) host rather than
   extending its Vinberg cutoff;
2. try a genuinely different rootless or low-root spatial host over
   \(\mathbf Q(\sqrt2)\) or the cubic-seven field, rather than another dense
   \(D_{24}^{+}\)/Leech variant;
3. retain the cubic-seven \(\Lambda_{22}\) rootless prefix as a control, but
   do not spend on full cone geometry until the first wall is reached;
4. for \(\mathbf Q(\sqrt5)\), use the AGCat cycle-packet/order/conductor
   architecture to enumerate the full rank-23--30 slice; isolated host
   failures do not prove nonexistence.

## Primary local records

- `results/silver_lambda22_h14.md`
- `results/silver_lambda22_h150_frontier.md`
- `results/silver_lambda22_h164_cutter.md`
- `results/silver_lambda22_h738_cutter.md`
- `results/silver_lambda22_h1598_cutter.md`
- `results/silver_lambda22_h5000_replacement.md`
- `results/silver_lambda22_h5000_replacement_independent_hostile_audit.md`
- `results/silver_lambda22_global_nonreflectivity.md`
- `results/silver_lambda22_global_nonreflectivity_independent_hostile_audit.md`
- `results/silver_p7_lambda22_first_packet.md`
- `results/silver_coxeter_todd_squared_first_packet.md`
- `results/silver_coxeter_todd_squared_h393_replacement.md`
- `results/silver_coxeter_todd_squared_h393_independent_hostile_audit.md`
- `results/silver_s26_h26_first_packet.md`
- `results/silver_odd_unimodular_h27_h29_h1000.md`
- `results/silver_dual_extremal_28_3_h362_replacement.md`
- `results/silver_dual_extremal_28_3_h362_replacement_independent_hostile_audit.md`
- `results/cubic7_coxeter_todd_squared_h100.md`
- `results/cubic7_coxeter_todd_squared_h33_cutter.md`
- `results/silver_a2_ternary_golay_h24_h17.md`
- `results/silver_a2_ternary_golay_h24_cutter.md`
- `results/silver_a2_ternary_golay_h200.md`
- `results/silver_niemeier_unscreened_triage_20260812.md`
- `results/silver_d12_squared_h24.md`
- `results/silver_a12_squared_h24_h150.md`
- `results/silver_a11_d7_e6_h24_h15.md`
- `results/silver_d10_e7_squared_h24_h15.md`
- `results/silver_odd_leech_section_h23_first_packets.md`
- `results/golden_odd_leech_section_h23_first_packets.md`
- `results/golden_odd_leech_section_h23_first_packets_independent_hostile_audit.md`
- `results/golden_path_h22_h10_prefix.md`
- `results/cubic7_odd_leech_section_h23_first_packet.md`
- `results/cubic7_lambda22_h100.md`
- `results/cubic7_bna1_h22_adaptive_search.md`
- `results/cubic7_bna1_rank6_and_target_prefixes.md`
- `results/silver_p7_d24_plus_h500.md`
- `results/cubic7_p7_d24_plus_h1000.md`
- `results/cubic7_leech_h24_h1000_continuation.md`
- `results/golden_dn_h22_h29_early_screen.md`
- `results/golden_agcat_h22_h29_decision_manifest.md`
- `results/golden_agcat_four_node_at_most_one_dotted.md`
- `results/golden_agcat_four_node_at_most_one_dotted_hostile_audit.md`
- `results/golden_agcat_two_dotted_four_node_trees.md`
- `results/golden_agcat_two_dotted_four_node_trees_independent_hostile_audit.md`
