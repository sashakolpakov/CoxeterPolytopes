# Independent review of the E8 batch, escape and local completion

20 October 2026. **PASS with the scope below.** Full written-source and
returned-certificate inspection; no scientific job, import, enumeration,
Lean or typesetting was performed by this reviewer.

Read in full, with SHA256:

* `compact_h9_sqrt3_e8_integral_splitting_20261020.md`:
  `0e8f62f9354fb3be65738db2bb6b8ffcbfdd864d535a70c10292f71d4210bacb`;
* `compact_h9_sqrt3_e8_first_vinberg_batch_20261020.md`:
  `9d71d939945ca043851ffe53b7eea64c7f7abd53067b94ece6f600c537c1b36b`;
* `compact_h9_sqrt3_twelve_wall_escape_and_completion_20261020.md`:
  `a0633db09d09f5de6bcd9d7b149ca05ed2773e801299e1bec5695e483a178f5f`.

These three files are under `results/`. I also inspected the actual A7
certificate and its generation/replay ledgers, as distinguished below.

## Integral model and actual returned A7 corner

The two quotient columns in the splitting have determinant-1. This proves
an integral E8 plus rank-one decomposition, rather than merely a rational
orthogonal decomposition or an unverified determinant-square identification.
The displayed v is future timelike of square1-s. Parity at eta forces every
normalized reflector's new time coefficient a to be divisible by eta.

The actual certificate
`results/compact_h9_sqrt3_a7_opposite_edge_20261020.json`, SHA256
`226a53992f880c06976bbd9b27ae7dd1ad2ef8b5ca7bbf3312201b6b01d40e76`,
selects exactly the displayed norm2 E8 root d. Its homogeneous endpoint is
9v and has square81(1-s)<0. The two other simultaneous candidates are not
alternative simple walls. Direct coordinate subtraction gives

    A8candidate =3d+(1,3,5,4,3,2,1).(n1,...,n7),
    D8candidate =2d+(1,2,3,2,1,0,0).(n1,...,n7).

All coefficients are nonnegative, so d dominates both in the actual A7
sector. The finite corner in M is E8; in the original chamber it is
E8 plus the orthogonal scalar wall. The ledgers record both generation
and separate replay PASS, in0.248/0.243 seconds and below56MB peak. This
review inspects that evidence and the exact returned algebra; it does not
claim a second independent enumeration.

## Both first-height orbits

I independently checked the coefficient bounds below height5+3s. For q2,
the only possible c=a/eta is1, which would require the unrepresented E8
norm2epsilon. For q2epsilon, c=epsilon is the sole possibility. Writing
the remaining E8 vector X+sY gives

    X^2+3Y^2=42,  X.Y=12,
    Y^2=2k,      k(7-k)>=12.

Thus k=3 or4, with Cauchy equality in both cases. These give respectively
u=epsilon*Y with Y^2=6, and u=(3+2s)*theta with theta^2=2. Coprimality of
2 and3 in the integral E8 module justifies theta's integrality. Omitting
the second case would be a genuine incompleteness; the current source
includes it throughout.

Adjoining beta to A7 gives Schur complement1/8 and inverse diagonal
(4,14,30,20,12,6,2,8). Its positive entries and the stated cross entries
leave exactly omega6 at norm6 and omega7 at norm2 in the dominant cone.
The two first walls P,Q have product-epsilon, hence label3; P's E8 contact
has label12 and Q's is dotted. Their distance-minimality proves their
actual facets and finite perpendicular-foot vertices. Zero-height mirrors
are handled by the E8 chamber inequalities; no positive-height mirror can
cut a shortest perpendicular foot.

The auxiliary codimension-two Gram has determinant1, so the N reduction
really is integral. I checked the projected generators x,y,w, the eta*w
relation, and the determinant-38-22s of D5+x+y. Division by eta^2 gives
epsilon^2(1-s). Keeping w is essential. Finite reflection-group angle
classification excludes an I2(12) reflection pair in an irreducible
rank>=3 component, also for arbitrary local mirrors; this supports the
full orthogonal-arrangement interpretation of N. No nonreflectivity of N
is established by this reduction.

## Twelve-wall escape and eighteen-to-twelve reduction

I recomputed all four new root norms, their E8 contacts and all six mutual
products. They give exactly the listed native labels or ultraparallel
pairs. A small inward perturbation of v is a strict interior point, and
acute projection makes every listed wall an actual facet of this prefix.

The vector X=(5+2s)v-2omega2 has square33-17s>0. Its four extra-wall
products are exactly -2s,-8-2s,0,-2eta. Thus the strict-escape perturbation
and open ideal-patch argument prove infinite volume for the twelve-wall
prefix. They do not prove that the full reflection chamber has this escape.

Deleting alpha2 leaves A1+A6 with the specified order. Finite extensions
give exactly the eighteen necessary profiles, including both endpoint
placements. I checked the inverse row values7 and(10,20,16,12,8,4), the
completed square and its outgoing sign Z<0. The six no-A6-contact rows
reduce modulo7 to norm values3,3,3,5,6,5. Because3 is nonsquare modulo7,
O/(7) is F49; each value is a nonsquare norm in F7. Exactly these six
profiles are excluded at all heights. The twelve others are not excluded.

## Relative-unit completion: precise interpretation

For the displayed bridge at position5, I recovered R=15,K=1 and
r^2=2 from Z^2-14(s-1)a^2=1 with Z=14h+15. The quartic field has
signature(2,1). Dirichlet's unit theorem, passage to the indicated order
and powering modulo14eta are explicit external existence inputs. They
provide infinitely many integral roots with a>0 and eventually h>0.

The limiting vector is the future null endpoint of the escaped edge. It
is strict on all nonwhite prefix walls: the null point lies strictly
between v and the displayed spacelike escape. The limiting argument
therefore makes the nonwhite products ultraparallel for large powers;
white products remain0 or-1. It also preserves any fixed finite set of
future timelike points. This establishes finite endpoint **candidates
compatible with the prefix**, at unbounded height. It does not make all
of them walls of the same full reflection chamber, whose first mirror
could already stop the edge earlier.

No compact H9 example, complete global chamber, or full-lattice
nonreflectivity conclusion follows from these three notes alone. The
prefix escape and the compatible endpoint family are mutually consistent.
