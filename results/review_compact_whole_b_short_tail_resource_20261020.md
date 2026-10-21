# Outside audit of the whole-B short-tail extension

20 October 2026. **PASS**, with the stated whole-vertex scope.
Full independent mathematical read of
`results/compact_whole_b_short_tail_extension_20261020.md`, SHA256
`b6adc03a5950a60a1e263aeb79869785268f1bb7905e74ec47a784446c989c5c`.
I did not coauthor this source. No scientific execution, worker import,
optimizer, Lean or TeX job, accepted-source edit, or status change.

I also read the unchanged inherited whole-B proof
`compact_whole_b_bulk_weak_packet_and_endpoint_pressure_20261020.md`,
SHA256 `e5a3c51858bb7b1e392b2828c874387be8d2044d5ae1be583a546b67eaddb54c`,
and the relevant isolated/actual-point arguments in the whole-A21 source,
SHA256 `3abe11bb70ce677069279462a91111ed6b8ab8889e3ec9e656ac71c1094200fb`.
The following B-specific extensions and final margins were rederived.

## Complete tail inventory and normalization

The unit convention is consistent throughout: the old long normals are
`(e_i-e_(i+1))/sqrt2`, the terminal normal is e_n, the long duals are
sqrt2 times prefixes, and the short dual is the unscaled whole prefix.
In particular a new simple3 contact at the old SHORT endpoint contributes
1/(2sqrt2), not1/2, to the released long-direction sigma.

At n-3, the retained B3 admits just the stated B4 long-end extension
and F4 short-end extension. A contact to the left component would extend
the entire F4 subdiagram and is impossible. With no left contact, own
zero gives the exact affine F4 path, and every positive own contact has
delta1/4 and sigma at least1/2+1/(2sqrt2); the strong bound follows.

At n-2, B2 admits exactly the two simple3 B3 extensions. For the new
short-end contact without a left contact, I obtain delta1/2 and the
same sigma correction. Own zero gives the displayed affine F4 path;
own4 gives the exact 4,3,3 star and own6 the affine6,3 path. Every
own label at least5, or dotted value, is strong. Only own3 remains.

If the new short-end branch also touches the left A path, the finite
connected white diagram has rank n>=23, contains the terminal4, and
must be B_n. The left contact is simple3 at one endpoint, with no
additional contact. I obtain delta=1/(2i) and the source's lower sigma.
This bound is strong EVEN WHEN THE OWN CONTACT IS ZERO. Thus no
zero-own branch was incorrectly discarded as an affine cycle.

The left-white arguments continue unchanged. I checked the potentially
delicate isolated5 boundary: the selected old A3 lies on the left,
and its unselected right dual omega_(i+1) remains long through i=n-2.
The displayed 3/4,1/4 interpolation is therefore valid. Isolated7
needs only a simple old neighbor; the >=8 and dotted estimates use
the old released node itself. The isolated3 terminal-4 affine fork,
isolated4 double-terminal4 path, and isolated6 affine subpath also
remain available in both new layers.

## Actual mate and helper inequalities

For the additional profile v=(omega_i+omega_n)/2, deleting n-1 gives
A_n and deleting n gives D_n. These are complete positive-definite
sets of ACTUAL normals. Positive acute projection gives the actual
endpoints; no merely abstract finite completion is substituted.

I independently obtain the mate residual (i+3)/(4(i+1)), its sigma,
and C=[(1+sqrt2)i+1]/sqrt(i+3). The polynomial proving C/q>79/50
has value17152 at i21 and increases thereafter. With C>9 and q<7,
the displayed bound using sqrt(C^2-1)>C-1/C gives a strictly strong
n-1 direction for every21<=i<=23.

All four actual-pressure constructions for the far packet remain
available through i=n-2. Concavity gives x1>sqrt2 times the displayed
bracket with the correct direction. For finite labels I checked
A-gamma*q>.079, using gamma>101/125 and A+gamma*q<14, and every
native/root interval in the four separate5,7,8,>=9 rows. Even the
weakest finite row exceeds.424 after the conservative sqrt2>7/5.

For dotted gamma, I rederived

    dA/dgamma-q=(1+q^2/2)/[A(q^2 gamma+1+qA)]
                 <(q^2+2)/(4q^3)<1/10.

The foot-root derivative is greater than1 on(1,6/5), so the bracket
decreases with gamma; its q derivative is negative because A>gamma*q.
The endpoint reduction is therefore global. Precisely the written
rational radical bounds give the lower value

    1.4142*(7.8981-1.2*6.480741+.732-.552771)
      =.42488196516>.424.

This verifies the tight helper demand without a sampled extremum.

## Prefix norm, holes and complete rank split

The accepted nonadjacency proof for two far packets also applies when
one is at n-2: its full contact form, entries C_ij=2min(i,j), and the
common LEFT helper projection are unchanged. The extra terminal
profile is treated separately and is never silently included in that
far-packet prohibition.

For n25 one of21,22 is strong. For n24, either21 or22 is strong,
or22 has the terminal profile and23 is strong; adjacent far21/22
is impossible. For n23 the same reference argument handles strong21
or the terminal profile at21, which makes22 strong. There can be
only one consecutive two-hole cluster, from the terminal profile and
its preceding far packet. Earlier far holes are separated from it
and from one another. The stated deliberately excessive loss
9/3200+1/800 is below.009.

The curvature and loss formulas check: an isolated omitted knot loses
less than1/3200, and two consecutive knots lose less than1/800 by
the variance of three slopes with successive gaps below1/40.
Disjoint cluster segment interiors make these losses additive.
The surviving reference stays concave with nonnegative dual weights.

In the remaining n23 case, far-weak21 makes20 strong and gives the
new x1>.424. There are at most seven nonadjacent holes in6..19.
Fixing the first coordinate changes the first five slopes to
.424,(f5-.424)/4 repeated four times. They remain positive and
strictly decrease across1 and5, including after any hole beginning
at6. The exact energy increment is (5/4)(.424-f5/5)^2. None of
these changed slopes participates in the later hole losses.

I checked E20>.99752, f5/5<.361, and the final rational sum exactly:

    .99752+(5/4)(.424-.361)^2-7/3200=1.00029375>1.

This completes n23..25; the immutable accepted n>=26 result supplies
the rest. I found no gap or correction to request. The conclusion
excludes WHOLE B_n farthest vertices for n>=23. It does not assert
the same rate for B components with spectators, a result below23,
an arithmetic-host exclusion, or a new ambient dimension ceiling.
