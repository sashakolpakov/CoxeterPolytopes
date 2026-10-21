# Incidence-surgery source and replay review

The root reviewer inspected the exact affine-cover recursion, modular
wrapper and C++ enumerator, and the separate Fraction checker for
right/dotted patterns. Fresh saved-result replays of the final modular
program and the 151-pattern Fraction checker passed. This is a source,
mathematics and replay review; it is not a second implementation of the
3,131,392-assignment modular enumeration.

The affine-cover condition is necessary for compactness with the stated
facets retained. The recursion's dominated-column replacement preserves
minimum cover size, while the subsequent direct combination enumeration
uses original finite-pair columns. The report correctly distinguishes
complete enumeration of those minimum finite-pair supports from merely
finding the minimum when old right pairs are also eligible.

The modular argument is valid over Z[1/2,√2,√3,√5]. At each supported
prime, the chosen images satisfy the three square-root relations and 2
is invertible. A determinant that is zero in characteristic zero must
therefore reduce to zero. Nonzero pivots throughout the eleven-column
elimination certify rank eleven. The code's row updates need not clear
earlier columns, since those columns never enter later pivot steps.
All tested principal subsets avoid free old-parallel entries. The loop
enumerates four different replacement labels per changed pair and
restores each support before moving to the next. Integer intermediates
stay within the asserted bounds. Every assignment is rejected already
at prime 241.

The right/dotted checker independently reconstructs the specified
minimum-dotted masks, verifies that each ten-node witness avoids every
free pair, and checks positive rational elimination pivots. All 151
witnesses are positive definite of rank ten, so they contradict positive
index at most nine regardless of the distances. The floating-point
proposal search supplies witness indices, not trusted signs.

The [main report](fresh_h9_incidence_surgery_20261013.md) states the
remaining restrictions accurately: minimum supports, the specified label
set and the specified right/dotted rule. This establishes no obstruction
to all compact H9 chambers or to all modifications of the three sources.
The [review manifest](fresh_h9_incidence_surgery_20261013_review.json)
binds the reviewed scripts and artifacts.
