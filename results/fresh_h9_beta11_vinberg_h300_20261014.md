# Full reflection enumeration through height 300

The complete height-300 prefix has **35 actual walls**: 31 of norm 2 and four
of norm 2(3+phi). It exhausts 80 positive height levels using 1,082,417 exact
enumeration nodes. All 46 levels and 21 walls of the frozen height-100 prefix
are reproduced literally. Generation and saved replay pass.

The lattice is E8 perpendicular to the binary form [[-2phi,-1],[-1,2]], with
phi²=phi+1. It is the same lattice as the repository's older H11 investigation,
after an explicit change of basis. The earlier search restricted itself to
norm-two reflections. Here the complete primitive norm alphabet is used and
the two classes are interleaved by geometric height
`(1+4phi)x²/(2q)`. This distinction changes the chamber from its first positive
wall onward.

Every pair of walls has a verified Coxeter relation: 278 right angles, 66
angles pi/3, 37 angles pi/5 and 214 ultraparallel pairs. A strict timelike
center and a relative interior point for each of the 35 facets are stored.
These checks establish an actual Coxeter chamber, but not finite volume or
termination of the reflection algorithm. A separate global geometry check is
needed.

The [wrapper](../scripts/fresh_h9_beta11_vinberg_h300_20261014.py) preserves the
[height-100 implementation](../scripts/fresh_h9_beta11_vinberg_20261014.py)
and writes a fresh [certificate](fresh_h9_beta11_vinberg_h300_20261014.json).
The [independent review](fresh_h9_beta11_vinberg_h300_20261014_review.md)
uses different arithmetic and reversed coordinate enumeration, checking all
80 levels rather than merely the accepted roots.

```sh
python scripts/fresh_h9_beta11_vinberg_h300_20261014.py --verify
```
