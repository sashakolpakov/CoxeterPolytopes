# Independent root review: a positive carrier-excess fraction in dimension12

18 October 2026. **PASS.** The universal proof and complete exact interval
checker in [the carrier-fraction packet](right_angled_mixed_carrier_excess_fraction_research_20261018.md)
were read independently of their author. A fresh root replay also passes.
For a finite-volume right-angled twelve-polytope with f finite vertices,

    E=sum_non-P3(facets(Q)+cusps(Q)-10) >= f/5390.

At least f/385 finite vertices meet a three-face with facet-plus-cusp sum
at least11. These are necessary incidence bounds, not an exclusion of12.

## Actual moves and the count conversion

The graph contains only actual finite edges of cosh length sqrt2 and
actual P3 opposite-vertex moves, weighted1 and1/20. Both relations are
undirected; their distances differ. At a finite vertex the right-angle
normals are orthonormal. For one global future unit timelike u, the
positive function g(q)=1/(-u.q) gives the stated reciprocal energy ratios.
The sum of row-minus-degree is nonnegative by x+1/x>=2 on each edge.
No componentwise invariant or unproved formal neighbor is required.

At a vertex all of whose three-faces are small, the inherited P3/Q8/R7
classification forces all finite edges to have cosh length sqrt2 and
every triple of ideal directions to be an actual P3. The earlier
independently reviewed R7 argument supplies the length. Hence the whole
small-vertex row, rather than merely selected formal contributions, has
the displayed finite-group and triple-group form.

Replacing signed frame coefficients by their absolute values and then
increasing their norm to one enlarges the row; all denominators remain
positive. The timelike frame initially has squared norm strictly below
one, which causes no boundary problem. The certificate proves a row
deficit at least1/8 at each such vertex.

At every other vertex the bounded move lengths give row excess below48.
The global reciprocal-pair inequality therefore gives b>=f/385. Each
exceptional three-face has at most14(nu-10) finite vertices: nu=11 uses
facets+2cusps>=12, and nu>=12 uses p<=2nu-4. Counting incidences bounds
b<=14E. This is an ordinary unweighted vertex count, with no remaining
energy weights or normalization.

## Calculus, intervals and signs

The mean/variance reciprocal bound follows from an exact identity for
positive denominators. Both derivative numerators in the report were
independently differentiated. Each is affine and decreasing in its group
mean, so clipping its unique critical point gives the maximum on the
whole group-mean interval. All enlarged denominators stay positive.

The comparison with the rational interval function must precede its
maximization on the enlarged domain. At every actual finite-group point,
y^2<=x/m<=h/m, so its inflated variance numerator is nonnegative.
Lowering sqrt2 and the denominator bounds is therefore valid there.
The ideal group has the same property, ky^2<=actual_norm<=inflated_norm.
Only after this pointwise comparison is the domain enlarged. Negative
variance values at added, nonphysical points cannot invalidate a maximum
over a superset. Thus this subtle sign issue does not create a gap.

The178 closed rational boxes cover x in[0,1] for every k=0,...,12,
with the degenerate norm splits at k=0,12 handled separately. Integer
square-root enclosures and the strict rational lower bound for sqrt2
are checked by squaring. All resulting gaps exceed0.125535, hence1/8.
This is interval verification of every real direction, not sampling.

## Fresh execution and scope

Root chunkb15b67 returned:

```
PASS 178 exact rational interval boxes: n12 carrier excess E>=f/5390; no optimizer or dimension exclusion
```

The run held the shared exclusive lock and used the signal-safe bounded
process wrapper,59-second limit and normal resource gates. It took0.023s;
the final gate reported196.88GiB disk free and8.61GiB available memory.
No optimizer, Lean process or temporary proof output was created.

The three-dimensional classification/attachment and R7 geometry remain
explicit inherited internal inputs. This review does not reprove all
earlier published premises. Full Lean verification and external
refereeing are pending. The constant is currently too small to close
the known dimension12 inequality gap. No dimension11 fraction is claimed.
Frozen original proof, checker and certificate bytes are preserved.
