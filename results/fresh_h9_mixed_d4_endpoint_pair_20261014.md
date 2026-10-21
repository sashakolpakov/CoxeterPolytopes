# A coupled endpoint pair for the mixed D4 chamber

Two exact caps extend the [twelve-wall mixed D4 chamber](fresh_h9_mixed_d4_coupling_20261014.md) to **fourteen actual facets**. They are mutually right, close two old open edges, and add three finite vertices connected by a four-edge path between old finite vertices. Both caps cut the previous strict escape. A new strict spacelike escape remains: **the resulting chamber still has infinite volume**.

The [producer](../scripts/fresh_h9_mixed_d4_endpoint_pair_20261014.py) and [exact data](fresh_h9_mixed_d4_endpoint_pair_20261014.json) are frozen. A separate [independent review](fresh_h9_mixed_d4_endpoint_pair_20261014_review.md) reconstructs all 2,002 nine-wall supports, rather than repeating the producer's cone-clipping implementation.

## Exact cap recipes

Both cap normals have square 2 and literal coordinates in `Q(sqrt5)`; the retained ambient geometry still uses its original quartic field. They are general endpoint caps, not polar normals of old projective extreme rays. Their prescribed contacts are:

| New wall | Old open edge | Rank-eight support | Nonright support contacts | Additional anchor |
|---:|---:|---|---|---|
| 12 | 6 | 0,2,3,4,5,7,9,10 | walls 2,9 at π/3 | wall 1 at π/5 |
| 13 | 78 | 0,1,3,4,5,7,9,10 | walls 5,7 at π/3 | wall 2 at π/5 |

All other contacts within each displayed support are right. The mutual product of walls 12 and 13 is exactly zero. The 25 new pair relations comprise 13 right, four π/3, two π/5 and six dotted pairs. No finite-angle identification relies on numerical fitting.

Both new normals are strict at every old finite vertex. An explicit future timelike common center and fourteen actual facet witnesses are saved. The old source walls are all retained, so the inherited conditional non-quasi-arithmetic obstruction remains available if a finite-volume completion is eventually found.

## New finite connections and the remaining gap

The complete projective cone has **15 finite and 91 spacelike extreme rays, with no ideal ray**. All fifteen finite vertices are simple. The three new vertex diagrams are `A8+A1`, `D7+I2(5)` and `E8+A1`; the last vertex lies on both caps.

Using the final cone's ray IDs, the new finite path is

`14 — 18 — 74 — 67 — 66`.

The endpoints 14 and 66 are old vertices. The intermediate three vertices are new; 74 is the joint-cap corner. This adds four finite edges and one cycle to the connected finite-edge graph. The two outer edges of the path close the original open directions 6 and 78.

The complete queue at the fifteen finite vertices contains 64 directed finite closures, representing **32 distinct finite edges**, and **71 open directions**. These counts cover edges incident to finite vertices, not vertex-free edges or every remaining boundary component.

The saved new escape has rational coordinates, positive Lorentz square and strict negative pairing with every one of the fourteen walls. Its pairing with the common future center is negative. It therefore certifies an open ideal-boundary region and infinite volume. Cutting the previous witness and creating a finite-neighbor cycle do not establish global closure.

Generation and saved replay pass for both the producer and independent review:

```sh
python scripts/fresh_h9_mixed_d4_endpoint_pair_20261014.py --verify
python scripts/fresh_h9_mixed_d4_endpoint_pair_20261014_review.py --verify
```

All prior chamber and frontier artifacts remain unchanged. This is a new compatible contact mechanism and exact intermediate chamber, not a compact H9 example or a new finite-volume dimension.
