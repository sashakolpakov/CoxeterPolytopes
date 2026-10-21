# Executed eleven-profile next-edge search

19 October 2026. Generation and separate replay **actually passed** in the
coordinator's serialized bounded wrapper. The preparation note and source
remain unchanged. No scientific job ran in this authoring lane.

The fixed28-wall search tested all11 retained short finite profiles at the
specified edge, for integral x=a+b phi with 0<x<=600 and |x'|<=2. Its1806
coefficient pairs contain1073 admissible x values, giving11803 exact
profile/x square tests. There were10 field-square solutions,20 signed
branches,8 integral primitive norm2 roots and2 surviving single-wall
extensions. This is a bounded necessary-profile search, not an all-height
root census or a chamber completion.

|Executed file|SHA256|
|---|---|
|`scripts/search_h19_next_edge_short_finite_profiles_20261019.py`|`23f496f097e35b2da90999ebebee78907122ec86fc1f62cc67cc1e0c83f5db57`|
|`results/h19_next_edge_short_finite_profiles_20261019.json`|`72a0873fbd6df501162444ef8608c572e51ef218fbac3cfe1254679a30313e0c`|
|`results/h19_next_edge_short_finite_profiles_generation_20261019.json`|`7f1b636719ea6d051e50a4638f00b2ea7a9a9e9c95d6640d97e4e532819cf3ac`|
|`results/h19_next_edge_short_finite_profiles_replay_20261019.json`|`42f6020e01031b076a7edb028057f3c6b0fd318bbac435448c3f70a637f94ad0`|

Generation/replay took0.412/0.404 seconds, with process-tree peak
RSS57,278,464/57,344,000 bytes. Both started with187.29GiB free disk;
available memory was6.12/6.14GiB. The wrapper owned the shared lock and
59-second/1GiB caps; the child used a45-second alarm. Earlier pre-launch
memory deferrals did not run the search body.

The two survivors are next_edge_profile_02 at x=12+20phi and
next_edge_profile_05 at x=6+10phi. Both have endpoint parameter
(-3+4phi)/19=(-1+2sqrt(5))/19, as encoded in the certificate in the
sqrt(5) basis. Their common finite endpoint has square
(869-389sqrt(5))/1444<0. Their products with the two full compact H8
facet normals r,w are strictly ultraparallel. The actual extension
checker will independently recheck all these identities and retention
of the prior named finite vertices; this search alone does not assert
that a completed chamber has finite volume.

The coordinator's hand identity is particularly useful. Write j for the
smaller survivor and J for the larger one, in the native E8 basis, and
put u=(0,0;0,0,1,2,3,2,1,0). Then

    J=2j-u,  u^2=6,  (j,u)=3,  j^2=J^2=2,
    (j,J)=1.

Thus both candidates are valid separately, but cannot be appended
together with the chosen inward halfspaces. A new fixed29-wall checker
selects j only. This exact identity is to be checked there; no separate
matrix computation was run while writing this note.
