# Record 31, `a=1`, source orbit 2: root-4 checkpoint

The light discovery proves that matrix root class 4 has exactly 17 depth-two
classes.  Classes 1 through 8 have terminal singleton archives under the
0.60 GiB watchdog:

| class | nodes | pure 120-facet targets | literal ridge targets | peak GiB | seconds |
|---:|---:|---:|---:|---:|---:|
| 1 | 732,404 | 0 | 0 | 0.033 | 13.457 |
| 2 | 3,995,899 | 0 | 0 | 0.175 | 126.581 |
| 3 | 4,682,090 | 12 | 0 | 0.237 | 172.927 |
| 4 | 3,421,616 | 109 | 0 | 0.137 | 159.177 |
| 5 | 4,017,160 | 33 | 0 | 0.168 | 173.817 |
| 6 | 282,129 | 517 | 8 | 0.011 | 8.491 |
| 7 | 1,143,783 | 0 | 0 | 0.053 | 42.989 |
| 8 | 4,974,329 | 35 | 0 | 0.234 | 182.896 |

Class 6 contributes eight literal ridge-valid targets; the other completed
classes contribute none.  Classes 9 through 17 have not been run.  Therefore
root class 4, source orbit 2, and the Record-31 `a=1` schedule all remain
**OPEN**.  No fragment-internal counts are combined, no timeout is promoted,
and no compact or finite-volume noncompact conclusion is drawn.

The terminal archives are
`data/enumerations/h9_14_record31_a1_orbit2_matrix_root4_d2_01_01.json`
through
`data/enumerations/h9_14_record31_a1_orbit2_matrix_root4_d2_08_08.json`.
Each archive pins the generated engine
`0bf66b0bd22827f90e4034ba80539eb1cc57b712cdd59e40070b4150e074ddfc`,
compiled binary, watchdog, stderr digest, normal `child_exit`, elapsed time,
and peak RSS.
