# Complete higher-face signatures of the original raw23 cones

All original faces of dimensions5 through9 have now been classified
exactly in each of the three frozen raw23 cones. This is an inventory,
not a new normal exclusion, capacity inequality or dimension theorem.

The [producer](../scripts/fresh_raw23_higher_face_signature_inventory_20261018.py),
[certificate](fresh_raw23_higher_face_signature_inventory_20261018.json)
and [saved checker](../scripts/check_raw23_higher_face_signature_saved_20261018.py)
use the complete50-ray lists and twelve known inequalities from the
[frozen polar packet](fresh_raw23_union_polar_completion_preflight_20261018.md).
They preserve the [four-face theorem](fresh_raw23_four_face_contact_preflight_20261018.md).

The counts per cone are identical:

|Face dimension|All faces|Entirely spacelike, positive-definite span|Other faces|
|---|---:|---:|---:|
|5|670|112|558|
|6|450|44|406|
|7|210|10|200|
|8|65|1|64|
|9|12|0|12|

The other faces contain timelike rays and have signature(d−1,1,0).
Every entirely-spacelike face found has signature(d,0,0); this is
computed rather than inferred from individual ray signs.

For dimension d, enumerate every(10−d)-subset of the twelve wall
rows. Retain a subset when its exact row rank is10−d and all frozen
rays annihilated by it span dimension d. Close it to the full common
active-wall set and deduplicate. This is complete: a d-face has an
independent(10−d)-element active basis, while any intersection with
supporting equality hyperplanes is a face. The frozen cone is pointed,
so all rays of that face are exactly the gathered original extreme
rays. No nonsimplicial ray is dropped. A greedily chosen independent
ray basis supplies the full span's Gram signature.

The standalone checker imports the earlier standalone four-face
checker's Q(√5) arithmetic, reverse-column rank elimination and
congruence inertia. It imports neither producer. It independently
reconstructs the complete active domains and all face memberships,
checks the stored ray bases and Grams and recomputes every signature.
The arithmetic design and original ray completeness are shared inputs.

Generation3159/cda7ef, separate producer replay81551/ee5047 and
independent saved checker6467/d27a72 terminated exit0/PASS. A peer
completed a full mathematical producer/report/helper read, checking
complete codimension1..5 active domains, full frozen50-ray membership,
closed-set deduplication, independent ray bases and span signatures,
including nonsimplicial and timelike-containing faces. No correction
was found. That reviewer executed no arithmetic or frozen-ray
completeness proof. The original ray-list completeness remains an
input. Source, certificate and report are final/frozen. No optimizer,
native-angle enumeration or new contact certificate is part of this packet.

The key next lead is the unique positive-definite eight-face containing
all32 spacelike rays. The separate contact-sum work tests its entire
face lattice; these signature counts alone do not exclude it.
