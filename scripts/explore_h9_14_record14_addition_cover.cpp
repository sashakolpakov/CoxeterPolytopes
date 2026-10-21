// Exact fixed-positive-A4 addition-cover driver for H9/14 record 14.
//
// The search implementation is shared with the independently exercised
// record-18 driver.  Only the graph-derived edge clutter and the two local
// base-g arrays differ.  Usage is unchanged:
//
//   program a A4mask ... A4mask [--fast-scheduled-gate]

#define H9_14_ADDITION_EDGE_INITIALIZER                                      \
  {                                                                          \
    (1u << 0) | (1u << 12), (1u << 1) | (1u << 12),                         \
        (1u << 0) | (1u << 13), (1u << 1) | (1u << 13),                     \
        (1u << 2) | (1u << 13)                                              \
  }
#define H9_14_ADDITION_BASE_G3_INITIALIZER \
  { 4, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 2, 1 }
#define H9_14_ADDITION_BASE_G4_INITIALIZER \
  { 5, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 2, 1 }
#define H9_14_ADDITION_KIND "h9_14_record14_fixed_positive_addition_cover"

#include "explore_h9_14_record18_addition_cover.cpp"
