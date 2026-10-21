// Exact ridge-aware fixed-positive-A4 addition census for H9/14 record 14.

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
#define H9_14_ADDITION_FIRST_ISOLATE 3
#define H9_14_RECORD14_RIDGE_KEY 1
#define H9_14_RIDGE_KIND "h9_14_record14_ridge_aware_additions"

#include "explore_h9_14_record14_ridge_aware_engine.cpp"
