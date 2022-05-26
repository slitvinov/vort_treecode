/**************************************************************/
/* REAL is the floating point type to be used */
#undef REAL
#define REAL double

/* what level run is to be done */
#define NONE 0
#define STARTUP 1
#define RUN_PRINT 2
#define DEBUG NONE

#define ABS(a) ((a) < 0 ? (-(a)) : (a))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

/* parameters governing dynamic arrays */
#define ARRAY_BLOCK 128

#define MAX_SHEET_CNT 4

typedef struct {
  REAL x1, x2, x3;
} vect;

typedef struct {
  vect pos, v_wght;
  REAL wght_nrm;
} work;

typedef struct {
  REAL th;
  vect pos, arg, new, vel;
  int ins_flag;
} node;

typedef struct {
  REAL alpha, dgamma;
  node *nodes;
  int node_cnt, node_alloc, ins_flag;
} fil;

typedef struct {
  fil *fils;
  int fil_cnt, fil_alloc;
} sheet;

extern const REAL pi, two_pi, four_pi;
