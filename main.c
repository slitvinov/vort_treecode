#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "proj.h"
#include "comp_f.h"
#include "refine.h"
#include "io.h"
#include "tree.h"
#include "disk_stats.h"

#define DT_STATS 1
#define DISK_STATS 1

const REAL pi = 3.1415926535897932384626433832795;
const REAL two_pi = 6.2831853071795864769252867665590;
const REAL four_pi = 12.566370614359172953850573533118;

static int get_args(char **argv, FILE **infile, REAL *dtmin, REAL *eps_alpha2,
                    REAL *eps_th2, REAL *stop, REAL *tol, int *node_min,
                    FILE **outfile);
static int rk_step(REAL dt, REAL del2, int sheet_cnt,
                   sheet sheets[MAX_SHEET_CNT], work *scratch, REAL tol,
                   int node_min, REAL *vel);
static vect v_add(vect a, REAL beta, vect b);
static void free_sheets(int sheet_cnt, sheet sheets[MAX_SHEET_CNT]);

int main(int argc, char **argv) {
  FILE *infile, *outfile;
  sheet sheets[MAX_SHEET_CNT];
  work *scratch;
  REAL dtmin, eps_alpha2, eps_th2, stop, tol, del2, t, dt0, vel0, vel;
  int node_min, sheet_cnt, status;

  if (argc != 9) {
    fprintf(stderr,
            "usage should be : prog <infile> <dtmin> <eps_alpha> <eps_th> "
            "<stop_time> <tol> <node_min> <outfile>\n");
    return EXIT_FAILURE;
  }
  if (get_args(argv, &infile, &dtmin, &eps_alpha2, &eps_th2, &stop, &tol,
               &node_min, &outfile))
    return EXIT_FAILURE;
  if (readin(infile, &del2, &t, &dt0, &vel0, &vel, &sheet_cnt, sheets,
             &scratch))
    return EXIT_FAILURE;

#if DEBUG == STARTUP
  status = dump(del2, t, dt0, vel0, vel, sheet_cnt, sheets, outfile);
  free_tree();
  free(scratch);
  free_sheets(sheet_cnt, sheets);
  return status ? EXIT_FAILURE : EXIT_SUCCESS;
#endif

  do {
    REAL dt, factor;

    if (t == 0) {
      dt = dt0;
    } else {
      if ((factor = vel / vel0) < 1)
        factor = 1;
      if ((dt = dt0 / factor) < dtmin)
        dt = dtmin;
    }
    if (stop - t < 1.1 * dt)
      dt = stop - t;

#if DT_STATS
    fprintf(stderr, "t = %11.5e, dt = %11.5e\n", t, dt);
    fflush(stderr);
#endif

#if DISK_STATS
    disk_stats(sheets);
    fflush(stdout);
#endif

    if (rk_step(dt, del2, sheet_cnt, sheets, scratch, tol, node_min, &vel))
      break;
    if (t == 0)
      vel0 = vel;
    if (refine(eps_alpha2, eps_th2, sheet_cnt, sheets, &scratch))
      break;

    t += dt;

  } while (t < stop);

  status = dump(del2, t, dt0, vel0, vel, sheet_cnt, sheets, outfile);
  free_tree();
  free(scratch);
  free_sheets(sheet_cnt, sheets);
  return status ? EXIT_FAILURE : EXIT_SUCCESS;
}

/**************************************************************/
/* order of arguments is :
 *    infile, dtmin, eps_alpha, eps_th, stop_time, tol, node_min, outfile
 */

static int get_args(char **argv, FILE **infile, REAL *dtmin, REAL *eps_alpha2,
                    REAL *eps_th2, REAL *stop, REAL *tol, int *node_min,
                    FILE **outfile) {
  REAL eps;

  if (!(*infile = fopen(*++argv, "rb"))) {
    fprintf(stderr, "fopen for infile on '%s' failed", *argv);
    return 1;
  }

  *dtmin = atof(*++argv);
  eps = atof(*++argv);
  *eps_alpha2 = eps * eps;
  eps = atof(*++argv);
  *eps_th2 = eps * eps;
  *stop = atof(*++argv);
  *tol = atof(*++argv);
  *node_min = atoi(*++argv);

  if (!(*outfile = fopen(*++argv, "wb"))) {
    fprintf(stderr, "fopen for outfile on '%s' failed", *argv);
    return 1;
  }

  return 0;
}

static int rk_step(REAL dt, REAL del2, int sheet_cnt,
                   sheet sheets[MAX_SHEET_CNT], work *scratch, REAL tol,
                   int node_min, REAL *vel) {
  REAL dt_2, dt_3, dt_6, vel2 = 0;
  int sheet_ind;

  dt_2 = dt / 2;
  dt_3 = dt / 3;
  dt_6 = dt / 6;

  for (sheet_ind = 0; sheet_ind < sheet_cnt; sheet_ind++) {
    int i, fil_cnt = sheets[sheet_ind].fil_cnt;
    fil *fil_ptr = sheets[sheet_ind].fils;
    for (i = 0; i < fil_cnt; fil_ptr++, i++) {
      int j, N = fil_ptr->node_cnt;
      node *n_ptr = fil_ptr->nodes;
      for (j = 0; j < N; n_ptr++, j++) {
        n_ptr->arg = n_ptr->pos;
      }
    }
  }

  /* stage 1 & compute velocity */
  if (f(del2, sheet_cnt, sheets, scratch, tol, node_min))
    return 1;
  for (sheet_ind = 0; sheet_ind < sheet_cnt; sheet_ind++) {
    int i, fil_cnt = sheets[sheet_ind].fil_cnt;
    fil *fil_ptr = sheets[sheet_ind].fils;
    for (i = 0; i < fil_cnt; fil_ptr++, i++) {
      int j, N = fil_ptr->node_cnt;
      node *n_ptr = fil_ptr->nodes;
      for (j = 0; j < N; n_ptr++, j++) {

        REAL tmp = n_ptr->vel.x1 * n_ptr->vel.x1 +
                   n_ptr->vel.x2 * n_ptr->vel.x2 +
                   n_ptr->vel.x3 * n_ptr->vel.x3;
        vel2 = MAX(vel2, tmp);

        n_ptr->new = v_add(n_ptr->pos, dt_6, n_ptr->vel);
        n_ptr->arg = v_add(n_ptr->pos, dt_2, n_ptr->vel);
      }
    }
  }
  *vel = sqrt(vel2);

  /* stage 2 */
  if (f(del2, sheet_cnt, sheets, scratch, tol, node_min))
    return 1;
  for (sheet_ind = 0; sheet_ind < sheet_cnt; sheet_ind++) {
    int i, fil_cnt = sheets[sheet_ind].fil_cnt;
    fil *fil_ptr = sheets[sheet_ind].fils;
    for (i = 0; i < fil_cnt; fil_ptr++, i++) {
      int j, N = fil_ptr->node_cnt;
      node *n_ptr = fil_ptr->nodes;
      for (j = 0; j < N; n_ptr++, j++) {
        n_ptr->new = v_add(n_ptr->new, dt_3, n_ptr->vel);
        n_ptr->arg = v_add(n_ptr->pos, dt_2, n_ptr->vel);
      }
    }
  }

  /* stage 3 */
  if (f(del2, sheet_cnt, sheets, scratch, tol, node_min))
    return 1;
  for (sheet_ind = 0; sheet_ind < sheet_cnt; sheet_ind++) {
    int i, fil_cnt = sheets[sheet_ind].fil_cnt;
    fil *fil_ptr = sheets[sheet_ind].fils;
    for (i = 0; i < fil_cnt; fil_ptr++, i++) {
      int j, N = fil_ptr->node_cnt;
      node *n_ptr = fil_ptr->nodes;
      for (j = 0; j < N; n_ptr++, j++) {
        n_ptr->new = v_add(n_ptr->new, dt_3, n_ptr->vel);
        n_ptr->arg = v_add(n_ptr->pos, dt, n_ptr->vel);
      }
    }
  }

  /* stage 4 */
  if (f(del2, sheet_cnt, sheets, scratch, tol, node_min))
    return 1;
  for (sheet_ind = 0; sheet_ind < sheet_cnt; sheet_ind++) {
    int i, fil_cnt = sheets[sheet_ind].fil_cnt;
    fil *fil_ptr = sheets[sheet_ind].fils;
    for (i = 0; i < fil_cnt; fil_ptr++, i++) {
      int j, N = fil_ptr->node_cnt;
      node *n_ptr = fil_ptr->nodes;
      for (j = 0; j < N; n_ptr++, j++) {
        n_ptr->pos = v_add(n_ptr->new, dt_6, n_ptr->vel);
      }
    }
  }

  return 0;
}

static vect v_add(vect a, REAL beta, vect b) {
  vect ans;

  ans.x1 = a.x1 + beta * b.x1;
  ans.x2 = a.x2 + beta * b.x2;
  ans.x3 = a.x3 + beta * b.x3;

  return ans;
}

static void free_sheets(int sheet_cnt, sheet sheets[MAX_SHEET_CNT]) {
  int sheet_ind;

  for (sheet_ind = 0; sheet_ind < sheet_cnt; sheet_ind++) {
    int fil_ind;
    for (fil_ind = 0; fil_ind < sheets[sheet_ind].fil_cnt; fil_ind++) {
      free(sheets[sheet_ind].fils[fil_ind].nodes);
    }
    free(sheets[sheet_ind].fils);
  }
}
