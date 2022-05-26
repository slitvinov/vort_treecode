#include <stdio.h>
#include <math.h>
#include "proj.h"
#include "comp_f.h"
#include "tree.h"

#define N_STAT 1
#define ORD_CNT 1
#define TREE_STATS 1
#define MIN_NODE_CNT 64

static const REAL CFF_0 = 0.0333, CFF_1 = 0.0063485;
static const REAL expansion_time[] = {0,    0,    0.11, 0.18, 0.30, 0.48,
                                      0.72, 1.03, 1.43, 1.93, 2.55};

#if ORD_CNT
static long ord_cnt[MAX_P], dir_cnt;
static long ord_level[MAX_P], ord_part_cnt[MAX_P], dir_level, dir_part_cnt;
#endif

#if TREE_STATS
static int v_cnt;
#endif

static void init_work(int sheet_cnt, sheet sheets[MAX_SHEET_CNT], work *scratch,
                      long *N_tot);
static vect comp_vel(vect pos, REAL del2, REAL tol, tree_cell *cell_ptr);
static vect split(vect pos, REAL del2, REAL tol, tree_cell *cell_ptr);
static vect direct_sum(vect pos, REAL del2, tree_cell *cell_ptr);
static void comp_mmnts(int p, tree_cell *cell_ptr);
static void scalar_mmnts(int p, vect diff, REAL mmnts[MAX_P][MAX_P][MAX_P]);
static vect expansion(vect diff, REAL del2, int p,
                      vect mmnts[MAX_P][MAX_P][MAX_P]);
static void comp_cff(REAL del2, vect diff, int p,
                     REAL ans[MAX_P + 1][MAX_P + 1][MAX_P + 1]);

int f(REAL del2, int sheet_cnt, sheet sheets[MAX_SHEET_CNT], work *scratch,
      REAL tol, int node_min) {
  tree_cell *root;
  int i, sheet_ind;
  long N_tot;

#if TREE_STATS
  long v_tot, v_sqr_tot;
  REAL v_avg;
  v_tot = v_sqr_tot = 0;
#endif

#if ORD_CNT
  dir_cnt = 0;
  dir_level = dir_part_cnt = 0.0;
  for (i = 0; i < MAX_P; i++) {
    ord_cnt[i] = 0;
    ord_level[i] = ord_part_cnt[i] = 0.0;
  }
#endif

#if DEBUG == RUN_PRINT
  fprintf(stderr, "inside f\n");
  fflush(stderr);
#endif

  init_work(sheet_cnt, sheets, scratch, &N_tot);

#if DEBUG == RUN_PRINT
  fprintf(stderr, "init_work complete\n");
  fflush(stderr);
#endif

#if N_STAT
  fprintf(stderr, "N_tot = %ld\n", N_tot);
#endif

  if (!(root = make_tree(0, scratch, N_tot, node_min)))
    return 1;

#if DEBUG == RUN_PRINT
  fprintf(stderr, "make_tree complete\n");
  fflush(stderr);
#endif

#if DEBUG == RUN_PRINT

  tree_dump(stderr, root);
  fprintf(stderr, "tree dump complete\n");
  fflush(stderr);
#endif

  for (sheet_ind = 0; sheet_ind < sheet_cnt; sheet_ind++) {
    int fil_ind, fil_cnt = sheets[sheet_ind].fil_cnt;
    fil *fil_ptr = sheets[sheet_ind].fils;
    for (fil_ind = 0; fil_ind < fil_cnt; fil_ptr++, fil_ind++) {
      int j, N = fil_ptr->node_cnt;
      node *n_ptr = fil_ptr->nodes;
      for (j = 0; j < N; n_ptr++, j++) {

#if TREE_STATS
        v_cnt = 0;
#endif

        n_ptr->vel = comp_vel(n_ptr->arg, del2, tol, root);

#if TREE_STATS
        v_tot += v_cnt;
        v_sqr_tot += (REAL)v_cnt * v_cnt;
#endif
      }
    }
  }

#if DEBUG == RUN_PRINT
  fprintf(stderr, "computations complete\n");
  fflush(stderr);
#endif

#if ORD_CNT
  fprintf(stderr,
          "dir_cnt = %ld dir_lev = %.4f dir_part_cnt = %.4f\n p cnt ord_lev "
          "ord_part_cnt\n",
          dir_cnt, (double)dir_level / (double)dir_cnt,
          (double)dir_part_cnt / (double)dir_cnt);
  for (i = 0; i < MAX_P; i++)
    if (!ord_cnt[i])
      fprintf(stderr, "%d %ld\n", i + 1, 0L);
    else
      fprintf(stderr, "%d %ld %.4f %.4f\n", i + 1, ord_cnt[i],
              (double)ord_level[i] / (double)ord_cnt[i],
              (double)ord_part_cnt[i] / (double)ord_cnt[i]);
  fflush(stderr);
#endif

#if TREE_STATS
  fprintf(stderr, "tree stats: depth = %d cells = %d cp_ratio = %f\n",
          depth(root), tree_cell_cnt(root), child_parent_ratio(root));
  v_avg = (double)v_tot / (double)N_tot;
  fprintf(stderr, "avg cells visits = %f std. dev. = %f\n", v_avg,
          sqrt((double)v_sqr_tot / (double)N_tot - v_avg * v_avg));
  fflush(stderr);
#endif

  return 0;
}

/**************************************************************/
/* compute tangent vectors and integration weights
 * takes advantage of the fact that the backward difference
 * for i is equal to the forward difference for i-1
 * compute the total number of particles
 */

static void init_work(int sheet_cnt, sheet sheets[MAX_SHEET_CNT], work *scratch,
                      long *N_tot) {
  int sheet_ind;

  *N_tot = 0;
  for (sheet_ind = 0; sheet_ind < sheet_cnt; sheet_ind++) {
    int fil_ind, fil_cnt = sheets[sheet_ind].fil_cnt;
    fil *fils = sheets[sheet_ind].fils;
    for (fil_ind = 0; fil_ind < fil_cnt; fil_ind++) {
      int N = fils[fil_ind].node_cnt;
      REAL kappa = fils[fil_ind].dgamma / four_pi;
      node *n_ptr = fils[fil_ind].nodes;
      vect tang_l, tang_r;
      REAL h_l, h_r, mult;
      int j, j1;

      *N_tot += N;

      /* multiply kappa by alpha integration weight */
      if (fil_cnt > 1) {
        if (fil_ind == 0)
          kappa *= (fils[1].alpha - fils[0].alpha) * 0.5;
        else if (fil_ind == fil_cnt - 1)
          kappa *= (fils[fil_cnt - 1].alpha - fils[fil_cnt - 2].alpha) * 0.5;
        else
          kappa *= (fils[fil_ind + 1].alpha - fils[fil_ind - 1].alpha) * 0.5;
      }

      /* set up so that j=0 doesn't access index -1 */
      h_r = n_ptr[0].th - (n_ptr[N - 1].th - two_pi);
      mult = 1 / h_r;
      tang_r.x1 = (n_ptr[0].arg.x1 - n_ptr[N - 1].arg.x1) * mult;
      tang_r.x2 = (n_ptr[0].arg.x2 - n_ptr[N - 1].arg.x2) * mult;
      tang_r.x3 = (n_ptr[0].arg.x3 - n_ptr[N - 1].arg.x3) * mult;

      for (j = 0; j < N; scratch++, j++) {
        j1 = (j == N - 1) ? 0 : j + 1;

        scratch->pos = n_ptr[j].arg;

        h_l = h_r;
        tang_l = tang_r;

        if (j == N - 1)
          h_r = n_ptr[j1].th + two_pi - n_ptr[j].th;
        else
          h_r = n_ptr[j1].th - n_ptr[j].th;
        mult = 1 / h_r;
        tang_r.x1 = (n_ptr[j1].arg.x1 - n_ptr[j].arg.x1) * mult;
        tang_r.x2 = (n_ptr[j1].arg.x2 - n_ptr[j].arg.x2) * mult;
        tang_r.x3 = (n_ptr[j1].arg.x3 - n_ptr[j].arg.x3) * mult;

        mult = kappa / 2;
        scratch->v_wght.x1 = mult * (h_r * tang_l.x1 + h_l * tang_r.x1);
        scratch->v_wght.x2 = mult * (h_r * tang_l.x2 + h_l * tang_r.x2);
        scratch->v_wght.x3 = mult * (h_r * tang_l.x3 + h_l * tang_r.x3);

        scratch->wght_nrm = sqrt(scratch->v_wght.x1 * scratch->v_wght.x1 +
                                 scratch->v_wght.x2 * scratch->v_wght.x2 +
                                 scratch->v_wght.x3 * scratch->v_wght.x3);
      }
    }
  }
}

static vect comp_vel(vect pos, REAL del2, REAL tol, tree_cell *cell_ptr) {
  vect diff;
  REAL dist, term, direct_time;
  int p;

#if TREE_STATS
  v_cnt++;
#endif

  if (cell_ptr->node_cnt < MIN_NODE_CNT)
    return direct_sum(pos, del2, cell_ptr);

  direct_time = CFF_0 + CFF_1 * cell_ptr->node_cnt;

  diff.x1 = pos.x1 - cell_ptr->exp_pt.x1;
  diff.x2 = pos.x2 - cell_ptr->exp_pt.x2;
  diff.x3 = pos.x3 - cell_ptr->exp_pt.x3;
  term = del2 + diff.x1 * diff.x1 + diff.x2 * diff.x2 + diff.x3 * diff.x3;
  dist = sqrt(term);

  if (dist < 1.5 * cell_ptr->r)
    return split(pos, del2, tol, cell_ptr);

#if BOUND == POTENTIAL
  term = 1 / (term * dist);
#else
  term = 1 / (term * term);
#endif

  for (p = 2; p <= MAX_P; term /= dist, p++) {
    if (expansion_time[p] > direct_time)
      break;
    if (cell_ptr->wght_nrm[p] * term < tol)
      break;
  }

  if ((p > MAX_P) || (expansion_time[p] > direct_time))
    return split(pos, del2, tol, cell_ptr);
  if (p > cell_ptr->mmnts_stat) {
    comp_mmnts(p, cell_ptr);
    cell_ptr->mmnts_stat = p;
  }

#if ORD_CNT
  ord_cnt[p - 1]++;
  ord_level[p - 1] += cell_ptr->level;
  ord_part_cnt[p - 1] += cell_ptr->node_cnt;
#endif

  return expansion(diff, del2, p, cell_ptr->mmnts);
}

static vect split(vect pos, REAL del2, REAL tol, tree_cell *cell_ptr) {
  vect vel = {0, 0, 0};
  REAL scale_tol = tol / cell_ptr->wght_nrm[0];
  int c_ind;

  if (!cell_ptr->child_cnt)
    return direct_sum(pos, del2, cell_ptr);
  for (c_ind = 0; c_ind < cell_ptr->child_cnt; c_ind++) {
    tree_cell *child_ptr = cell_ptr->child[c_ind];
    REAL child_tol = scale_tol * child_ptr->wght_nrm[0];
    vect child_vel = comp_vel(pos, del2, child_tol, child_ptr);
    vel.x1 += child_vel.x1;
    vel.x2 += child_vel.x2;
    vel.x3 += child_vel.x3;
  }
  return vel;
}

static vect direct_sum(vect pos, REAL del2, tree_cell *cell_ptr) {
  work *scratch = cell_ptr->nodes;
  vect ans = {0, 0, 0};
  long i, lim = cell_ptr->node_cnt;

#if ORD_CNT
  dir_cnt++;
  dir_level += cell_ptr->level;
  dir_part_cnt += cell_ptr->node_cnt;
#endif

  for (i = 0; i < lim; scratch++, i++) {
    vect diff;
    REAL nrm, mult;

    diff.x1 = pos.x1 - scratch->pos.x1;
    diff.x2 = pos.x2 - scratch->pos.x2;
    diff.x3 = pos.x3 - scratch->pos.x3;

    nrm =
        sqrt(del2 + diff.x1 * diff.x1 + diff.x2 * diff.x2 + diff.x3 * diff.x3);
    mult = -1 / (nrm * nrm * nrm);

    /* crss = cross(diff, scratch->v_wght);
     * ans.x1 += mult * crss.x1;
     * ans.x2 += mult * crss.x2;
     * ans.x3 += mult * crss.x3;
     */
    ans.x1 +=
        mult * (diff.x2 * scratch->v_wght.x3 - diff.x3 * scratch->v_wght.x2);
    ans.x2 +=
        mult * (diff.x3 * scratch->v_wght.x1 - diff.x1 * scratch->v_wght.x3);
    ans.x3 +=
        mult * (diff.x1 * scratch->v_wght.x2 - diff.x2 * scratch->v_wght.x1);
  }
  return ans;
}

static void comp_mmnts(int p, tree_cell *cell_ptr) {
  work *scratch = cell_ptr->nodes;
  REAL mmnts[MAX_P][MAX_P][MAX_P];
  vect exp_pt = cell_ptr->exp_pt, diff;
  long i, lim = cell_ptr->node_cnt;
  int m1, m2, m3;

  /* i == 0 */
  diff.x1 = exp_pt.x1 - scratch->pos.x1;
  diff.x2 = exp_pt.x2 - scratch->pos.x2;
  diff.x3 = exp_pt.x3 - scratch->pos.x3;
  scalar_mmnts(p, diff, mmnts);
  for (m1 = 0; m1 < p; m1++)
    for (m2 = 0; m1 + m2 < p; m2++)
      for (m3 = 0; m1 + m2 + m3 < p; m3++) {
        cell_ptr->mmnts[m1][m2][m3].x1 = scratch->v_wght.x1 * mmnts[m1][m2][m3];
        cell_ptr->mmnts[m1][m2][m3].x2 = scratch->v_wght.x2 * mmnts[m1][m2][m3];
        cell_ptr->mmnts[m1][m2][m3].x3 = scratch->v_wght.x3 * mmnts[m1][m2][m3];
      }

  for (scratch++, i = 1; i < lim; scratch++, i++) {
    diff.x1 = exp_pt.x1 - scratch->pos.x1;
    diff.x2 = exp_pt.x2 - scratch->pos.x2;
    diff.x3 = exp_pt.x3 - scratch->pos.x3;
    scalar_mmnts(p, diff, mmnts);
    for (m1 = 0; m1 < p; m1++)
      for (m2 = 0; m1 + m2 < p; m2++)
        for (m3 = 0; m1 + m2 + m3 < p; m3++) {
          cell_ptr->mmnts[m1][m2][m3].x1 +=
              scratch->v_wght.x1 * mmnts[m1][m2][m3];
          cell_ptr->mmnts[m1][m2][m3].x2 +=
              scratch->v_wght.x2 * mmnts[m1][m2][m3];
          cell_ptr->mmnts[m1][m2][m3].x3 +=
              scratch->v_wght.x3 * mmnts[m1][m2][m3];
        }
  }
}

/**************************************************************/
/* compute the weighted scalar moments of a single node in a single cell */

static void scalar_mmnts(int p, vect diff, REAL mmnts[MAX_P][MAX_P][MAX_P]) {
  int m1, m2, m3;

  /* m1 = 0, m2 = 0 */
  mmnts[0][0][0] = 1;
  for (m3 = 1; m3 < p; m3++)
    mmnts[0][0][m3] = mmnts[0][0][m3 - 1] * diff.x3;

  /* m1 = 0, m2 > 0 */
  for (m2 = 1; m2 < p; m2++) {
    mmnts[0][m2][0] = mmnts[0][m2 - 1][0] * diff.x2;
    for (m3 = 1; m2 + m3 < p; m3++)
      mmnts[0][m2][m3] = mmnts[0][m2][m3 - 1] * diff.x3;
  }

  /* m1 > 0 */
  for (m1 = 1; m1 < p; m1++) {
    /* m2 = 0 */
    mmnts[m1][0][0] = mmnts[m1 - 1][0][0] * diff.x1;
    for (m3 = 1; m1 + m3 < p; m3++)
      mmnts[m1][0][m3] = mmnts[m1][0][m3 - 1] * diff.x3;

    /* m2 > 0 */
    for (m2 = 1; m1 + m2 < p; m2++) {
      mmnts[m1][m2][0] = mmnts[m1][m2 - 1][0] * diff.x2;
      for (m3 = 1; m1 + m2 + m3 < p; m3++)
        mmnts[m1][m2][m3] = mmnts[m1][m2][m3 - 1] * diff.x3;
    }
  }
}

static vect expansion(vect diff, REAL del2, int p,
                      vect mmnts[MAX_P][MAX_P][MAX_P]) {
  vect ans = {0, 0, 0};
  REAL cff[MAX_P + 1][MAX_P + 1][MAX_P + 1];
  vect vect_a;
  int m1, m1_1, m2, m2_1, m3, m3_1;

  comp_cff(del2, diff, p, cff);
  for (m1 = 0, m1_1 = 1; m1 < p; m1_1++, m1++)
    for (m2 = 0, m2_1 = 1; m1 + m2 < p; m2_1++, m2++)
      for (m3 = 0, m3_1 = 1; m1 + m2 + m3 < p; m3_1++, m3++) {
        vect_a.x1 = m1_1 * cff[m1_1][m2][m3];
        vect_a.x2 = m2_1 * cff[m1][m2_1][m3];
        vect_a.x3 = m3_1 * cff[m1][m2][m3_1];

        /* add = cross(vect_a, mmnts[m1][m2][m3]);
         * ans.x1 += add.x1;
         * ans.x2 += add.x2;
         * ans.x3 += add.x3;
         */
        ans.x1 +=
            vect_a.x2 * mmnts[m1][m2][m3].x3 - vect_a.x3 * mmnts[m1][m2][m3].x2;
        ans.x2 +=
            vect_a.x3 * mmnts[m1][m2][m3].x1 - vect_a.x1 * mmnts[m1][m2][m3].x3;
        ans.x3 +=
            vect_a.x1 * mmnts[m1][m2][m3].x2 - vect_a.x2 * mmnts[m1][m2][m3].x1;
      }
  return ans;
}

static void comp_cff(REAL del2, vect diff, int p,
                     REAL ans[MAX_P + 1][MAX_P + 1][MAX_P + 1]) {
  REAL cff1[MAX_P + 1], cff2[MAX_P + 1], x12, x22, x32, mult;
  int i, j, k;

  for (i = 1; i <= p; i++) {
    cff1[i] = 1 - (REAL)1 / (2 * i);
    cff2[i] = 1 - (REAL)1 / i;
  }
  x12 = diff.x1 * 2;
  x22 = diff.x2 * 2;
  x32 = diff.x3 * 2;
  mult =
      -1 / (del2 + diff.x1 * diff.x1 + diff.x2 * diff.x2 + diff.x3 * diff.x3);

  /* base case */
  ans[0][0][0] = sqrt(-mult);

  /* two of the indices are zero */
  ans[0][0][1] = diff.x3 * (ans[0][0][0] * mult);
  ans[0][1][0] = diff.x2 * (ans[0][0][0] * mult);
  ans[1][0][0] = diff.x1 * (ans[0][0][0] * mult);
  for (i = 2; i <= p; i++) {
    ans[0][0][i] =
        (x32 * cff1[i] * ans[0][0][i - 1] + cff2[i] * ans[0][0][i - 2]) * mult;
    ans[0][i][0] =
        (x22 * cff1[i] * ans[0][i - 1][0] + cff2[i] * ans[0][i - 2][0]) * mult;
    ans[i][0][0] =
        (x12 * cff1[i] * ans[i - 1][0][0] + cff2[i] * ans[i - 2][0][0]) * mult;
  }

  /* one index = 0, one index = 1, other index >= 1 */
  ans[0][1][1] = 3 * diff.x3 * (ans[0][1][0] * mult);
  ans[1][0][1] = 3 * diff.x3 * (ans[1][0][0] * mult);
  ans[1][1][0] = 3 * diff.x2 * (ans[1][0][0] * mult);
  for (i = 2; i < p; i++) {
    ans[0][1][i] =
        (diff.x2 * ans[0][0][i] + x32 * ans[0][1][i - 1] + ans[0][1][i - 2]) *
        mult;
    ans[1][0][i] =
        (diff.x1 * ans[0][0][i] + x32 * ans[1][0][i - 1] + ans[1][0][i - 2]) *
        mult;
    ans[0][i][1] =
        (x22 * ans[0][i - 1][1] + ans[0][i - 2][1] + diff.x3 * ans[0][i][0]) *
        mult;
    ans[1][i][0] =
        (diff.x1 * ans[0][i][0] + x22 * ans[1][i - 1][0] + ans[1][i - 2][0]) *
        mult;
    ans[i][1][0] =
        (x12 * ans[i - 1][1][0] + ans[i - 2][1][0] + diff.x2 * ans[i][0][0]) *
        mult;
    ans[i][0][1] =
        (x12 * ans[i - 1][0][1] + ans[i - 2][0][1] + diff.x3 * ans[i][0][0]) *
        mult;
  }

  /* one index = 0, other indices >= 2 */
  for (i = 2; i <= p - 2; i++)
    for (j = 2; i + j <= p; j++) {
      ans[0][i][j] =
          (x22 * cff1[i] * ans[0][i - 1][j] + cff2[i] * ans[0][i - 2][j] +
           x32 * ans[0][i][j - 1] + ans[0][i][j - 2]) *
          mult;
      ans[i][0][j] =
          (x12 * cff1[i] * ans[i - 1][0][j] + cff2[i] * ans[i - 2][0][j] +
           x32 * ans[i][0][j - 1] + ans[i][0][j - 2]) *
          mult;
      ans[i][j][0] =
          (x12 * cff1[i] * ans[i - 1][j][0] + cff2[i] * ans[i - 2][j][0] +
           x22 * ans[i][j - 1][0] + ans[i][j - 2][0]) *
          mult;
    }

  /* two indices = 1, other index >= 1 */
  ans[1][1][1] = 5 * diff.x3 * mult * ans[1][1][0];
  for (i = 2; i <= p - 2; i++) {
    ans[1][1][i] = (diff.x1 * ans[0][1][i] + x22 * ans[1][0][i] +
                    x32 * ans[1][1][i - 1] + ans[1][1][i - 2]) *
                   mult;
    ans[1][i][1] = (diff.x1 * ans[0][i][1] + x22 * ans[1][i - 1][1] +
                    ans[1][i - 2][1] + x32 * ans[1][i][0]) *
                   mult;
    ans[i][1][1] = (x12 * ans[i - 1][1][1] + ans[i - 2][1][1] +
                    diff.x2 * ans[i][0][1] + x32 * ans[i][1][0]) *
                   mult;
  }

  /* one index = 1, other indices >= 2 */
  for (i = 2; i <= p - 3; i++)
    for (j = 2; i + j < p; j++) {
      ans[1][i][j] =
          (diff.x1 * ans[0][i][j] + x22 * ans[1][i - 1][j] + ans[1][i - 2][j] +
           x32 * ans[1][i][j - 1] + ans[1][i][j - 2]) *
          mult;
      ans[i][1][j] =
          (x12 * ans[i - 1][1][j] + ans[i - 2][1][j] + diff.x2 * ans[i][0][j] +
           x32 * ans[i][1][j - 1] + ans[i][1][j - 2]) *
          mult;
      ans[i][j][1] =
          (x12 * ans[i - 1][j][1] + ans[i - 2][j][1] + x22 * ans[i][j - 1][1] +
           ans[i][j - 2][1] + diff.x3 * ans[i][j][0]) *
          mult;
    }

  /* all indices >= 2 */
  for (i = 2; i <= p - 4; i++)
    for (j = 2; i + j <= p - 2; j++)
      for (k = 2; i + j + k <= p; k++)
        ans[i][j][k] =
            (x12 * cff1[i] * ans[i - 1][j][k] + cff2[i] * ans[i - 2][j][k] +
             x22 * ans[i][j - 1][k] + ans[i][j - 2][k] +
             x32 * ans[i][j][k - 1] + ans[i][j][k - 2]) *
            mult;
}
