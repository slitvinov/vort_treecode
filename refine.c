/*  CVS:$Id: refine.c,v 1.2 2002/01/21 18:38:39 klindsay Exp $
    CVS:$Name: vort_treecode_1_0 $

    Copyright (C) 2002  Keith Lindsay (klindsay@ucar.edu)

    This file is part of vort_treecode, a program to compute the
    evolution of a three-dimensional vortex sheet using an adaptive
    treecode to evaluate the velocity field.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#if DEBUG == RUN_PRINT
#include <stdio.h>
#endif

#include <stdlib.h>
#include <math.h>
#include "proj.h"
#include "refine.h"

/**************************************************************/
/* variables local to this module */
static long M_tot;
static int realloc_flag;

/**************************************************************/
/* local prototypes */

static int refine_sheet (REAL eps_alpha2, REAL eps_th2, sheet *sheet_ptr);
static int refine_fil (REAL eps2, fil *fil_ptr);
static REAL fil_dist_sqr (fil *fil_ptr1, fil *fil_ptr2);
static int make_fil (fil *fil_ptr);
static int more_nodes (int new_cnt, fil *fil_ptr);
static int more_fils (sheet *sheet_ptr);
static vect fil_pos1 (REAL th, fil *fil_ptr);
static vect fil_pos2 (REAL th_0, int i, int N, node *n_ptr);
static vect vect_intrp (REAL r_0, REAL r[4], vect pos[4]);

/**************************************************************/

int
refine (REAL eps_alpha2, REAL eps_th2, int sheet_cnt, sheet sheets[MAX_SHEET_CNT], work **scratch)
{
   int i;
   sheet *sheet_ptr = sheets;

   M_tot = 0;
   realloc_flag = 0;

#if DEBUG == RUN_PRINT
   fprintf (stderr, "inside refine\n");
   fflush (stderr);
#endif

   for (i = 0; i < sheet_cnt; sheet_ptr++, i++)
      if (refine_sheet (eps_alpha2, eps_th2, sheet_ptr)) {
	 fprintf (stderr, "error from refine_sheet\n");
	 return 1;
      }
#if DEBUG == RUN_PRINT
   fprintf (stderr, "sheets refined\n");
   fflush (stderr);
#endif

   /* reallocate scratch if any other reallocation was done */
   if (realloc_flag) {
      void *ptr;

#if DEBUG == RUN_PRINT
      fprintf (stderr, "allocating more scratch\n");
      fflush (stderr);
#endif

      if (!(ptr = realloc (*scratch, sizeof (work) * M_tot))) {
	 fprintf (stderr, "realloc for scratch failed");
	 return 1;
      }
      *scratch = ptr;

#if DEBUG == RUN_PRINT
      fprintf (stderr, "more scratch allocated\n");
      fflush (stderr);
#endif

   }

   return 0;
}

/**************************************************************/

static int
refine_sheet (REAL eps_alpha2, REAL eps_th2, sheet *sheet_ptr)
{
   int i, j, added;
   fil *fil_ptr;

   /* refine filaments */
   fil_ptr = sheet_ptr->fils;
   for (i = 0; i < sheet_ptr->fil_cnt; fil_ptr++, i++)
      if (refine_fil (eps_th2, fil_ptr)) {
	 fprintf (stderr, "error from refine_fil");
	 return 1;
      }
#if DEBUG == RUN_PRINT
   fprintf (stderr, "fils refined\n");
   fflush (stderr);
#endif

   /* first add empty space between filaments */
   for (i = sheet_ptr->fil_cnt - 1; i > 0; i--)
      sheet_ptr->fils[2 * i] = sheet_ptr->fils[i];

   /* put new filaments in place where they are needed */
   added = 0;
   for (i = 1; i + 2 < sheet_ptr->fil_cnt; i++) {
      fil_ptr = sheet_ptr->fils + 2 * i;
      fil_ptr->ins_flag = fil_dist_sqr (fil_ptr, fil_ptr + 2) > eps_alpha2;
      if (!fil_ptr->ins_flag)
	 continue;

      if (make_fil (fil_ptr)) {
	 fprintf (stderr, "error from make_fil");
	 return 1;
      }
      added++;
      realloc_flag = 1;
      M_tot += fil_ptr[1].node_alloc;
   }

#if DEBUG == RUN_PRINT
   fprintf (stderr, "%d fils added\n", added);
   fflush (stderr);
#endif

   /* remove empty fil slots */
   sheet_ptr->fils[1] = sheet_ptr->fils[2];
   for (j = 2, i = 1; i + 2 < sheet_ptr->fil_cnt; i++) {
      if (sheet_ptr->fils[2 * i].ins_flag)
	 sheet_ptr->fils[j++] = sheet_ptr->fils[2 * i + 1];
      sheet_ptr->fils[j++] = sheet_ptr->fils[2 * i + 2];
   }
   sheet_ptr->fils[j] = sheet_ptr->fils[2 * (sheet_ptr->fil_cnt - 1)];

   sheet_ptr->fil_cnt += added;
   if (2 * sheet_ptr->fil_cnt > sheet_ptr->fil_alloc) {

#if DEBUG == RUN_PRINT
      fprintf (stderr, "allocating more fils, fil_cnt = %d, fil_alloc = %d\n", sheet_ptr->fil_cnt, sheet_ptr->fil_alloc);
      fflush (stderr);
#endif

      if (more_fils (sheet_ptr)) {
	 fprintf (stderr, "error from more_fils");
	 return 1;
      }
#if DEBUG == RUN_PRINT
      fprintf (stderr, "more fils allocated\n");
      fflush (stderr);
#endif

   }

   return 0;
}

/**************************************************************/

static int
refine_fil (REAL eps2, fil *fil_ptr)
{
   node *n_ptr = fil_ptr->nodes;
   int N, i, i1, j, add = 0;

   N = fil_ptr->node_cnt;

   for (i = 0; i < N; i++) {
      vect diff;

      i1 = (i == N - 1) ? 0 : i + 1;
      diff.x1 = n_ptr[i1].pos.x1 - n_ptr[i].pos.x1;
      diff.x2 = n_ptr[i1].pos.x2 - n_ptr[i].pos.x2;
      diff.x3 = n_ptr[i1].pos.x3 - n_ptr[i].pos.x3;

      n_ptr[i].ins_flag = diff.x1 * diff.x1 + diff.x2 * diff.x2 + diff.x3 * diff.x3 > eps2;

      if (n_ptr[i].ins_flag) {
	 REAL th;
	 if (i == N - 1)
	    th = (n_ptr[i].th + two_pi) / 2;
	 else
	    th = (n_ptr[i].th + n_ptr[i + 1].th) / 2;
	 n_ptr[i].new = fil_pos2 (th, i, N, n_ptr);
	 add++;
      }
   }

   /* allocate more space if necessary */
   if (N + add > fil_ptr->node_alloc) {
      if (more_nodes (N + add, fil_ptr))
	 return 1;
      realloc_flag = 1;
      n_ptr = fil_ptr->nodes;
   }

   /* insert new positions in */
   for (i = N - 1, j = N - 1 + add; i < j; --i, --j) {
      if (n_ptr[i].ins_flag) {
	 if (i == N - 1)
	    n_ptr[j].th = (n_ptr[i].th + two_pi) / 2;
	 else
	    n_ptr[j].th = (n_ptr[i].th + n_ptr[i + 1].th) / 2;
	 n_ptr[j].pos = n_ptr[i].new;
	 --j;
      }
      n_ptr[j].th = n_ptr[i].th;
      n_ptr[j].pos = n_ptr[i].pos;
   }

   fil_ptr->node_cnt += add;
   M_tot += fil_ptr->node_alloc;

   return 0;
}

/**************************************************************/
/* Approximate the square of the distance from the first
 * filament to the second one. Works by maximizing distance
 * from each node on first filament to corresponding point
 * on second one.
 */

static REAL
fil_dist_sqr (fil *fil_ptr1, fil *fil_ptr2)
{
   REAL ans = 0;
   node *n_ptr1 = fil_ptr1->nodes;
   int N1 = fil_ptr1->node_cnt, ind;

   for (ind = 0; ind < N1; ind++) {
      vect pos2, diff;
      REAL dist_sqr;

      pos2 = fil_pos1 (n_ptr1[ind].th, fil_ptr2);
      diff.x1 = n_ptr1[ind].pos.x1 - pos2.x1;
      diff.x2 = n_ptr1[ind].pos.x2 - pos2.x2;
      diff.x3 = n_ptr1[ind].pos.x3 - pos2.x3;
      dist_sqr = diff.x1 * diff.x1 + diff.x2 * diff.x2 + diff.x3 * diff.x3;
      ans = MAX (ans, dist_sqr);
   }

   return ans;
}

/**************************************************************/

static int
make_fil (fil *fil_ptr)
{
   void *ptr;
   REAL new_alpha, alpha[4];
   int N, i;

#if DEBUG == RUN_PRINT
   fprintf (stderr, "adding new fil\n");
   fflush (stderr);
#endif

   fil_ptr[1].node_cnt = N = fil_ptr[2].node_cnt;
   fil_ptr[1].node_alloc = fil_ptr[2].node_alloc;
   if (!(ptr = malloc (sizeof (node) * fil_ptr[1].node_alloc))) {
      fprintf (stderr, "malloc for nodes for new fil failed");
      return 1;
   }
   fil_ptr[1].nodes = ptr;

   /* compute alpha */
   alpha[0] = fil_ptr[-2].alpha;
   alpha[1] = fil_ptr[0].alpha;
   alpha[2] = fil_ptr[2].alpha;
   alpha[3] = fil_ptr[4].alpha;
   fil_ptr[1].alpha = new_alpha = (fil_ptr[0].alpha + fil_ptr[2].alpha) / 2;

   /* compute dgamma analytically */
   fil_ptr[1].dgamma = -sin (new_alpha);

#if 0
   /* interpolate to find dgamma */
   dgamma[0] = fil_ptr[-2].dgamma;
   dgamma[1] = fil_ptr[0].dgamma;
   dgamma[2] = fil_ptr[2].dgamma;
   dgamma[3] = fil_ptr[4].dgamma;

   for (step = 1; step < 4; step++)
      for (j = 3; j >= step; --j)
	 dgamma[j] = (dgamma[j] - dgamma[j - 1]) / (alpha[j] - alpha[j - step]);

   new_dgamma = (new_alpha - alpha[2]) * dgamma[3] + dgamma[2];
   new_dgamma = (new_alpha - alpha[1]) * new_dgamma + dgamma[1];
   fil_ptr[1].dgamma = (new_alpha - alpha[0]) * new_dgamma + dgamma[0];
#endif

   /* interpolate to find the point positions, using next
    * outer filaments theta placements
    */
   for (i = 0; i < N; i++) {
      vect pos[4];
      REAL th = fil_ptr[2].nodes[i].th;

      pos[0] = fil_pos1 (th, fil_ptr - 2);
      pos[1] = fil_pos1 (th, fil_ptr);
      pos[2] = fil_ptr[2].nodes[i].pos;
      pos[3] = fil_pos1 (th, fil_ptr + 4);

      fil_ptr[1].nodes[i].th = th;
      fil_ptr[1].nodes[i].pos = vect_intrp (new_alpha, alpha, pos);
   }

   return 0;
}

/**************************************************************/

static int
more_nodes (int new_cnt, fil *fil_ptr)
{
   void *ptr;

   fil_ptr->node_alloc = (new_cnt / ARRAY_BLOCK + 1) * ARRAY_BLOCK;
   if (!(ptr = realloc (fil_ptr->nodes, sizeof (node) * fil_ptr->node_alloc))) {
      fprintf (stderr, "realloc for nodes failed");
      return 1;
   }
   fil_ptr->nodes = ptr;

   return 0;
}

/**************************************************************/

static int
more_fils (sheet *sheet_ptr)
{
   void *ptr;

   sheet_ptr->fil_alloc = (2 * sheet_ptr->fil_cnt / ARRAY_BLOCK + 1) * ARRAY_BLOCK;
   if (!(ptr = realloc (sheet_ptr->fils, sizeof (fil) * sheet_ptr->fil_alloc))) {
      fprintf (stderr, "realloc for fils failed");
      return 1;
   }
   sheet_ptr->fils = ptr;

   return 0;
}

/**************************************************************/
/* wrapper for fil_pos2 for cubic interpolation along a filament */
/* this routine finds the index of the interval containing the point */

static vect
fil_pos1 (REAL th, fil *fil_ptr)
{
   node *n_ptr = fil_ptr->nodes;
   int N = fil_ptr->node_cnt, left, right, mid;

   if (th >= n_ptr[N - 1].th) {
      left = N - 1;
   } else {
      left = 0;
      right = N - 1;
      while (left + 1 < right) {
	 mid = left + (right - left) / 2;
	 if (n_ptr[mid].th < th)
	    left = mid;
	 else if (n_ptr[mid].th > th)
	    right = mid;
	 else
	    return n_ptr[mid].pos;
      }
   }

   if (n_ptr[left].th == th)
      return n_ptr[left].pos;
   else
      return fil_pos2 (th, left, N, n_ptr);
}

/**************************************************************/
/* cubic interpolation w.r.t. theta along a filament */

static vect
fil_pos2 (REAL th_0, int i, int N, node *n_ptr)
{
   vect pos[4];
   REAL th[4];
   int j;

   for (j = 0; j < 4; j++) {
      if (i - 1 + j < 0) {
	 th[j] = n_ptr[i - 1 + j + N].th - two_pi;
	 pos[j] = n_ptr[i - 1 + j + N].pos;
      } else if (i - 1 + j >= N) {
	 th[j] = n_ptr[i - 1 + j - N].th + two_pi;
	 pos[j] = n_ptr[i - 1 + j - N].pos;
      } else {
	 th[j] = n_ptr[i - 1 + j].th;
	 pos[j] = n_ptr[i - 1 + j].pos;
      }
   }

   return vect_intrp (th_0, th, pos);
}

/**************************************************************/

static vect
vect_intrp (REAL r_0, REAL r[4], vect pos[4])
{
   vect ans;
   int step, j;

   for (step = 1; step < 4; step++) {
      for (j = 3; j >= step; --j) {
	 REAL mult = 1 / (r[j] - r[j - step]);
	 pos[j].x1 = (pos[j].x1 - pos[j - 1].x1) * mult;
	 pos[j].x2 = (pos[j].x2 - pos[j - 1].x2) * mult;
	 pos[j].x3 = (pos[j].x3 - pos[j - 1].x3) * mult;
      }
   }

   ans = pos[3];
   for (j = 2; j >= 0; --j) {
      ans.x1 = (r_0 - r[j]) * ans.x1 + pos[j].x1;
      ans.x2 = (r_0 - r[j]) * ans.x2 + pos[j].x2;
      ans.x3 = (r_0 - r[j]) * ans.x3 + pos[j].x3;
   }

   return ans;
}
