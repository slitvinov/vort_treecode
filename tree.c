/*  CVS:$Id: tree.c,v 1.2 2002/01/21 18:38:39 klindsay Exp $
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

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include "proj.h"
#include "tree.h"

#define C_OF_MASS_EXP_PT 0
#define CELL_BLOCK_SIZE 16
#define CELL_SPLIT_FACTOR 0.707

/**************************************************************/
/* local structure definitions */

typedef struct cell_block {
   tree_cell cells[CELL_BLOCK_SIZE];
   struct cell_block *next_block;
   int cells_used;
} cell_block;

/**************************************************************/
/* local variables */

static cell_block *head_block, *curr_block;

/**************************************************************/
/* local prototypes */

static tree_cell *split_node (int level, tree_cell *t_ptr, vect cut_pos, int node_min);
static REAL coord (int dir, vect pos);
static void swap (work *nodes, long i, long j);
static int init_tree (void);
static tree_cell *cell_alloc (void);
static cell_block *new_block (void);

/**************************************************************/
/* recursively construct the spatial-tree
 * The routine is given a cell which has has been initialized
 * and the node_list has been constructed.
 * It computes the children the cell will have, initializes them
 * and calls itself on each child if the child is not a leaf and
 * the child has enough nodes
 */

tree_cell *
make_tree (int level, work *nodes, long node_cnt, int node_min)
{
   work *w_ptr;
   tree_cell *t_ptr;
   vect pos_min, pos_max, exp_pt, width_2;
   REAL max_r, wght_sum;
   long i;
   int p;

   if (level == 0) {
      if (init_tree ())
	 return NULL;
   }

   if (!(t_ptr = cell_alloc ()))
      return NULL;

   for (i = 0; i < MAX_CHILD_CNT; i++) {
      t_ptr->child[i] = NULL;
   }

   /* compute bounding box and exp_pt */
   w_ptr = nodes;
   pos_min = pos_max = w_ptr->pos;
   wght_sum = 0;
   exp_pt.x1 = exp_pt.x2 = exp_pt.x3 = 0;
   for (i = 0; i < node_cnt; w_ptr++, i++) {
      if (w_ptr->pos.x1 < pos_min.x1)
	 pos_min.x1 = w_ptr->pos.x1;
      if (w_ptr->pos.x2 < pos_min.x2)
	 pos_min.x2 = w_ptr->pos.x2;
      if (w_ptr->pos.x3 < pos_min.x3)
	 pos_min.x3 = w_ptr->pos.x3;

      if (w_ptr->pos.x1 > pos_max.x1)
	 pos_max.x1 = w_ptr->pos.x1;
      if (w_ptr->pos.x3 > pos_max.x3)
	 pos_max.x3 = w_ptr->pos.x3;
      if (w_ptr->pos.x2 > pos_max.x2)
	 pos_max.x2 = w_ptr->pos.x2;

#if C_OF_MASS_EXP_PT
      wght_sum += w_ptr->wght_nrm;
      exp_pt.x1 += w_ptr->wght_nrm * w_ptr->pos.x1;
      exp_pt.x2 += w_ptr->wght_nrm * w_ptr->pos.x2;
      exp_pt.x3 += w_ptr->wght_nrm * w_ptr->pos.x3;
   }
   exp_pt.x1 /= wght_sum;
   exp_pt.x2 /= wght_sum;
   exp_pt.x3 /= wght_sum;
#else
   }
   exp_pt.x1 = (pos_max.x1 + pos_min.x1) / 2;
   exp_pt.x2 = (pos_max.x2 + pos_min.x2) / 2;
   exp_pt.x3 = (pos_max.x3 + pos_min.x3) / 2;
#endif

   width_2.x1 = (pos_max.x1 - pos_min.x1) / 2;
   width_2.x2 = (pos_max.x2 - pos_min.x2) / 2;
   width_2.x3 = (pos_max.x3 - pos_min.x3) / 2;

   /* compute wght_nrm's and max_r */
   for (p = 0; p <= MAX_P; p++) {
      t_ptr->wght_nrm[p] = 0;
   }
   max_r = 0;
   for (w_ptr = nodes, i = 0; i < node_cnt; w_ptr++, i++) {
      vect diff;
      REAL r, term;

      diff.x1 = w_ptr->pos.x1 - exp_pt.x1;
      diff.x2 = w_ptr->pos.x2 - exp_pt.x2;
      diff.x3 = w_ptr->pos.x3 - exp_pt.x3;
      r = sqrt (diff.x1 * diff.x1 + diff.x2 * diff.x2 + diff.x3 * diff.x3);
      max_r = MAX (max_r, r);
      term = w_ptr->wght_nrm;
      for (p = 0; p <= MAX_P; term *= r, p++) {

#if BOUND == POTENTIAL
	 t_ptr->wght_nrm[p] += term;
#else
	 t_ptr->wght_nrm[p] += (p + 1) * (p + 1) * term;
#endif

      }
   }

   t_ptr->exp_pt = exp_pt;
   t_ptr->width_2 = width_2;
   t_ptr->r = max_r;
   t_ptr->level = level;
   t_ptr->mmnts_stat = 0;
   t_ptr->nodes = nodes;
   t_ptr->node_cnt = node_cnt;

   if (node_cnt > node_min) {
      vect cut_pos;
      cut_pos.x1 = (pos_max.x1 + pos_min.x1) / 2;
      cut_pos.x2 = (pos_max.x2 + pos_min.x2) / 2;
      cut_pos.x3 = (pos_max.x3 + pos_min.x3) / 2;
      return split_node (level, t_ptr, cut_pos, node_min);
   } else {
      t_ptr->child_cnt = 0;
      return t_ptr;
   }
}

/**************************************************************/

static tree_cell *
split_node (int level, tree_cell *t_ptr, vect cut_pos, int node_min)
{
   work *nodes = t_ptr->nodes;
   REAL max_width_2;
   int i, child_cnt, true_child_cnt, dir;
   long node_cnt[MAX_CHILD_CNT], c_node_ind[MAX_CHILD_CNT];

   max_width_2 = MAX (t_ptr->width_2.x1, t_ptr->width_2.x2);
   max_width_2 = MAX (max_width_2, t_ptr->width_2.x3);

   child_cnt = 1;
   node_cnt[0] = t_ptr->node_cnt;
   c_node_ind[0] = 0;

   /* partition nodes in direction dir if it is wide enough */
   for (dir = 1; dir <= 3; dir++) {
      REAL cut;

      if (coord (dir, t_ptr->width_2) < CELL_SPLIT_FACTOR * max_width_2)
	 continue;

      for (i = child_cnt - 1; i >= 0; --i) {
	 node_cnt[2 * i] = node_cnt[i];
	 c_node_ind[2 * i] = c_node_ind[i];
      }

      cut = coord (dir, cut_pos);
      for (i = 0; i < child_cnt; i++) {
	 long l, l_lim, r, r_lim;

	 /* nothing to partition */
	 if (node_cnt[2 * i] == 0) {
	    node_cnt[2 * i + 1] = 0;
	    continue;
	 }

	 l = r_lim = c_node_ind[2 * i];
	 r = l_lim = l + node_cnt[2 * i] - 1;

	 do {
	    swap (nodes, l, r);
	    while (l <= l_lim && coord (dir, nodes[l].pos) < cut)
	       l++;
	    while (r >= r_lim && coord (dir, nodes[r].pos) >= cut)
	       r--;
	 }
	 while (l < r);

	 node_cnt[2 * i + 1] = node_cnt[2 * i] - (l - c_node_ind[2 * i]);
	 node_cnt[2 * i] = l - c_node_ind[2 * i];
	 c_node_ind[2 * i + 1] = l;
      }
      child_cnt *= 2;
   }

   true_child_cnt = 0;
   for (i = 0; i < child_cnt; i++) {
      tree_cell *child_ptr;

      if (node_cnt[i] == 0)
	 continue;
      child_ptr = make_tree (level + 1, nodes + c_node_ind[i], node_cnt[i], node_min);
      if (!child_ptr)
	 return NULL;
      t_ptr->child[true_child_cnt++] = child_ptr;
   }
   t_ptr->child_cnt = true_child_cnt;
   return t_ptr;
}

/**************************************************************/

static REAL
coord (int dir, vect pos)
{
   switch (dir) {
   case 1:
      return pos.x1;
   case 2:
      return pos.x2;
   case 3:
      return pos.x3;
   }
   fprintf (stderr, "bad dir %d in coord", dir);
   return 0;
}

/**************************************************************/

static void
swap (work *nodes, long i, long j)
{
   work temp = nodes[i];
   nodes[i] = nodes[j];
   nodes[j] = temp;
}

/**************************************************************/

int
depth (tree_cell *cell_ptr)
{
   int c_ind, c_depth = 0;

   if (cell_ptr->child_cnt == 0)
      return cell_ptr->level;
   for (c_ind = 0; c_ind < cell_ptr->child_cnt; c_ind++) {
      int tmp = depth (cell_ptr->child[c_ind]);
      c_depth = MAX (c_depth, tmp);
   }
   return c_depth;
}

/**************************************************************/

int
tree_cell_cnt (tree_cell *cell_ptr)
{
   int c_ind, ans = 0;

   for (c_ind = 0; c_ind < cell_ptr->child_cnt; c_ind++)
      ans += tree_cell_cnt (cell_ptr->child[c_ind]);
   return ans + 1;
}

/**************************************************************/

REAL
child_parent_ratio (tree_cell *cell_ptr)
{
   REAL ans = 0;
   int c_ind;

   if (cell_ptr->child_cnt == 0)
      return 0.0;
   for (c_ind = 0; c_ind < cell_ptr->child_cnt; c_ind++) {
      REAL val1 = child_parent_ratio (cell_ptr->child[c_ind]);
      REAL val2 = (REAL) cell_ptr->child[c_ind]->node_cnt / cell_ptr->node_cnt;
      REAL val = MAX (val1, val2);
      ans = MAX (ans, val);
   }
   return ans;
}

/**************************************************************/

void
tree_dump (FILE *f_ptr, tree_cell *cell_ptr)
{
   int i, c_ind;

   for (i = 0; i < cell_ptr->level; i++)
      fprintf (f_ptr, "  ");
   fprintf (f_ptr, "node_cnt = %ld, c_cnt = %d\n", cell_ptr->node_cnt, cell_ptr->child_cnt);
   for (c_ind = 0; c_ind < cell_ptr->child_cnt; c_ind++)
      tree_dump (f_ptr, cell_ptr->child[c_ind]);
}

/**************************************************************/

static int
init_tree (void)
{
   if (head_block == NULL) {
      if ((head_block = new_block ()) == NULL) {
	 fprintf (stderr, "new_block failed for init_tree");
	 return 0;
      }
   }
   curr_block = head_block;
   return 0;
}

/**************************************************************/

static tree_cell *
cell_alloc (void)
{
   if (curr_block->cells_used == CELL_BLOCK_SIZE) {
      if (curr_block->next_block == NULL) {
	 if ((curr_block->next_block = new_block ()) == NULL) {
	    fprintf (stderr, "new_block failed for cell_alloc");
	    return NULL;
	 }
      }
      curr_block = curr_block->next_block;
      curr_block->cells_used = 0;
   }
   return curr_block->cells + curr_block->cells_used++;
}

/**************************************************************/

static cell_block *
new_block (void)
{
   cell_block *ans;
   if (!(ans = malloc (sizeof (cell_block)))) {
      fprintf (stderr, "malloc for cell_block failed in new_block");
      return NULL;
   }
   ans->next_block = NULL;
   ans->cells_used = 0;
   return ans;
}

/**************************************************************/

void
free_tree (void)
{
   cell_block *a, *b;

   a = head_block;
   while (a) {
      b = a->next_block;
      free (a);
      a = b;
   }
}
