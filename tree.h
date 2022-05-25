/*  CVS:$Id: tree.h,v 1.2 2002/01/21 18:38:39 klindsay Exp $
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

#ifndef TREE_INCL
#define TREE_INCL

#include "proj.h"

/* various maximums & minimums */
#define MAX_P 8
#define MAX_CHILD_CNT 8

/* what type of error bound is to be used */
#define POTENTIAL 0
#define VELOCITY 1
#define BOUND POTENTIAL

/**************************************************************/
/* structure definitions */

struct str_tree_cell {
   struct str_tree_cell *child[MAX_CHILD_CNT];
   vect exp_pt, width_2, mmnts[MAX_P][MAX_P][MAX_P];
   REAL r, wght_nrm[MAX_P + 1];
   work *nodes;
   long node_cnt;
   int level, mmnts_stat, child_cnt;
};

typedef struct str_tree_cell tree_cell;

/**************************************************************/
/* external function prototypes */

tree_cell *make_tree (int level, work *nodes, long node_cnt, int node_min);
int depth (tree_cell *cell_ptr);
int tree_cell_cnt (tree_cell *cell_ptr);
REAL child_parent_ratio (tree_cell *cell_ptr);
void tree_dump (FILE *f_ptr, tree_cell *cell_ptr);
void free_tree (void);

#endif /* if !defined(TREE_INCL) */
