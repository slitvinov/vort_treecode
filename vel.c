/*  CVS:$Id: vel.c,v 1.2 2002/01/21 18:38:39 klindsay Exp $
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
#include <stdlib.h>
#include <stdio.h>
#include "proj.h"
#include "comp_f.h"
#include "io.h"
#include "tree.h"

#define DT_STATS 1

/**************************************************************/
/* global variable definitions */

const REAL pi = 3.1415926535897932384626433832795;
const REAL two_pi = 6.2831853071795864769252867665590;
const REAL four_pi = 12.566370614359172953850573533118;

/**************************************************************/
/* local prototypes */

static int get_args (char **argv, FILE **infile, REAL *tol, int *node_min, FILE **outfile);
static vect v_add (vect a, REAL beta, vect b);
static void free_sheets (int sheet_cnt, sheet sheets[MAX_SHEET_CNT]);

/**************************************************************/

int
main (int argc, char *argv[])
{
   FILE *infile, *outfile;
   sheet sheets[MAX_SHEET_CNT];
   work *scratch;
   node *n_ptr;
   REAL tol, del2, t, dt0, vel0, vel;
   int node_min, sheet_cnt, sheet_ind, status;

   if (argc != 5) {
      fprintf (stderr, "usage should be : prog <infile> <tol> <node_min> <outfile>");
      return EXIT_FAILURE;
   }
   if (get_args (argv, &infile, &tol, &node_min, &outfile))
      return EXIT_FAILURE;
   if (readin (infile, &del2, &t, &dt0, &vel0, &vel, &sheet_cnt, sheets, &scratch))
      return EXIT_FAILURE;

#if DEBUG == STARTUP
   status = dump (del2, t, dt0, vel0, vel, sheet_cnt, sheets, outfile);
   free_tree ();
   free (scratch);
   free_sheets (sheet_cnt, sheets);
   return status ? EXIT_FAILURE : EXIT_SUCCESS;
#endif

   for (sheet_ind = 0; sheet_ind < sheet_cnt; sheet_ind++) {
      int i, fil_cnt = sheets[sheet_ind].fil_cnt;
      fil *fil_ptr = sheets[sheet_ind].fils;
      for (i = 0; i < fil_cnt; fil_ptr++, i++) {
	 int j, N = fil_ptr->node_cnt;
	 for (n_ptr = fil_ptr->nodes, j = 0; j < N; n_ptr++, j++) {
	    n_ptr->arg = n_ptr->pos;
	 }
      }
   }
   f (del2, sheet_cnt, sheets, scratch, tol, node_min);
   for (sheet_ind = 0; sheet_ind < sheet_cnt; sheet_ind++) {
      int i, fil_cnt = sheets[sheet_ind].fil_cnt;
      fil *fil_ptr = sheets[sheet_ind].fils;
      for (i = 0; i < fil_cnt; fil_ptr++, i++) {
	 int j, N = fil_ptr->node_cnt;
	 for (n_ptr = fil_ptr->nodes, j = 0; j < N; n_ptr++, j++) {
	    n_ptr->pos = v_add (n_ptr->pos, dt0, n_ptr->vel);
	 }
      }
   }

   status = dump (del2, t, dt0, vel0, vel, sheet_cnt, sheets, outfile);
   free_tree ();
   free (scratch);
   free_sheets (sheet_cnt, sheets);
   return status ? EXIT_FAILURE : EXIT_SUCCESS;
}

/**************************************************************/
/* order of arguments is :
 *    infile, tol, node_min, outfile
 */

static int
get_args (char **argv, FILE **infile, REAL *tol, int *node_min, FILE **outfile)
{
   if (!(*infile = fopen (*++argv, "rb"))) {
      fprintf (stderr, "fopen for infile on '%s' failed", *argv);
      return 1;
   }

   *tol = atof (*++argv);
   *node_min = atoi (*++argv);

   if (!(*outfile = fopen (*++argv, "wb"))) {
      fprintf (stderr, "fopen for outfile on '%s' failed", *argv);
      return 1;
   }

   return 0;
}

/**************************************************************/

static vect
v_add (vect a, REAL beta, vect b)
{
   vect ans;

   ans.x1 = a.x1 + beta * b.x1;
   ans.x2 = a.x2 + beta * b.x2;
   ans.x3 = a.x3 + beta * b.x3;

   return ans;
}

/**************************************************************/

static void
free_sheets (int sheet_cnt, sheet sheets[MAX_SHEET_CNT])
{
   int sheet_ind;

   for (sheet_ind = 0; sheet_ind < sheet_cnt; sheet_ind++) {
      int fil_ind;
      for (fil_ind = 0; fil_ind < sheets[sheet_ind].fil_cnt; fil_ind++) {
	 free (sheets[sheet_ind].fils[fil_ind].nodes);
      }
      free (sheets[sheet_ind].fils);
   }
}
