/*  CVS:$Id: io.c,v 1.2 2002/01/21 18:38:38 klindsay Exp $
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

#include <stdlib.h>
#include <stdio.h>
#include "proj.h"
#include "io.h"

/**************************************************************/

/* headers in data file info */
#define R_HDR_CNT 8
#define I_HDR_CNT 8

typedef struct {
   REAL th;
   vect pos;
} io_struct;

/**************************************************************/
/*
 * File Header
 *   REAL : del2, t, dt0, vel0, dt, vel
 *   int  : i/o version, sheet_cnt
 *
 * Sheet Header
 *   int  : fil_cnt
 *
 * Filament Header
 *   REAL : alpha, dgamma
 *   int  : node_cnt
 *
 * Node
 *   REAL : theta, pos.x1, pos.x2, pos.x3
 *
 */

/**************************************************************/
/* local prototypes */

static int read_sheet (FILE *infile, sheet *a, int *node_alloc_tot);
static int read_fil (FILE *infile, fil *a, int *node_alloc_tot);
static int write_sheet (sheet a, FILE *outfile);
static int write_fil (fil a, FILE *outfile);

/**************************************************************/

int
readin (FILE *infile, REAL *del2, REAL *t0, REAL *dt0, REAL *vel0, REAL *vel, int *sheet_cnt, sheet sheets[MAX_SHEET_CNT],
	work **scratch)
{
   REAL r_hdr[R_HDR_CNT];
   int i_hdr[I_HDR_CNT], sheet_ind, node_alloc_tot = 0;

   if (fread (r_hdr, sizeof (REAL), R_HDR_CNT, infile) != R_HDR_CNT) {
      fprintf (stderr, "fread from infile for r_hdr failed in readin");
      return 1;
   }
   if (fread (i_hdr, sizeof (int), I_HDR_CNT, infile) != I_HDR_CNT) {
      fprintf (stderr, "fread from infile for i_hdr failed in readin");
      return 1;
   }

   *del2 = r_hdr[0];
   *t0 = r_hdr[1];
   *dt0 = r_hdr[2];
   *vel0 = r_hdr[3];
   *vel = r_hdr[4];

   if (i_hdr[0] != 1) {
      fprintf (stderr, "i/o version error in data file");
      return 1;
   }
   *sheet_cnt = i_hdr[1];

   for (sheet_ind = 0; sheet_ind < *sheet_cnt; sheet_ind++)
      if (read_sheet (infile, &sheets[sheet_ind], &node_alloc_tot))
	 return 1;

   fclose (infile);

   if (!(*scratch = malloc (sizeof (work) * node_alloc_tot))) {
      fprintf (stderr, "malloc for scratch failed");
      return 1;
   }

   return 0;
}

/**************************************************************/

int
read_sheet (FILE *infile, sheet *a, int *node_alloc_tot)
{
   REAL r_hdr[R_HDR_CNT];
   int i_hdr[I_HDR_CNT], fil_ind;

   if (fread (r_hdr, sizeof (REAL), R_HDR_CNT, infile) != R_HDR_CNT) {
      fprintf (stderr, "fread from infile for r_hdr failed in read_sheet");
      return 1;
   }
   if (fread (i_hdr, sizeof (int), I_HDR_CNT, infile) != I_HDR_CNT) {
      fprintf (stderr, "fread from infile for i_hdr failed in read_sheet");
      return 1;
   }

   a->fil_cnt = i_hdr[0];
   a->fil_alloc = (2 * a->fil_cnt / ARRAY_BLOCK + 1) * ARRAY_BLOCK;

#if DEBUG == RUN_PRINT
   fprintf (stderr, "fil_cnt = %d, fil_alloc = %d\n", a->fil_cnt, a->fil_alloc);
   fflush (stderr);
#endif

   if (!(a->fils = malloc (sizeof (fil) * a->fil_alloc))) {
      fprintf (stderr, "malloc for fils failed");
      return 1;
   }

   for (fil_ind = 0; fil_ind < a->fil_cnt; fil_ind++)
      if (read_fil (infile, &a->fils[fil_ind], node_alloc_tot))
	 return 1;

   return 0;
}

/**************************************************************/

int
read_fil (FILE *infile, fil *a, int *node_alloc_tot)
{
   REAL r_hdr[R_HDR_CNT];
   int i_hdr[I_HDR_CNT], node_ind;

   if (fread (r_hdr, sizeof (REAL), R_HDR_CNT, infile) != R_HDR_CNT) {
      fprintf (stderr, "fread from infile for r_hdr failed in read_fil");
      return 1;
   }
   if (fread (i_hdr, sizeof (int), I_HDR_CNT, infile) != I_HDR_CNT) {
      fprintf (stderr, "fread from infile for i_hdr failed in read_fil");
      return 1;
   }

   a->alpha = r_hdr[0];
   a->dgamma = r_hdr[1];
   a->node_cnt = i_hdr[0];
   a->node_alloc = (a->node_cnt / ARRAY_BLOCK + 1) * ARRAY_BLOCK;
   *node_alloc_tot += a->node_alloc;

   if (!(a->nodes = malloc (sizeof (node) * a->node_alloc))) {
      fprintf (stderr, "malloc for nodes failed");
      return 1;
   }

   for (node_ind = 0; node_ind < a->node_cnt; node_ind++) {
      io_struct to_read;
      if (fread (&to_read, sizeof (io_struct), 1, infile) != 1) {
	 fprintf (stderr, "fread from infile for io_struct failed");
	 return 1;
      }
      a->nodes[node_ind].th = to_read.th;
      a->nodes[node_ind].pos = to_read.pos;
   }

   return 0;
}

/**************************************************************/

int
dump (REAL del2, REAL t, REAL dt0, REAL vel0, REAL vel, int sheet_cnt, sheet sheets[MAX_SHEET_CNT], FILE *outfile)
{
   REAL r_hdr[R_HDR_CNT];
   int i_hdr[I_HDR_CNT];
   int i, sheet_ind;

   for (i = 0; i < R_HDR_CNT; i++)
      r_hdr[i] = 0;
   for (i = 0; i < I_HDR_CNT; i++)
      i_hdr[i] = 0;
   r_hdr[0] = del2;
   r_hdr[1] = t;
   r_hdr[2] = dt0;
   r_hdr[3] = vel0;
   r_hdr[4] = vel;
   i_hdr[0] = 1;
   i_hdr[1] = sheet_cnt;

   if (fwrite (r_hdr, sizeof (REAL), R_HDR_CNT, outfile) != R_HDR_CNT) {
      fprintf (stderr, "fwrite to outfile for r_hdr failed");
      return 1;
   }
   if (fwrite (i_hdr, sizeof (int), I_HDR_CNT, outfile) != I_HDR_CNT) {
      fprintf (stderr, "fwrite to outfile for i_hdr failed");
      return 1;
   }

   for (sheet_ind = 0; sheet_ind < sheet_cnt; sheet_ind++)
      if (write_sheet (sheets[sheet_ind], outfile))
	 return 1;

   return 0;
}

/**************************************************************/

static int
write_sheet (sheet a, FILE *outfile)
{
   REAL r_hdr[R_HDR_CNT];
   int i_hdr[I_HDR_CNT];
   int i, fil_ind;

   for (i = 0; i < R_HDR_CNT; i++)
      r_hdr[i] = 0;
   for (i = 0; i < I_HDR_CNT; i++)
      i_hdr[i] = 0;
   i_hdr[0] = a.fil_cnt;

   if (fwrite (r_hdr, sizeof (REAL), R_HDR_CNT, outfile) != R_HDR_CNT) {
      fprintf (stderr, "fwrite to outfile for r_hdr failed");
      return 1;
   }
   if (fwrite (i_hdr, sizeof (int), I_HDR_CNT, outfile) != I_HDR_CNT) {
      fprintf (stderr, "fwrite to outfile for i_hdr failed");
      return 1;
   }

   for (fil_ind = 0; fil_ind < a.fil_cnt; fil_ind++)
      if (write_fil (a.fils[fil_ind], outfile))
	 return 1;

   return 0;
}

/**************************************************************/

static int
write_fil (fil a, FILE *outfile)
{
   REAL r_hdr[R_HDR_CNT];
   int i_hdr[I_HDR_CNT];
   int i, node_ind;

   for (i = 0; i < R_HDR_CNT; i++)
      r_hdr[i] = 0;
   for (i = 0; i < I_HDR_CNT; i++)
      i_hdr[i] = 0;
   r_hdr[0] = a.alpha;
   r_hdr[1] = a.dgamma;
   i_hdr[0] = a.node_cnt;

   if (fwrite (r_hdr, sizeof (REAL), R_HDR_CNT, outfile) != R_HDR_CNT) {
      fprintf (stderr, "fwrite to outfile for r_hdr failed");
      return 1;
   }
   if (fwrite (i_hdr, sizeof (int), I_HDR_CNT, outfile) != I_HDR_CNT) {
      fprintf (stderr, "fwrite to outfile for i_hdr failed");
      return 1;
   }

   for (node_ind = 0; node_ind < a.node_cnt; node_ind++) {
      io_struct to_write;
      to_write.th = a.nodes[node_ind].th;
      to_write.pos = a.nodes[node_ind].pos;
      if (fwrite (&to_write, sizeof (io_struct), 1, outfile) != 1) {
	 fprintf (stderr, "fwrite to outfile for io_struct failed");
	 return 1;
      }
   }

   return 0;
}
