/* generate data for a two disks of radius 1 and fil_cnt filaments each */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "proj.h"
#include "io.h"

const REAL pi = 3.1415926535897932384626433832795;
const REAL two_pi = 6.2831853071795864769252867665590;
const REAL four_pi = 12.566370614359172953850573533118;

/**************************************************************/
/* local prototypes */

static int init (int fil_cnt, sheet sheets[MAX_SHEET_CNT], REAL ang);

/**************************************************************/

int
main (int argc, char **argv)
{
   FILE *outfile;
   sheet sheets[MAX_SHEET_CNT];
   REAL del, dt0, ang;
   int fil_cnt;

   if (argc != 6) {
      fprintf (stderr, "usage should be : gendat <del> <dt0> <fil_cnt> <ang> <outfile>");
      return EXIT_FAILURE;
   }

   del = atof (*++argv);
   dt0 = atof (*++argv);
   fil_cnt = atoi (*++argv);
   ang = atof (*++argv);
   printf ("ang = %f\n", ang);
   if (!(outfile = fopen (*++argv, "wb"))) {
      fprintf (stderr, "fopen for outfile on '%s' failed", *argv);
      return EXIT_FAILURE;
   }

   if (init (fil_cnt, sheets, ang))
      return EXIT_FAILURE;

   return dump (del * del, 0.0, dt0, 0.0, 0.0, 2, sheets, outfile) ? EXIT_FAILURE : EXIT_SUCCESS;
}

/**************************************************************/

static int
init (int fil_cnt, sheet sheets[MAX_SHEET_CNT], REAL ang)
{
   REAL cosang = cos (ang), sinang = sin (ang);
   void *ptr;
   int fil_ind;

   sheets[0].fil_cnt = sheets[1].fil_cnt = fil_cnt;
   if (!(ptr = malloc (sizeof (fil) * 2 * fil_cnt))) {
      fprintf (stderr, "malloc for fils failed");
      return 1;
   }
   sheets[0].fils = ptr;
   sheets[1].fils = sheets[0].fils + fil_cnt;

   for (fil_ind = 0; fil_ind < fil_cnt; fil_ind++) {
      fil *f0 = sheets[0].fils + fil_ind, *f1 = sheets[1].fils + fil_ind;
      node *n0, *n1;
      REAL alpha, r;
      int N, i;

      f0->alpha = f1->alpha = alpha = pi / 2 * fil_ind / (fil_cnt - 1);
      f0->dgamma = f1->dgamma = -sin (alpha);
      r = sin (alpha);
#if 1
      N = 128 * (1 + r);
      if (N % 8)
	 N += 8 - N % 8;
#else
      N = 256;
#endif
      if (fil_ind == 0)
	 N = 4;
      f0->node_cnt = f1->node_cnt = N;
      if (!(ptr = malloc (sizeof (node) * 2 * N))) {
	 fprintf (stderr, "malloc for nodes failed");
	 return 1;
      }
      f0->nodes = n0 = ptr;
      f1->nodes = n1 = f0->nodes + N;

      /* compute starting positions */
      for (i = 0; i < N; i++) {
	 REAL th = (two_pi * i) / N;

	 n0[i].th = th;
	 n0[i].pos.x1 = 1 + r * cos (th) * cosang;
	 n0[i].pos.x2 = r * sin (th);
	 n0[i].pos.x3 = -r * cos (th) * sinang;

	 n1[i].th = th;
	 n1[i].pos.x1 = -1 - r * cos (th) * cosang;
	 n1[i].pos.x2 = -r * sin (th);
	 n1[i].pos.x3 = -r * cos (th) * sinang;
      }
   }

   return 0;
}
