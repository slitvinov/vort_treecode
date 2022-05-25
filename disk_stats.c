#include <stdio.h>
#include <math.h>
#include "proj.h"
#include "disk_stats.h"

/* print statistics for a single sheet
 * only intended for a circular disk
 * whose filament are mostly parallel to the xy-plane
 */
void
disk_stats (sheet *sheet_ptr)
{
   int fil_ind;

   fprintf (stdout, "fil_cnt = %d\n", sheet_ptr->fil_cnt);
   for (fil_ind = 0; fil_ind < sheet_ptr->fil_cnt; fil_ind++) {
      fil *fil_ptr = sheet_ptr->fils + fil_ind;
      node *n_ptr = fil_ptr->nodes;
      int N = fil_ptr->node_cnt, i;
      REAL wght, r_avg, z_avg, d1, d2, d3, var;

      r_avg = z_avg = 0;
      for (i = 0; i < N; i++) {
	 if (i == 0) {
	    wght = 0.5 * (n_ptr[1].th - (n_ptr[N - 1].th - two_pi));
	 } else if (i == N - 1) {
	    wght = 0.5 * ((n_ptr[0].th + two_pi) - n_ptr[N - 2].th);
	 } else {
	    wght = 0.5 * (n_ptr[i + 1].th - n_ptr[i - 1].th);
	 }
	 r_avg += wght * sqrt (n_ptr[i].pos.x1 * n_ptr[i].pos.x1 + n_ptr[i].pos.x2 * n_ptr[i].pos.x2);
	 z_avg += wght * n_ptr[i].pos.x3;
      }
      r_avg /= two_pi;
      z_avg /= two_pi;

      var = 0;
      for (i = 0; i < N; i++) {
	 if (i == 0) {
	    wght = 0.5 * (n_ptr[1].th - (n_ptr[N - 1].th - two_pi));
	 } else if (i == N - 1) {
	    wght = 0.5 * ((n_ptr[0].th + two_pi) - n_ptr[N - 2].th);
	 } else {
	    wght = 0.5 * (n_ptr[i + 1].th - n_ptr[i - 1].th);
	 }
	 d1 = n_ptr[i].pos.x1 - r_avg * cos (n_ptr[i].th);
	 d2 = n_ptr[i].pos.x2 - r_avg * sin (n_ptr[i].th);
	 d3 = n_ptr[i].pos.x3 - z_avg;
	 var += wght * (d1 * d1 + d2 * d2 + d3 * d3);
      }
      var = sqrt (var / two_pi);
      fprintf (stdout, "d_s: %.3e %.3e %.3e %.3e %.3e\n", fil_ptr->alpha, fil_ptr->dgamma, r_avg, z_avg, var);
   }
}
