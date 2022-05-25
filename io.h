/*  CVS:$Id: io.h,v 1.2 2002/01/21 18:38:39 klindsay Exp $
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

#ifndef IO_H
#define IO_H

#include "proj.h"

int readin (FILE *infile, REAL *del2, REAL *t0, REAL *dt0, REAL *vel0, REAL *vel, int *sheet_cnt, sheet sheets[MAX_SHEET_CNT],
	    work **scratch);
int dump (REAL del2, REAL t, REAL dt0, REAL vel0, REAL vel, int sheet_cnt, sheet sheets[MAX_SHEET_CNT], FILE *outfile);

#endif /* !defined(IO_H) */