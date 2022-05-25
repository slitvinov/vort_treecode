/*  CVS:$Id: proj.h,v 1.2 2002/01/21 18:38:39 klindsay Exp $
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

#ifndef PROJ_INCL
#define PROJ_INCL

/**************************************************************/
/* REAL is the floating point type to be used */
#undef REAL
#define REAL double

/* what level run is to be done */
#define NONE 0
#define STARTUP 1
#define RUN_PRINT 2
#define DEBUG NONE

#define ABS(a) ((a) < 0 ? (-(a)) : (a))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* parameters governing dynamic arrays */
#define ARRAY_BLOCK 128

#define MAX_SHEET_CNT 4

/**************************************************************/
/* structure definitions */

typedef struct {
   REAL x1, x2, x3;
} vect;

typedef struct {
   vect pos, v_wght;
   REAL wght_nrm;
} work;

typedef struct {
   REAL th;
   vect pos, arg, new, vel;
   int ins_flag;
} node;

typedef struct {
   REAL alpha, dgamma;
   node *nodes;
   int node_cnt, node_alloc, ins_flag;
} fil;

typedef struct {
   fil *fils;
   int fil_cnt, fil_alloc;
} sheet;

/**************************************************************/
/* global variable declarations */

extern const REAL pi, two_pi, four_pi;

#endif /* if !defined(PROJ_INCL) */
