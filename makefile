#   CVS:$Id: makefile,v 1.2 2002/01/21 18:38:39 klindsay Exp $
#   CVS:$Name: vort_treecode_1_0 $
#
#   Copyright (C) 2002  Keith Lindsay (klindsay@ucar.edu)
#
#   This file is part of vort_treecode, a program to compute the
#   evolution of a three-dimensional vortex sheet using an adaptive
#   treecode to evaluate the velocity field.
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

CC = gcc
CFLAGS = -ansi -Wall -O3

all: prog gendat vel

prog_objects = main.o io.o comp_f.o refine.o tree.o disk_stats.o
gendat_objects = gendat.o io.o
vel_objects = vel.o io.o comp_f.o tree.o

prog: $(prog_objects)
	$(CC) $(CFLAGS) $(prog_objects) -lm -o prog

gendat: $(gendat_objects)
	$(CC) $(CFLAGS) $(gendat_objects) -lm -o gendat

vel: $(vel_objects)
	$(CC) $(CFLAGS) $(vel_objects) -lm -o vel

clean:
	rm -f *.o

comp_f.o: comp_f.c proj.h comp_f.h tree.h
disk_stats.o: disk_stats.c proj.h disk_stats.h
gendat.o: gendat.c proj.h io.h
io.o: io.c proj.h io.h
main.o: main.c proj.h comp_f.h refine.h io.h tree.h disk_stats.h
refine.o: refine.c proj.h refine.h
tree.o: tree.c proj.h tree.h
vel.o: vel.c proj.h comp_f.h io.h tree.h
