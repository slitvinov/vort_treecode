.POSIX:
.SUFFIXES: .o
.SUFFIXES: .c
CC = c99
LINK = $(CC)
CFLAGS = -O2 -g

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

.c.o:; $(LINK) -c $(LDFLAGS) $<
clean:; rm -f *.o

include dep.mk
