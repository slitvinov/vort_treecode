.POSIX:
.SUFFIXES: .o
.SUFFIXES: .c
CC = c99
LINK = $(CC)
CFLAGS = -O2 -g
ALL = prog gendat vel

all: $(ALL)

prog_objects = main.o io.o comp_f.o refine.o tree.o disk_stats.o
gendat_objects = gendat.o io.o
vel_objects = vel.o io.o comp_f.o tree.o

prog: $(prog_objects)
	$(LINK) $(LDFLAGS) $(prog_objects) -lm -o prog

gendat: $(gendat_objects)
	$(LINK) $(LDFLAGS) $(gendat_objects) -lm -o gendat

vel: $(vel_objects)
	$(LINK) $(LDFLAGS) $(vel_objects) -lm -o vel

.c.o:; $(CC) -c $(CFLAGS) $<
clean:; rm -f $(prog_objects) $(gendat_objects) $(vel_objects) $(ALL)

include dep.mk
