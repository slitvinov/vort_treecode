CC = c99
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
