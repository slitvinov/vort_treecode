import sys
import struct
import numpy as np

REAL = 'd'


def read_real(f, n):
    buf = f.read(struct.calcsize(REAL * n))
    return struct.unpack(REAL * n, buf)


def read_int(f, n):
    buf = f.read(struct.calcsize('i' * n))
    return struct.unpack('i' * n, buf)


with open(sys.argv[1], 'rb') as f:
    del2, t, dt0, vel0, dt, vel, *rest = read_real(f, 8)
    version, sheet_cnt, *rest = read_int(f, 8)
    for sheet_ind in range(sheet_cnt):
        rest = read_real(f, 8)
        fil_cnt, *rest = read_int(f, 8)
        for fil_ind in range(fil_cnt):
            alpha, dgamma, *rest = read_real(f, 8)
            node_cnt, *rest = read_int(f, 8)
            sys.stderr.write("node_cnt: %d\n" % node_cnt)
            buf = f.read(struct.calcsize(4 * REAL * node_cnt))
            node = np.ndarray((node_cnt, 4), REAL, buf)
            # for t, x, y, z in node:
            #    print(x, y, z)
