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


nodes = []
with open(sys.argv[1], 'rb') as f:
    del2, t, dt0, vel0, dt, vel, *rest = read_real(f, 8)
    version, sheet_cnt, *rest = read_int(f, 8)
    for sheet_ind in range(sheet_cnt):
        rest = read_real(f, 8)
        fil_cnt, *rest = read_int(f, 8)
        for fil_ind in range(fil_cnt):
            alpha, dgamma, *rest = read_real(f, 8)
            node_cnt, *rest = read_int(f, 8)
            buf = f.read(struct.calcsize(4 * REAL * node_cnt))
            node = np.ndarray((node_cnt, 4), REAL, buf)
            nodes.append(node)

nnodes = sum(len(node) for node in nodes)
nbytes = sum(3 * node.nbytes // 4 for node in nodes)

with open(sys.argv[2], 'wb+') as f:
    f.write(b"""\
# vtk DataFile Version 2.0
generate with tool/read.py
BINARY
DATASET POLYDATA
POINTS %d double
""" % nnodes)
    offset = f.tell()
    f.seek(nbytes - 1, 1)
    f.write(b'\0')

pos = np.memmap(sys.argv[2], '>' + REAL, 'r+', offset, (nnodes, 3))
shift = 0
for node in nodes:
    np.copyto(pos[shift:shift + len(node), :], node[:, 1:])
    shift += len(node)
