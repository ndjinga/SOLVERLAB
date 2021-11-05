#!/usr/bin/env python3
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Tests of the library mpi4py from MPI4PY tutorial
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2021
# Description : https://mpi4py.readthedocs.io/en/stable/tutorial.html
#================================================================================================================================

from mpi4py import MPI
import numpy as np

# Tests from MPI4PY tutorial https://mpi4py.readthedocs.io/en/stable/tutorial.html

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

print("My rank is ", rank, " among ", size, "processors ")

###Point-to-Point Communication

#Python objects (pickle under the hood):

if rank == 0:
    data = {'a': 7, 'b': 3.14}
    comm.send(data, dest=1, tag=11)
elif rank == 1:
    data = comm.recv(source=0, tag=11)

#Python objects with non-blocking communication:

if rank == 0:
    data = {'a': 7, 'b': 3.14}
    req = comm.isend(data, dest=1, tag=11)
    req.wait()
elif rank == 1:
    req = comm.irecv(source=0, tag=11)
    data = req.wait()

# passing MPI datatypes explicitly
if rank == 0:
    data = np.arange(1000, dtype='i')
    comm.Send([data, MPI.INT], dest=1, tag=77)
elif rank == 1:
    data = np.empty(1000, dtype='i')
    comm.Recv([data, MPI.INT], source=0, tag=77)

# automatic MPI datatype discovery
if rank == 0:
    data = np.arange(100, dtype=np.float64)
    comm.Send(data, dest=1, tag=13)
elif rank == 1:
    data = np.empty(100, dtype=np.float64)
    comm.Recv(data, source=0, tag=13)

###Collective Communication

#Broadcasting a Python dictionary:

if rank == 0:
    data = {'key1' : [7, 2.72, 2+3j],
            'key2' : ( 'abc', 'xyz')}
else:
    data = None
data = comm.bcast(data, root=0)

#Scattering Python objects:

if rank == 0:
    data = [(i+1)**2 for i in range(size)]
else:
    data = None
data = comm.scatter(data, root=0)
assert data == (rank+1)**2

#Gathering Python objects:

data = (rank+1)**2
data = comm.gather(data, root=0)
if rank == 0:
    for i in range(size):
        assert data[i] == (i+1)**2
else:
    assert data is None

# Broadcasting a NumPy array:

if rank == 0:
    data = np.arange(100, dtype='i')
else:
    data = np.empty(100, dtype='i')
comm.Bcast(data, root=0)
for i in range(100):
    assert data[i] == i

#Scattering NumPy arrays:

sendbuf = None
if rank == 0:
    sendbuf = np.empty([size, 100], dtype='i')
    sendbuf.T[:,:] = range(size)
recvbuf = np.empty(100, dtype='i')
comm.Scatter(sendbuf, recvbuf, root=0)
assert np.allclose(recvbuf, rank)

#Gathering NumPy arrays:

sendbuf = np.zeros(100, dtype='i') + rank
recvbuf = None
if rank == 0:
    recvbuf = np.empty([size, 100], dtype='i')
comm.Gather(sendbuf, recvbuf, root=0)
if rank == 0:
    for i in range(size):
        assert np.allclose(recvbuf[i,:], i)

#Parallel matrix-vector product:


