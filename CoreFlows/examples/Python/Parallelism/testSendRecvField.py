#!/usr/bin/env python3
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Tests of sending and receiving a MEDCoupling field lying on the same mesh between two processors
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2021
# Description : Use of the parallel Data Exchange Channel of MEDCoupling
#================================================================================================================================

from mpi4py import MPI
import medcoupling as mc
from medcoupling import *
from math import sin, pi

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if(size!=2):
	raise ValueError("Processor ", rank, " : aborting.\n Simulation should done on two processors.\n", size, " processors given")
	
print("My rank is ", rank, " among ", size, "processors ")

procs_source = [0]
procs_target = [1]

interface = mc.CommInterface()
source_group = mc.MPIProcessorGroup(interface, procs_source)
target_group = mc.MPIProcessorGroup(interface, procs_target)
dec = mc.InterpKernelDEC(source_group, target_group)
dec = mc.StructuredCoincidentDEC(source_group, target_group)

# Create a MEDCouplingUMesh from a 3D cartesian mesh
xarr=mc.DataArrayDouble.New(11,1)
xarr.iota(0.)
cmesh=mc.MEDCouplingCMesh.New()
cmesh.setCoords(xarr,xarr,xarr)
mesh=cmesh.buildUnstructured()

#Create a field by application of an analytic function 

if source_group.containsMyRank():
	field=mesh.fillFromAnalytic(ON_CELLS,1,"(x-5.)*(x-5.)+(y-5.)*(y-5.)+(z-5.)*(z-5.)")
	field.setName("SourceField")
else:
	field=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
	field.setMesh(mesh)
	field.setName("TargetField")
	field.fillFromAnalytic(1,"0")
	
dec.attachLocalField(field)
dec.synchronize()

if source_group.containsMyRank():
    dec.sendData()
else:
    dec.recvData()
    #mc.WriteField("toto.med", field, True)
