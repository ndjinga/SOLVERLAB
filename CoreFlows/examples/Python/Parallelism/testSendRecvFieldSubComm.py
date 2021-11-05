#!/usr/bin/env python3
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Tests of using a subcommnicator for sending and receiving a 3D MEDCoupling field on cells (P0) lying on the same mesh between two processors
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2021
# Description : Use of the parallel Data Exchange Channel StructuredCoincidentDEC of MEDCoupling
#================================================================================================================================

from mpi4py import MPI
import medcoupling as mc

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if(size!=3):
	raise ValueError("Processor ", rank, " : aborting.\n Simulation should done on three processors.\n", size, " processors given")
	
print("My rank is ", rank, " among ", size, "processors")

procs_source = [0]
procs_target = [1]
procs_idle = [2]

interface = mc.CommInterface()
source_group = mc.MPIProcessorGroup(interface, procs_source)
target_group = mc.MPIProcessorGroup(interface, procs_target)
dec = mc.StructuredCoincidentDEC(source_group, target_group)

# Create a MEDCouplingUMesh from a 3D cartesian mesh
xarr=mc.DataArrayDouble.New(11,1)
xarr.iota(0.)
cmesh=mc.MEDCouplingCMesh.New()
cmesh.setCoords(xarr,xarr,xarr)
mesh=cmesh.buildUnstructured()
mesh.setName("RegularSquare")

#Create a field by application of an analytic function 
if source_group.containsMyRank():
	field=mesh.fillFromAnalytic(mc.ON_CELLS,1,"(x-5.)*(x-5.)+(y-5.)*(y-5.)+(z-5.)*(z-5.)")
	field.setName("SourceField")
	mc.WriteField("source_field.med", field, True)
	print("Processor ", rank, " has created and saved the source field")
else:
	field=mesh.fillFromAnalytic(mc.ON_CELLS,1,"0")
	field.setName("TargetField")
	print("Processor ", rank, " has created the target field")
	
dec.attachLocalField(field)
dec.synchronize()

if source_group.containsMyRank():
	dec.sendData()
	print("Processor ", rank, " has sent the source field")
elif target_group.containsMyRank():
	dec.recvData()
	print("Processor ", rank, " has received the source field on the target mesh")
	exact_field=mesh.fillFromAnalytic(mc.ON_CELLS,1,"(x-5.)*(x-5.)+(y-5.)*(y-5.)+(z-5.)*(z-5.)")
	exact_field.setName("ExactField")
	error=(field-exact_field).normL2()[0]
	print("Processor ", rank, " received source field differs from theoretical value by ", error, " (L2 norm on cells)" )
	assert abs(error)<1.e-6
	mc.WriteField("target_field.med", field, True)
	mc.WriteField("exact_field.med", exact_field, True)
else:
	print("Processor ", rank, " did nothing" )
