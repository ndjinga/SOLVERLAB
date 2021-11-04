#!/usr/bin/env python3
# -*-coding:utf-8 -*

#===============================================================================================================================
# Name        : Tests of sending and receiving 2D MEDCoupling fields on nodes (P1) lying on different meshes between two processors
# Author      : Michaël Ndjinga
# Copyright   : CEA Saclay 2021
# Description : Use of the parallel Data Exchange Channel InterpKernelDEC of MEDCoupling
#================================================================================================================================

from mpi4py import MPI
import medcoupling as mc

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if(size!=2):
	raise ValueError("Processor ", rank, " : aborting.\n Simulation should done on two processors.\n", size, " processors given")
	
print("My rank is ", rank, " among ", size, "processors")

procs_source = [0]
procs_target = [1]

interface = mc.CommInterface()
source_group = mc.MPIProcessorGroup(interface, procs_source)
target_group = mc.MPIProcessorGroup(interface, procs_target)
dec = mc.InterpKernelDEC(source_group, target_group)

# Create a MEDCouplingUMesh from a 3D cartesian mesh
xarr=mc.DataArrayDouble.New(11,1)
xarr.iota(0.)
cmesh=mc.MEDCouplingCMesh.New()
cmesh.setCoords(xarr,xarr)
mesh=cmesh.buildUnstructured()
mesh.setName("RegularSquare")
mesh.simplexize(rank)#The squares are cut in two right triangles according to one of the two possible diagonals

#Create a field by application of an analytic function 
if source_group.containsMyRank():
	field=mesh.fillFromAnalytic(mc.ON_NODES,1,"(x-5.)*(x-5.)+(y-5.)*(y-5.)")
	field.setName("SourceField")
	#field.setNature(mc.ExtensiveConservation)
	field.setNature(mc.IntensiveMaximum)
	mc.WriteField("source_field.med", field, True)
	print("Processor ", rank, " has created and saved the source field")
else:
	field=mesh.fillFromAnalytic(mc.ON_NODES,1,"0")
	field.setName("TargetField")
	#field.setNature(mc.ExtensiveConservation)
	field.setNature(mc.IntensiveMaximum)
	print("Processor ", rank, " has created the target field")
	
dec.setMethod("P1")
dec.attachLocalField(field)
dec.synchronize()

if source_group.containsMyRank():
	dec.sendData()
	print("Processor ", rank, " has sent the source field")
else:
	dec.recvData()
	print("Processor ", rank, " has received the source field on the target mesh")
	exact_field=mesh.fillFromAnalytic(mc.ON_NODES,1,"(x-5.)*(x-5.)+(y-5.)*(y-5.)")
	exact_field.setName("ExactField")
	#Computing maximum error
	coordsArr = mesh.getCoords()
	values=field.getValueOnMulti(coordsArr)
	exact_values=exact_field.getValueOnMulti(coordsArr)
	error=(values-exact_values).normMax()/exact_values.normMax()
	print("Processor ", rank, " received source field differs from theoretical value by less than", error, " (maximum relative norm on node values)" )
	assert error<1.e-1
	mc.WriteField("target_field.med", field, True)
	mc.WriteField("exact_field.med", exact_field, True)
