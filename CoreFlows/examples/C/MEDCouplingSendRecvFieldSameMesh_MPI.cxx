//============================================================================
// Name        : Tests of using a subcommnicator for sending and receiving a 3D MEDCoupling field on cells (P0) lying on the same mesh between two groups of two processors
// Author      : Michael NDJINGA
// Date        : November 2021
// Description : Use of the parallel Data Exchange Channel StructuredCoincidentDEC of MEDCoupling
//============================================================================

#include <iostream>
#include <string>
#include <set>

#include "StructuredCoincidentDEC.hxx"
#include "MPIProcessorGroup.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDLoader.hxx"

#include <mpi.h>
#include <cassert>

using namespace std;

 
int main(int argc, char *argv[])
{
	/* MPI initialisation */
	MPI_Init(&argc, &argv);
	int    size;        /* size of communicator */
	int    rank;        /* processor rank */
	int sub_rank, sub_size;/* rank in subcommunicator */
	int color;/* tells if I belong to the sub_communicator */
	MPI_Comm sub_comm ;/*subcommunicator will be used in exchanges */
	std::set<int> procs_source, procs_target;/* group of procs that will send or receive data */
	MEDCoupling::MEDCouplingFieldDouble * field;/*field used to send or receive data */
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	if(size!=4)
		printf("Processor %d : aborting.\n Simulation should done on four processors.\n %d processors given\n",rank,size);
	
	color=rank/2;
		
	printf("My rank is %d among %d processors, my color is %d\n",rank, size,color);
	
	MPI_Comm_split(MPI_COMM_WORLD, color, rank, &sub_comm);	/* two groups (0,1) and (2,3) */
	MPI_Comm_rank(sub_comm, &sub_rank);
	MPI_Comm_size(sub_comm, &sub_size);
	
	printf("WORLD RANK/SIZE: %d/%d \t subcommunicator RANK/SIZE: %d/%d\n",	rank, size, sub_rank, sub_size);

	
	procs_source.insert(0);/* sub rank 0 will send data */
	procs_target.insert(1);/* sub rank 1 will receive data */

	MEDCoupling::CommInterface interface = MEDCoupling::CommInterface();
	MEDCoupling::MPIProcessorGroup source_group = MEDCoupling::MPIProcessorGroup(interface, procs_source,sub_comm);
	MEDCoupling::MPIProcessorGroup target_group = MEDCoupling::MPIProcessorGroup(interface, procs_target,sub_comm);
	MEDCoupling::StructuredCoincidentDEC dec = MEDCoupling::StructuredCoincidentDEC(source_group, target_group);

	//Create a MEDCouplingUMesh from a 2D cartesian mesh
	MEDCoupling::DataArrayDouble * xarr=MEDCoupling::DataArrayDouble::New();
	xarr->alloc(11,1);
	xarr->iota(0.);
	MEDCoupling::MEDCouplingCMesh * cmesh=MEDCoupling::MEDCouplingCMesh::New("My2D_CMesh");
	cmesh->setCoords(xarr,xarr,xarr);
	MEDCoupling::MEDCouplingUMesh * mesh=cmesh->buildUnstructured();
	mesh->setName("RegularSquare");

	if(sub_rank == 0)
	{
		field = mesh->fillFromAnalytic(MEDCoupling::ON_CELLS,1,"(x-5.)*(x-5.)+(y-5.)*(y-5.)+(z-5.)*(z-5.)");
		field->setName("SourceField");
		MEDCoupling::WriteField("source_field"+to_string(rank)+".med", field, true);
		printf("Processor with global rank %d has created and saved the source field\n", rank);
	}
	else
	{
		field=mesh->fillFromAnalytic(MEDCoupling::ON_CELLS,1,"0");
		field->setName("TargetField");
		printf("Processor with global rank %d has created the target field\n", rank);
	}
	
	dec.attachLocalField(field);
	dec.synchronize();
	
	if(sub_rank == 0)
	{
		dec.sendData();
		printf("Processor with global rank %d has sent the source field\n", rank);
	}
	else
	{
		dec.recvData();
		printf("Processor with global rank %d has received the source field on the target mesh\n", rank);
		/* Solve the bug in StructuredCoincidentDEC then uncomment the lines below to check the result */
		//MEDCoupling::MEDCouplingFieldDouble * exact_field=mesh->fillFromAnalytic(MEDCoupling::ON_CELLS,1,"(x-5.)*(x-5.)+(y-5.)*(y-5.)+(z-5.)*(z-5.)");
		//exact_field->setName("ExactField");
  		//double error=((*field)-(*exact_field))->normMax(0)/exact_field->normMax(0);
		//printf("Processor with global rank %d received source field that differs from theoretical value by %d (maximum relative norm on cells)\n", rank, error );
		//assert( fabs(error)<1.e-6 );
		//MEDCoupling::WriteField("target_field"+to_string(rank)+".med", field, true);
		//MEDCoupling::WriteField("exact_field"+to_string(rank)+".med", exact_field, true);	
	}

	MPI_Comm_free(&sub_comm);
	MPI_Finalize();
    return 0;
}
