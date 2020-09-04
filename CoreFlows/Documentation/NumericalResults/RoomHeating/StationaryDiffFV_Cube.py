# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution VF de l'équation de Laplace 3D -\Delta T = 0 avec conditions aux limites de Dirichlet u non nulle (fenetre et radiateur)
# Authors     : Michaël Ndjinga, Sédrick Kameni Ngwamou
# Copyright   : CEA Saclay 2019
# Description : Utilisation de la méthode des volumes finis avec champs u discrétisés aux cellules d'un maillage de cubes
#               Conditions limites correspondant au refroidissement dû à une fenêtre et au chauffage dû à un radiateur
#				Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant la librairie CDMATH
#================================================================================================================================

import CoreFlows as cf
import cdmath
import PV_routines
import VTK_routines

def StationaryDiffusionEquation_3DVF_RoomCooling(resolution):
	spaceDim = 3;
	
	#Chargement du maillage cartésien du domaine
	#==============================================
	my_mesh = cdmath.Mesh('RoomCoulingCubeBon.med')#Rom_CoulingCube2.med
	
	nbCells = my_mesh.getNumberOfCells()
	print "Loaded Structured 3D mesh"
	
	#Conditions limites
	Tmur=20
	Tfenetre=0
	Tradiateur=40

	FEComputation=False
	myProblem = cf.StationaryDiffusionEquation(spaceDim,FEComputation);
	myProblem.setMesh(my_mesh);
	
	myProblem.setDirichletBoundaryCondition("Fenetre",Tfenetre)
	myProblem.setDirichletBoundaryCondition("Radiateur_Dessous",Tradiateur)
	myProblem.setDirichletBoundaryCondition("Radiateur_Devant",Tmur)
	myProblem.setDirichletBoundaryCondition("Radiateur_Droit",Tmur)
	myProblem.setDirichletBoundaryCondition("Murs",Tmur)

    # name of result file
	fileName = "3DVF_StructuredCubes"+str(nbCells);

    # computation parameters
	myProblem.setFileName(fileName);

    # Run the computation
	myProblem.initialize();
	print("Running python "+ fileName );

	ok = myProblem.solveStationaryProblem();
	if (not ok):
		print( "Python simulation of " + fileName + "  failed ! " );
		pass
	else:
		print( "Python simulation of " + fileName + "  successful ! " );
		pass

	myProblem.terminate();
	
    # save 2D picture
	VTK_routines.Clip_VTK_data_to_VTK("StationaryDiffusionEquation_3DVF_StructuredCubes"+str(nbCells)+'_0.vtu',"Clip_VTK_data_to_VTK_"+ "StationaryDiffusionEquation_3DVF_StructuredCubes_"+str(nbCells)+'_0.vtu',[0.5,0.5,0.5], [-0.5,-0.5,-0.5],resolution )
	PV_routines.Save_PV_data_to_picture_file("Clip_VTK_data_to_VTK_"+"StationaryDiffusionEquation_3DVF_StructuredCubes_"+str(nbCells)+'_0.vtu',"Temperature",'NODES',"Clip_VTK_data_to_VTK_"+"StationaryDiffusionEquation_3DVF_StructuredCubes_"+str(nbCells))
	return ok

if __name__ == """__main__""":
    StationaryDiffusionEquation_3DVF_RoomCooling(100)
