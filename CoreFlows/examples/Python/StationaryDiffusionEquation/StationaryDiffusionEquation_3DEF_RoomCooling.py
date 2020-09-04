#!/usr/bin/env python
# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution EF de l'équation de Laplace 3D -\Delta T = 0 avec conditions aux limites de Dirichlet u non nulle (fenetre et radiateur)
# Authors     : Michaël Ndjinga, Sédrick Kameni Ngwamou
# Copyright   : CEA Saclay 2019
# Description : Utilisation de la méthode des éléménts finis P1 avec champs u discrétisés aux noeuds d'un maillage tétraédrique
#               Condition limites correspondant au refroidissement dû à une fenêtre et au chauffage dû à un radiateur
#				Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant la librairie CDMATH
#================================================================================================================================

import CoreFlows as cf
import cdmath

def StationaryDiffusionEquation_3DEF_RoomCooling():
	spaceDim = 3;
	
	#Chargement du maillage tétraédrique du domaine
	#==============================================
	my_mesh = cdmath.Mesh("../resources/RoomWithTetras2488.med")
	
	print "Loaded unstructured 3D mesh"
	
	#Conditions limites
	Tmur=20
	Tfenetre=0
	Tradiateur=40

	FEComputation=True
	myProblem = cf.StationaryDiffusionEquation(spaceDim,FEComputation);
	myProblem.setMesh(my_mesh);
	
	myProblem.setDirichletBoundaryCondition("Fenetre",Tfenetre)
	myProblem.setDirichletBoundaryCondition("Radiateur_sous_fenetre",Tradiateur)
	myProblem.setDirichletBoundaryCondition("Radiateur_Devant",Tmur)
	myProblem.setDirichletBoundaryCondition("Radiateur_droit",Tmur)
	myProblem.setDirichletBoundaryCondition("Murs",Tmur)

    # name of result file
	fileName = "StationnaryDiffusion_3DEF_UnstructuredTetrahedra";

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
	return ok

if __name__ == """__main__""":
    StationaryDiffusionEquation_3DEF_RoomCooling()
