#!/usr/bin/env python
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

def StationaryDiffusionEquation_3DVF_RoomCooling_StructuredCubes():
	spaceDim = 3;
	
	#Chargement du maillage cartésien du domaine
	#==============================================
	my_mesh = cdmath.Mesh("../resources/RoomWithCubes480.med")
	
	print "Loaded Structured 3D mesh"
	
	#Conditions limites
	Tmur=20
	Tfenetre=0
	Tradiateur=40

	FEComputation=False
	myProblem = cf.StationaryDiffusionEquation(spaceDim,FEComputation);
	myProblem.setMesh(my_mesh);
	
	myProblem.setDirichletBoundaryCondition("Fenetre",Tfenetre)
	myProblem.setDirichletBoundaryCondition("Radiateur_sous_fenetre",Tradiateur)
	myProblem.setDirichletBoundaryCondition("Radiateur_devant",Tmur)
	myProblem.setDirichletBoundaryCondition("Radiateur_droite",Tmur)
	myProblem.setDirichletBoundaryCondition("Mur",Tmur)

    # name of result file
	fileName = "StationnaryDiffusion_3DVF_StructuredCubes";

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
    StationaryDiffusionEquation_3DVF_RoomCooling_StructuredCubes()
