# -*-coding:utf-8 -*
#===============================================================================================================================
# Name        : Résolution EF de l'équation de Laplace 3D -\Delta T = 0 avec conditions aux limites de Dirichlet u non nulle
# Authors     : Michaël Ndjinga, Sédrick Kameni Ngwamou
# Copyright   : CEA Saclay 2019
# Description : Utilisation de la méthode des éléménts finis P1 avec champs u discrétisés aux noeuds d'un maillage tétraédrique
#               Condition limites correspondant au refroidissement dû à une fenêtre et au chauffage dû à un radiateur
#				Création et sauvegarde du champ résultant ainsi que du champ second membre en utilisant la librairie CDMATH
#================================================================================================================================

import cdmath
import numpy as np

# Fonction qui remplace successivement les colonnes d'une matrices par un vecteur donné et retourne la liste des déterminants
def gradientNodal(M, values):
	matrices=[0]*(len(values)-1)
	for i in range(len(values)-1):
		matrices[i] = M.deepCopy()        
		for j in range(len(values)):
			matrices[i][j,i] = values[j]

	result = cdmath.Vector(len(values)-1)    
	for i in range(len(values)-1):
		result[i] = matrices[i].determinant()

	return result

# Fonction qui calcule la valeur de la condition limite en un noeud donné
def boundaryValue(nodeId): 
	Ni=my_mesh.getNode(nodeId)

	# 4 groupes sont considérés sur le bord
	if boundaryNodes.count(nodeId)==0:
		return 0
	elif Ni.getGroupName()=='Fenetre':
		return Tfenetre;
	elif Ni.getGroupName()=='Radiateur_sous-fenetre':
		return Tradiateur
	#elif Ni.getGroupName()=='Radiateur_Devant':
		#return Tradiateur;
	#elif Ni.getGroupName()=='Radiateur_droit':
		#return Tradiateur
	else: 	
		return Tmur;

#Chargement du maillage tétraédrique du domaine
#==============================================
my_mesh = cdmath.Mesh("Mesh_RadiatorAndWindow.med")
if( my_mesh.getSpaceDimension()!=3 or my_mesh.getMeshDimension()!=3) :
    raise ValueError("Wrong space or mesh dimension : space and mesh dimensions should be 3")
if(not my_mesh.isTetrahedral()) :
	raise ValueError("Wrong cell types : mesh is not made of tetrahedra")

nbNodes = my_mesh.getNumberOfNodes()
nbCells = my_mesh.getNumberOfCells()

print("Mesh loading done")
print("nb of nodes=", nbNodes)
print("nb of cells=", nbCells)

#Conditions limites
Tmur=20
Tfenetre=0
Tradiateur=40

#Détermination des noeuds intérieurs
#======================================================================
nbInteriorNodes = 0
nbBoundaryNodes = 0
maxNbNeighbours = 0#This is to determine the number of non zero coefficients in the sparse finite element rigidity matrix
interiorNodes = []
boundaryNodes = []

#parcours des noeuds pour discrétisation du second membre et extraction 1) des noeuds intérieur 2) des noeuds frontière 3) du nb max voisins d'un noeud
for i in range(nbNodes):
	Ni=my_mesh.getNode(i)
    
	if my_mesh.isBorderNode(i): # Détection des noeuds frontière getGroupName my_mesh.getNode(i)
		boundaryNodes.append(i)
	else: # Détection des noeuds intérieurs
		interiorNodes.append(i)
		maxNbNeighbours= max(1+Ni.getNumberOfEdges(),maxNbNeighbours) 

nbInteriorNodes=len(interiorNodes)
nbBoundaryNodes=len(boundaryNodes)


print("nb of interior nodes=", nbInteriorNodes)
print("nb of Boundary nodes=", nbBoundaryNodes)
print("Max nb of neighbours=", maxNbNeighbours)


# Construction de la matrice de rigidité et du vecteur second membre du système linéaire
#=======================================================================================
Rigidite=cdmath.SparseMatrixPetsc(nbInteriorNodes,nbInteriorNodes,maxNbNeighbours)
RHS=cdmath.Vector(nbInteriorNodes)

# Vecteurs gradient de la fonction de forme associée à chaque noeud d'un tétrèdre (hypothèse 3D)
M=cdmath.Matrix(4,4)
GradShapeFunc0=cdmath.Vector(3)
GradShapeFunc1=cdmath.Vector(3)
GradShapeFunc2=cdmath.Vector(3)
GradShapeFunc3=cdmath.Vector(3)

#On parcourt les tétraèdres du domaine pour remplir la matrice
for i in range(nbCells):

	Ci=my_mesh.getCell(i)

	#Extraction des noeuds de la cellule
	nodeId0=Ci.getNodeId(0)
	nodeId1=Ci.getNodeId(1)
	nodeId2=Ci.getNodeId(2)
	nodeId3=Ci.getNodeId(3)
	N0=my_mesh.getNode(nodeId0)
	N1=my_mesh.getNode(nodeId1)
	N2=my_mesh.getNode(nodeId2)
	N3=my_mesh.getNode(nodeId3)

	M[0,0]=N0.x()
	M[0,1]=N0.y()
	M[0,2]=N0.z()
	M[0,3]=1
	M[1,0]=N1.x()
	M[1,1]=N1.y()
	M[1,2]=N1.z()
	M[1,3]=1
	M[2,0]=N2.x()
	M[2,1]=N2.y()
	M[2,2]=N2.z()
	M[2,3]=1
	M[3,0]=N3.x()
	M[3,1]=N3.y()
	M[3,2]=N3.z()
	M[3,3]=1

	#Values of each shape function at each node
	values0=[1,0,0,0]
	values1=[0,1,0,0]
	values2=[0,0,1,0]
	values3=[0,0,0,1]

	GradShapeFunc0 = gradientNodal(M,values0)/6
	GradShapeFunc1 = gradientNodal(M,values1)/6
	GradShapeFunc2 = gradientNodal(M,values2)/6
	GradShapeFunc3 = gradientNodal(M,values3)/6
	
	#Création d'un tableau (numéro du noeud, gradient de la fonction de forme)
	GradShapeFuncs={nodeId0 : GradShapeFunc0}
	GradShapeFuncs[nodeId1]=GradShapeFunc1
	GradShapeFuncs[nodeId2]=GradShapeFunc2
	GradShapeFuncs[nodeId3]=GradShapeFunc3

	# Remplissage de  la matrice de rigidité et du second membre
	for j in [nodeId0,nodeId1,nodeId2,nodeId3] : 
		if boundaryNodes.count(j)==0: #seuls les noeuds intérieurs contribuent au système linéaire (matrice de rigidité et second membre)
			j_int=interiorNodes.index(j)#indice du noeud j en tant que noeud intérieur
			boundaryContributionAdded=False#Needed in case j is a border cell
			#Ajout de la contribution de la cellule ttétraédrique i au second membre du noeud j 
			for k in [nodeId0,nodeId1,nodeId2,nodeId3] : 
				if boundaryNodes.count(k)==0 : #seuls les noeuds intérieurs contribuent à la matrice du système linéaire
					k_int = interiorNodes.index(k)#indice du noeud k en tant que noeud intérieur
					Rigidite.addValue(j_int,k_int,GradShapeFuncs[j]*GradShapeFuncs[k]/Ci.getMeasure())
				elif boundaryContributionAdded == False: #si condition limite non nulle au bord (ou maillage non recouvrant), ajouter la contribution du bord au second membre de la cellule j
					# Valeurs de g_h aux noeuds du tétraèdre
					T0 = boundaryValue(nodeId0)
					T1 = boundaryValue(nodeId1)
					T2 = boundaryValue(nodeId2)
					T3 = boundaryValue(nodeId3)
					boundaryContributionAdded=True#Contribution from the boundary to matrix line j is done in one step
					GradGh = gradientNodal(M,[T0,T1,T2,T3])/6
					RHS[j_int] += -(GradGh*GradShapeFuncs[j])/Ci.getMeasure()
            

    
print("Linear system matrix building done")

LS=cdmath.LinearSolver(Rigidite,RHS,100,1.E-6,"CG","ILU")#,"ILU" Remplacer CG par CHOLESKY pour solveur direct

SolSyst=LS.solve()

# Création du champ résultat
#===========================
my_Temperature = cdmath.Field("Temperature", cdmath.NODES, my_mesh, 1)
for j in range(nbInteriorNodes):
   	my_Temperature[interiorNodes[j]]=SolSyst[j];#remplissage des valeurs pour les noeuds intérieurs SolSyst[j]
#Remplissage des valeurs pour les noeuds frontière (condition limite)
for j in boundaryNodes:
    my_Temperature[j]=boundaryValue(j)


#sauvegarde sur le disque dur du résultat dans un fichier paraview
my_Temperature.writeVTK("FiniteElements3DTemperature")

print( "Minimum temperature= ", my_Temperature.min(), ", maximum temperature= ", my_Temperature.max() )
assert my_Temperature.min()>= min(Tmur,Tfenetre,Tradiateur) and my_Temperature.max()<= max(Tmur,Tfenetre,Tradiateur)

print( "Numerical solution of 3D Laplace equation using finite elements done" )

