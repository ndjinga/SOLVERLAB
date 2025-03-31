import salome
from salome.geom import geomBuilder
import GEOM
import SMESH
from salome.smesh import smeshBuilder
import math

salome.salome_init()
theStudy = salome.myStudy
geompy = geomBuilder.New(theStudy)

# Créer les points du grand parallélogramme
P0 = geompy.MakeVertex(0, 0, 0)
P1 = geompy.MakeVertex(0.5, 1, 0)
P2 = geompy.MakeVertex(1.5, 1, 0)
P3 = geompy.MakeVertex(1, 0, 0)
P4 = geompy.MakeVertex(0.25, 0.5, 0)
P5 = geompy.MakeVertex(1.25, 0.5, 0)
""" P6 = geompy.MakeVertex(0.5, 0, 0)
P7 = geompy.MakeVertex(1, 1, 0)
P8 = geompy.MakeVertex(0.75, 0.5, 0) """

# Créer les arêtes du grand parallélogramme
L0 = geompy.MakeLineTwoPnt(P0, P1)
L1 = geompy.MakeLineTwoPnt(P1, P2)
L2 = geompy.MakeLineTwoPnt(P2, P3)
L3 = geompy.MakeLineTwoPnt(P3, P0)

# Créer les lignes internes pour la subdivision
L4 = geompy.MakeLineTwoPnt(P4, P5)
L5 = geompy.MakeLineTwoPnt(P0, P4)
L6 = geompy.MakeLineTwoPnt(P3, P5)
L7 = geompy.MakeLineTwoPnt(P5, P2)
L8 = geompy.MakeLineTwoPnt(P4, P1)
""" L9 = geompy.MakeLineTwoPnt(P6, P8)
L10 = geompy.MakeLineTwoPnt(P7, P8) """

# Créer des sous-parallélogrammes (contours fermés)
SubFace1 = geompy.MakeFaceWires([L5, L4, L6, L3], 1)  # Contour fermé : L0, L7, L4, L5
SubFace2 = geompy.MakeFaceWires([L7, L1, L8, L4], 1)  # Contour fermé : L7, L1, L8, L4
""" SubFace3 = geompy.MakeFaceWires([L5, L4, L6, L3], 1)  # Contour fermé : L5, L4, L6, L3
SubFace4 = geompy.MakeFaceWires([L4, L8, L2, L6], 1)  # Contour fermé : L4, L8, L2, L6 """

# Fusionner toutes les sous-faces en un seul composé
CompoundFaces = geompy.MakeCompound([SubFace1, SubFace2]) #

# Ajouter la géométrie à l'étude
geompy.addToStudy(P0, "P0")
geompy.addToStudy(P1, "P1")
geompy.addToStudy(P2, "P2")
geompy.addToStudy(P3, "P3")
geompy.addToStudy(P4, "P4")
geompy.addToStudy(P5, "P5")
geompy.addToStudy(L0, "L0")
geompy.addToStudy(L1, "L1")
geompy.addToStudy(L2, "L2")
geompy.addToStudy(L3, "L3")
geompy.addToStudy(L4, "L4")
geompy.addToStudy(L5, "L5")
geompy.addToStudy(L6, "L6")
geompy.addToStudy(L7, "L7")
geompy.addToStudy(L8, "L8")
geompy.addToStudy(SubFace1, "SubFace1")
geompy.addToStudy(SubFace2, "SubFace2")
""" geompy.addToStudy(SubFace3, "SubFace3")
geompy.addToStudy(SubFace4, "SubFace4") """
geompy.addToStudy(CompoundFaces, "Sous-Parallelogrammes")

# Créer le maillage
smesh = smeshBuilder.New()
Mesh = smesh.Mesh(CompoundFaces)

# Définir les hypothèses de maillage
algo = Mesh.Segment()
algo.NumberOfSegments(2)  # Divise chaque côté en 2 segments pour un maillage uniforme

# Créer des éléments quadrangulaires
algo2D = Mesh.Quadrangle()
#algo2D.Quadrangle()

# Générer le maillage
isDone = Mesh.Compute()

# Exporter le maillage pour vérification
smesh.SetName(Mesh, 'Parallelogram_Mesh')
smesh.SetName(algo, 'Segment_Hypothesis')

# Exportation du maillage en format .med
Mesh.ExportMED('./Parallelogram_Mesh.med')