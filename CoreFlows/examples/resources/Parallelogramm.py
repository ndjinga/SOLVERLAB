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
# Créer les arêtes du grand parallélogramme
L0 = geompy.MakeLineTwoPnt(P0, P1)
L1 = geompy.MakeLineTwoPnt(P1, P2)
L2 = geompy.MakeLineTwoPnt(P2, P3)
L3 = geompy.MakeLineTwoPnt(P3, P0)
# Créer de parallélogramme (contours fermés)
Face = geompy.MakeFaceWires([L0, L1, L2, L3], 1) 

# Ajouter la géométrie à l'étude
geompy.addToStudy(P0, "P0")
geompy.addToStudy(P1, "P1")
geompy.addToStudy(P2, "P2")
geompy.addToStudy(P3, "P3")

geompy.addToStudy(L0, "L0")
geompy.addToStudy(L1, "L1")
geompy.addToStudy(L2, "L2")
geompy.addToStudy(L3, "L3")
geompy.addToStudy(Face, "Parallelogramme")

# Créer le maillage
smesh = smeshBuilder.New()
Mesh = smesh.Mesh(Face)

# Définir les hypothèses de maillage
algo = Mesh.Segment()
algo.NumberOfSegments(40)  # Divise chaque côté en 2 segments pour un maillage uniforme

# Créer des éléments quadrangulaires
algo2D = Mesh.Quadrangle()
# Générer le maillage
isDone = Mesh.Compute()

# Exporter le maillage pour vérification
smesh.SetName(Mesh, 'Parallelogram_Mesh')
smesh.SetName(algo, 'Segment_Hypothesis')

# Exportation du maillage en format .med
Mesh.ExportMED('./Parallelogram_Mesh40.med')