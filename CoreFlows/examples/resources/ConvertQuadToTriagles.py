from medcoupling import *
from mechanic import *  
# Charger le fichier .med
mesh_file = "AnnulusSpiderWeb40x128.med"
mesh_name = "spider_web"  # Nom correct du maillage

# Charger le maillage depuis le fichier .med
mesh = ReadMeshFromFile(mesh_file, mesh_name, 0)  # 0 = premier pas de temps


mesh.QuadToTri([1,2], SMESH.FT_MinimumAngle)


















# Vérifier le type des éléments présents dans le maillage
cell_types = mesh.getAllGeoTypes()  # Récupère tous les types d'éléments

# Vérifier s'il y a des quadrilatères (type 4)
if 4 not in cell_types:
    raise ValueError("Le maillage ne contient pas de quadrilatères.")

# Extraire les cellules quadrilatères
quad_cells = mesh.getAllCellsOfType(4)  # Type 4 = quadrilatères

# Convertir les quadrilatères en triangles
tri_cells = quad_cells.convertToTriangles()

# Créer un nouveau maillage triangulaire
tri_mesh = mesh.deepCopy()
tri_mesh.allocateCells(len(tri_cells))

for i, tri in enumerate(tri_cells):
    tri_mesh.insertNextCell(tri)

tri_mesh.finishInsertingCells()

# Sauvegarder le maillage triangulaire
output_file = "triangular_mesh.med"
tri_mesh.writeToFile(output_file)
print(f" Maillage triangulaire sauvegardé sous {output_file}")