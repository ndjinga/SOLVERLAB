import MEDLoader as ml
import os, sys

def read_typ2(fichier, nom_med):
    with open(fichier, "r") as fic: lines = fic.readlines()
    nb_som = int(lines[1].strip())
    nb_cel = int(lines[3 + nb_som].strip())

    som = [map(float, line.split()) for line in lines[2:2 + nb_som]]

    # les connectivites ne partent pas de 0...
    tmp_cel = [map(int, line.split()[1:]) for line in lines[4 + nb_som:]]
    cel = [[p - 1 for p in c] for c in tmp_cel]

    mesh = ml.MEDCouplingUMesh("mesh", 2)
    mesh.allocateCells(len(cel))
    for p in cel: mesh.insertNextCell(ml.NORM_POLYGON, len(p), p)
    mesh.finishInsertingCells()

    pts = []
    for p in som: pts.extend(p)
    co = ml.DataArrayDouble(pts, len(pts) / 2, 2)
    mesh.setCoords(co)

    mf, d, di, r, ri = mesh.buildDescendingConnectivity()
    mf.setName("mesh")
    mm = ml.MEDFileUMesh.New()
    mm.setMeshAtLevel(0, mesh)
    mm.setMeshAtLevel(-1, mf)

    g = []
    nb_vois = ri.deltaShiftIndex()
    for i in range(mf.getNumberOfCells()):
        if nb_vois[i] == 1: g.append(i)
    grp = ml.DataArrayIdType.New(g)
    grp.setName("boundary")
    mm.addGroup(-1, grp)

    mm.write("{}/mesh.med".format(nom_med), 2)

def read_between(infile, patern1, patern2):
    with open(infile) as fic:
        copy = False
        string = ''
        for line in fic:
            if any([line.strip().lower().startswith(patern.lower()) for patern in patern1]):
                copy = True
            elif any([line.strip().lower().startswith(patern.lower()) for patern in patern2]):
                copy = False
            elif copy:
                string += line
    return string

def connectivity_from_string(string):
    indicies = map(int, string.split())
    tab = []
    off = 0
    while off < len(indicies):
        off_old = off + 1
        off += indicies[off] + 1
        tab.append(indicies[off_old:off])
    return tab

def read_typ3(fichier, nom_med):
    with open(fichier, "r") as fic: lines = fic.readlines()
    nb_som = int(lines[9].strip())
    nb_cel = int(lines[11].strip())
    nb_fac = int(lines[13].strip())

    som = [map(float, line.split()) for line in lines[17:17 + nb_som]]

    # les connectivites ne partent pas de 0...
    s = read_between(fichier, ["Volumes->Faces"], ["Volumes->Vertices"])
    tmp_cel = connectivity_from_string(s)
    # tmp_cel = [map(int, line.split()[1:]) for line in lines[18 + nb_som:18 + nb_som + nb_cel]]
    cel = [[p - 1 for p in c] for c in tmp_cel]

    s = read_between(fichier, ["Faces->Vertices"], ["Faces->Control volumes", "Faces->volumes"])
    tmp_fac = connectivity_from_string(s)
    # tmp_fac = [map(int, line.split()[1:]) for line in lines[21 + nb_som + 2 * nb_cel + nb_fac:21 + nb_som + 2 * nb_cel + 2 * nb_fac]]
    fac = [[p - 1 for p in c] for c in tmp_fac]

    mesh = ml.MEDCouplingUMesh("mesh", 3)
    mesh.allocateCells(len(cel))

    for e in range(len(cel)):
        con = []
        for face in cel[e]:
            con.extend(fac[face])
            con.append(-1)
        mesh.insertNextCell(ml.NORM_POLYHED, con[:-1])
    mesh.finishInsertingCells()

    pts = []
    for p in som: pts.extend(p)
    co = ml.DataArrayDouble(pts, len(pts) / 3, 3)
    mesh.setCoords(co)

    mf, d, di, r, ri = mesh.buildDescendingConnectivity()
    mf.setName("mesh")
    mm = ml.MEDFileUMesh.New()
    mm.setMeshAtLevel(0, mesh)
    mm.setMeshAtLevel(-1, mf)

    g = []
    nb_vois = ri.deltaShiftIndex()
    for i in range(mf.getNumberOfCells()):
        if nb_vois[i] == 1: g.append(i)
    grp = ml.DataArrayIdType.New(g)
    grp.setName("boundary")
    mm.addGroup(-1, grp)

    mm.write("{}/mesh.med".format(nom_med), 2)

# maillages 3D du benchmark FVCA6
# meshes = (("meshAA-random",      "RandMesh",     ("4", "8", "16", "32")),
          # #("meshBB_well",        "WellMesh_",     ("1", "2", "3", "4", "5", "6", "7")),
          # #("meshB_tetra",        "tet.",          ("00", "0", "1", "2", "3", "4", "5", "6")),
          # ("meshC_voro",         "vmesh_",        ("1", "2", "3", "4", "5")),
          # ("meshD_kershaw",      "dkershaw",      ("08", "16", "32", "64")),
          # ("meshF_dbls",         "dbls_",         ("10", "20", "30", "40")))
          # # ("meshH_locrafgrid",   "locrafgrid_",   ("1", "2", "3", "4", "5")),
          # # ("meshI_checkerboard", "checkerboard_", ("2x2x2", "4x4x4", "8x8x8", "16x16x16", "32x32x32")))
# for t, m, d in meshes:
    # for n in d:
        # print( t, n)
        # folder = "{}/jdd_{}".format(t, n)
        # os.system("mkdir -p {}".format(folder))
        # read_typ3("Meshes_3D/{}/{}{}.msh".format(t, m, n), folder)

if __name__ == "__main__":

	if len(sys.argv) != 2:
	  print("USAGE: convert_gmsh_to_med.py file.typ")
	  sys.exit(-1)
	
	filename = sys.argv[1]
	print("Converting ", filename)

	l=len(filename)
	name=filename[:l-5]
	extension=filename[l-5:]
	if extension==".typ2":
		read_typ2(filename, name+".med")
	elif extension==".typ3":
		read_typ3(filename, name+".med")
	else :
		raise ValueError("File "+filename+" has unknown file extension "+extension)
