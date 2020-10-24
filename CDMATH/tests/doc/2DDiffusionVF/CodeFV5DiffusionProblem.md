## The script

```python
#Discrétisation du second membre et extraction du nb max de voisins d'une cellule
#================================================================================
my_RHSfield = cdmath.Field("RHS_field", cdmath.CELLS, my_mesh, 1)
maxNbNeighbours=0#This is to determine the number of non zero coefficients in the sparse finite element rigidity matrix

for i in range(nbCells): 
	Ci = my_mesh.getCell(i)
	x = Ci.x()
	y = Ci.y()

	my_RHSfield[i]=2*pi*pi*sin(pi*x)*sin(pi*y)#mettre la fonction definie au second membre de l edp
	# compute maximum number of neighbours
	maxNbNeighbours= max(1+Ci.getNumberOfFaces(),maxNbNeighbours)

# Construction de la matrice et du vecteur second membre du système linéaire
#===========================================================================
Rigidite=cdmath.SparseMatrixPetsc(nbCells,nbCells,maxNbNeighbours)# warning : third argument is max number of non zero coefficients per line of the matrix
RHS=cdmath.Vector(nbCells)
normal=cdmath.Vector(dim)
#Parcours des cellules du domaine
for i in range(nbCells):
	RHS[i]=my_RHSfield[i] #la valeur moyenne du second membre f dans la cellule i
	Ci=my_mesh.getCell(i)
	for j in range(Ci.getNumberOfFaces()):# parcours des faces voisinnes
		Fj=my_mesh.getFace(Ci.getFaceId(j))
            		for idim in range(dim) :
                		normal[idim] = Ci.getNormalVector(j, idim);#normale sortante
		if not Fj.isBorder():
			k=Fj.getCellId(0)
			if k==i :
				k=Fj.getCellId(1)
			Ck=my_mesh.getCell(k)
			distance=Ci.getBarryCenter().distance(Ck.getBarryCenter())
			coeff=Fj.getMeasure()/Ci.getMeasure()/distance*(normal[0]*normal[0] + K*normal[1]*normal[1])
			Rigidite.addValue(i,k,-coeff) # terme extradiagonal
		else:
			coeff=Fj.getMeasure()/Ci.getMeasure()/Ci.getBarryCenter().distance(Fj.getBarryCenter())*(normal[0]*normal[0] + K*normal[1]*normal[1])
			#For the particular case where the mesh boundary does not coincide with the domain boundary
			x=Fj.getBarryCenter().x()
			y=Fj.getBarryCenter().y()
			RHS[i]+=coeff*sin(pi*x)*sin(pi*y)#mettre ici la condition limite du problème de Dirichlet
		Rigidite.addValue(i,i,coeff) # terme diagonal


# Résolution du système linéaire
#=================================
LS=cdmath.LinearSolver(Rigidite,RHS,500,1.E-6,"GMRES","ILU")
SolSyst=LS.solve()

# Automatic postprocessing :  save 2D picture and plot diagonal data
#===========================
PV_routines.Save_PV_data_to_picture_file("my_ResultField_0.vtu',"ResultField",'CELLS',"my_ResultField")
diag_data=VTK_routines.Extract_field_data_over_line_to_numpyArray(my_ResultField,[0,1,0],[1,0,0], resolution)
plt.plot(curv_abs, diag_data, label= str(nbCells)+ ' cells mesh')
plt.savefig("FV5_on_square_PlotOverDiagonalLine.png")

```
