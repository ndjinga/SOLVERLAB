file(GLOB WaveSystem_EXAMPLES_TO_INSTALL 

  WaveSystem1DUpwind WaveSystem1DUpwind_RiemannProblem WaveSystem1Staggered_RiemannProblem WaveSystem2DUpwind_RiemannProblem
  WaveSystem_Stationary/WaveSystemUpwind WaveSystem_Stationary/WaveSystemCentered WaveSystem_Stationary/WaveSystemPStag WaveSystem_Stationary/WaveSystemStaggered 
  WaveSystem_Shock/WaveSystemUpwind WaveSystem_Shock/WaveSystemCentered WaveSystem_Shock/WaveSystemPStag WaveSystem_Shock/WaveSystemStaggered 
)

install(DIRECTORY ${WaveSystem_EXAMPLES_TO_INSTALL} DESTINATION share/examples/WaveSystem)

SET(MESH_MED
  ${MED_MESHES}/squareWithTriangles.med
  ${MED_MESHES}/squareWithSquares.med
  ${MED_MESHES}/squareWithBrickWall.med
  ${MED_MESHES}/squareWithCheckerboardSquares.med
  ${MED_MESHES}/squareWithDeformedQuadrangles.med
  ${MED_MESHES}/squareWithHexagons.med
  )

IF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)

    ADD_SUBDIRECTORY(WaveSystem1DUpwind)
    ADD_SUBDIRECTORY(WaveSystem1DStaggered_RiemannProblem)
    ADD_SUBDIRECTORY(WaveSystem1DUpwind_RiemannProblem)
    ADD_SUBDIRECTORY(WaveSystem2DUpwind_RiemannProblem)
    ADD_SUBDIRECTORY(WaveSystem_Stationary/WaveSystemUpwind)
    ADD_SUBDIRECTORY(WaveSystem_Stationary/WaveSystemPStag)
    ADD_SUBDIRECTORY(WaveSystem_Stationary/WaveSystemStaggered)
    ADD_SUBDIRECTORY(WaveSystem_Stationary/WaveSystemCentered)
    ADD_SUBDIRECTORY(WaveSystem_Shock/WaveSystemUpwind)
    ADD_SUBDIRECTORY(WaveSystem_Shock/WaveSystemPStag)
    ADD_SUBDIRECTORY(WaveSystem_Shock/WaveSystemStaggered)
    ADD_SUBDIRECTORY(WaveSystem_Shock/WaveSystemCentered)

ENDIF (CDMATH_WITH_PYTHON AND CDMATH_WITH_PETSC AND CDMATH_WITH_POSTPRO)


