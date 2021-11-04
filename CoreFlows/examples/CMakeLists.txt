project(test)

file(GLOB MESHES_TO_INSTALL resources )
install(DIRECTORY ${MESHES_TO_INSTALL} DESTINATION share/examples)

file(COPY resources DESTINATION ${CMAKE_CURRENT_BINARY_DIR} 
    FILE_PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE
)


add_subdirectory (C)

if (COREFLOWS_WITH_PYTHON )

  add_subdirectory (Python)

endif (COREFLOWS_WITH_PYTHON )

add_custom_target (CoreFlows COMMAND ctest -O testsCoreFlows.log)

add_custom_target (cpp COMMAND ctest -R .exe)# may be replace ctest -R with ctest -L
add_custom_target (mpi COMMAND ctest -R Procs)# may be replace ctest -R with ctest -L
add_custom_target (eos COMMAND ctest -R EOS)# may be replace ctest -R with ctest -L
add_custom_target (coupled COMMAND ctest -R Coupled)# may be replace ctest -R with ctest -L

add_custom_target (SinglePhase COMMAND ctest -R SinglePhase)# may be replace ctest -R with ctest -L
add_custom_target (DriftModel  COMMAND ctest -R DriftModel) # may be replace ctest -R with ctest -L
add_custom_target (IsothermalTwoFluid COMMAND ctest -R IsothermalTwoFluid)# may be replace ctest -R with ctest -L
add_custom_target (   FiveEqsTwoFluid COMMAND ctest -R    FiveEqsTwoFluid)# may be replace ctest -R with ctest -L

add_custom_target (DiffusionEquation COMMAND ctest -R DiffusionEquation -E StationaryDiffusionEquation)# may be replace ctest -R with ctest -L
add_custom_target (diffusion_tests         COMMAND ctest -R DiffusionEquation)# may be replace ctest -R with ctest -L
add_custom_target (TransportEquation COMMAND ctest -R TransportEquation)# may be replace ctest -R with ctest -L
add_custom_target (transport_tests         COMMAND ctest -R TransportEquation)# may be replace ctest -R with ctest -L
add_custom_target (StationaryDiffusionEquation COMMAND ctest -R StationaryDiffusionEquation)# may be replace ctest -R with ctest -L

add_custom_target (convergence COMMAND ctest -R convergence)# may be replace ctest -R with ctest -L

add_custom_target (fv_tests         COMMAND ctest -R FV)# may be replace ctest -R with ctest -L
add_custom_target (FV_tests         COMMAND ctest -R FV)# may be replace ctest -R with ctest -L
add_custom_target (fe_tests         COMMAND ctest -R FE)# may be replace ctest -R with ctest -L
add_custom_target (FE_tests         COMMAND ctest -R FE)# may be replace ctest -R with ctest -L
add_custom_target (1D_tests         COMMAND ctest -R 1D)# may be replace ctest -R with ctest -L
add_custom_target (2D_tests         COMMAND ctest -R 2D)# may be replace ctest -R with ctest -L
add_custom_target (3D_tests         COMMAND ctest -R 3D)# may be replace ctest -R with ctest -L
add_custom_target (Dirichlet  COMMAND ctest -R Dirichlet)# may be replace ctest -R with ctest -L
add_custom_target (Neumann    COMMAND ctest -R Neumann)# may be replace ctest -R with ctest -L

