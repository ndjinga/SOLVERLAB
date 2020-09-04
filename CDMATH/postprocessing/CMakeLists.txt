cmake_minimum_required (VERSION 2.8)
PROJECT(postprocessing)

file(COPY analyticalFunction_discretizer.py VTK_routines.py PV_routines.py MED_routines.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

INSTALL(FILES analyticalFunction_discretizer.py VTK_routines.py PV_routines.py MED_routines.py DESTINATION bin/cdmath/postprocessing)

