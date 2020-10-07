
if (CDMATH_WITH_PYTHON AND CDMATH_WITH_POSTPRO)#The environment is probably too rich for the current test contain
  ADD_TEST(Test_Boys_Surface ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/BoysSurfaceVTK.py)
endif ()

INSTALL(FILES BoysSurfaceVTK.py VTKnewReader.py DESTINATION bin/BoysSurface)

