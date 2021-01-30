# Copyright (C) 2012-2021  CEA/DEN, EDF R&D, OPEN CASCADE
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#



SET(GUI_bin_SCRIPTS
    ${CMAKE_CURRENT_SOURCE_DIR}/CoreFlows_Standalone.py
)

SET(GUI_lib_SCRIPTS
    ${CMAKE_CURRENT_SOURCE_DIR}/CFDesktop.py
)


IF   ( (EXISTS $ENV{CONFIGURATION_ROOT_DIR}) AND (EXISTS $ENV{KERNEL_ROOT_DIR}) AND (EXISTS $ENV{GUI_ROOT_DIR}) )
  ADD_SUBDIRECTORY(salome)
ELSE ( (EXISTS $ENV{CONFIGURATION_ROOT_DIR}) AND (EXISTS $ENV{KERNEL_ROOT_DIR}) AND (EXISTS $ENV{GUI_ROOT_DIR}) )
  INSTALL(FILES ${GUI_bin_SCRIPTS} DESTINATION bin/salome)
  INSTALL(FILES ${GUI_lib_SCRIPTS} DESTINATION lib/python)

  # Configure file SalomePyQt_MockUp.py
  SET(SGPYQT_RES_DIR "share/salome/resources")
  configure_file(
      SalomePyQt_MockUp.py.in
      SalomePyQt_MockUp.py
      @ONLY
  )
  INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/SalomePyQt_MockUp.py DESTINATION lib/python)
  # Configure file utils.py
  SET(SALOME_INSTALL_PYTHON "lib/python")#This is to use the same variable as salome
  configure_file(
      utils.py.in
      utils.py
      @ONLY
  )
  INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/utils.py DESTINATION lib/python)
ENDIF( (EXISTS $ENV{CONFIGURATION_ROOT_DIR}) AND (EXISTS $ENV{KERNEL_ROOT_DIR}) AND (EXISTS $ENV{GUI_ROOT_DIR}) )

ADD_SUBDIRECTORY(ui)

