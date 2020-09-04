# Copyright (C) 2012-2014  CEA/DEN, EDF R&D, OPEN CASCADE
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

# --- resources ---

# uic files 
SET(_pyuic_files
  MainCFWidget.ui
)

# --- scripts ---
SET(_all_lib_SCRIPTS
    MainCFWidget.py      
)

# --- rules ---
SALOME_INSTALL_SCRIPTS("${_pyuic_SCRIPTS}" ${SALOME_INSTALL_PYTHON})
SALOME_INSTALL_SCRIPTS("${_all_lib_SCRIPTS}" ${SALOME_INSTALL_PYTHON})

INSTALL(FILES ${_pyuic_files} DESTINATION ${SALOME_INSTALL_PYTHON})
