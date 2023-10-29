#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# Copyright (C) 2008-2023  CEA
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
# 
# See http://www.salome-platform.org or email : webmaster.salome@opencascade.com
# %% LICENSE_END

import sys

"""
***********************************************
* prerequisites environment variables 
* to use solverlabGUI
***********************************************

# example of set environment to use solverlabGUI

bash

export CONDADIR=/volatile/wambeke/miniconda3
export PATH=$CONDADIR/bin:$PATH
conda activate pyscientific # as correct python3 PyQt5
  
# useful environ var
SOLVERLABGUI_ROOT_DIR=.../solverlabGUI

# optionally where solverlab CODE .de is located
# default is $SOLVERLABGUI_ROOT_DIR/solverlabCode
SOLVERLABCODE_ROOT_DIR=YourChoiceForAnotherDevelopmentDirectory


${SOLVERLABGUI_ROOT_DIR}/solverlabGUI --help    # help
${SOLVERLABGUI_ROOT_DIR}/solverlabGUI --gui     # launch GUI
"""



#####################################
_importErrorMessage = r"""
***********************************************
* solverlabGUI configuration is incorrect.
***********************************************

in 'classical solverlabGUI installation' configuration.
user usually have minimum environ variable set.

  >> bash
  >> SOLVERLABGUI_ROOT_DIR=...    # where solverlabGUI is located
  >> ${SOLVERLABGUI_ROOT_DIR}/solverlabGUI --gui  # launch GUI

Usual python configuration is (2020):
  python 3.7
  PyQt5, PyQt5.QtWebEngine
  xml
  numpy
  matplotlib
  pandas

Tricks:
  Easily get local python3 configuration (as no superuser Linux or Windows) 
  
  see:
  ${SOLVERLABGUI_ROOT_DIR}/sandbox/conda/README_pySolverlabGUI3.md
  
  and create correct python environment for solverlabGUI
 
"""

#####################################
_rootFeaturesErrorMessage = r"""
***********************************************
* solverlab CODE configuration is incorrect.
***********************************************

in 'classical solverlab CODE compilation/installation' configuration.
user usually have minimum environ variable set.

  >>> export SOLVERLABCODE_ROOT_DIR=...       # where solverlab .de CODE is located

In solverlabGUI defaut value for SOLVERLABCODE_ROOT_DIR is ${SOLVERLABGUI_ROOT_DIR}/solverlabCode
"""

import subprocess as SP

def error(message):
  print("ERROR: %s" % message)
  return

def TestImports():
  """
  test of prerequisites import python for SolverlabGui
  message for problem(s), aborting immediatly. 
  """
  # TODO imports = "PyQt5 PyQt5.QtWebEngineWidgets salomepy xyzpy numpy pandas".split()
  imports = "PyQt5 salomepy xyzpy zlib numpy pandas".split()

  res = "OK"
  for ii in imports:
    try:
      exec("import %s" % ii)
    except:
      error("problem with 'import %s'" % ii)
      res = "KO"

  if res == "KO":
    print(_importErrorMessage)
    # sys.exit("------> Fix it. Aborting.\n")
    # sys.exit("------> Fix it. Now it's at your owns risks.\n")
  return res

def TestSolverlabFeatures():
  """
  test of solverlab CODE features
  """
  cmd = "solverlab.exe -h"
  res = SP.Popen(cmd, shell=True, stdout=SP.PIPE, stderr=SP.PIPE).communicate()
  if res[1] != "":
    error("problem with '%s': %s" % (cmd, res[1])) #stderr
  stdout = res[0].split()
  prereq = "qt python".split()
  res = "OK"
  for ii in prereq:
    if ii not in stdout:
      error("problem with necessary ROOT compilation feature '%s'" % ii)
      res = "KO"

  if res == "KO":
    print(_rootFeaturesErrorMessage)
    print(_importErrorMessage)
    #print("------> Fix it. Now it's at your owns risks.\n")
    sys.exit("------> Fix it. Aborting.\n")
  return res


# TestImports()
# TestSolverlabFeatures()

