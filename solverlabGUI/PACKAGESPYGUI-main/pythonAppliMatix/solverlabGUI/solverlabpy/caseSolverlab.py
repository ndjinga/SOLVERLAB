#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# Copyright (C) 2008-2015  CEA/DEN
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


import os
import platform
import math
from PyQt5 import QtGui, QtCore, QtWidgets as QTW
import subprocess as SP
import pprint as PP

from xyzpy.baseXyz import _XyzConstrainBase, ListOfBaseXyz

import xyzpy.utilsXyz as UXYZ
import xyzpy.classFactoryXyz as CLFX
import xyzpy.loggingXyz as LOG
import widgetpy.messageBoxDialog as MBD

logger = LOG.getLogger()

# set classes Xyz in factory xyzpy.utilsXyz
import xyzpy.intFloatListXyz as IFLX

import solverlabpy.abcdExpression as ABCD
import solverlabpy.configSvl as CFGSVL
import debogpy.debug as DBG

_Mole = 6.02e23

verbose = False
verboseEvent = False

###############################################################
# common methods
###############################################################

# ...

###############################################################
# classes
###############################################################

###############################################################
class BoolFalseSvl(IFLX.BoolXyz):
  """
  ['False', 'True'] with write config file solverlab strCfg [0, 1]
  default False
  """
  _defaultValue = "False"
  _toCfg = {'False': 0, 'True': 1}

  def strCfg(self):
    res = self._toCfg[str(self)]
    return res

###############################################################
class BoolTrueSvl(BoolFalseSvl):
  """
  ['False', 'True'] with write config file solverlab strCfg [0, 1]
  default True
  """
  _defaultValue = "True"
  pass


class SelectEquation(IFLX.StrNoEditionXyz):

  def createEditor(self, parent):
    equationList = self.getRoot().getListEquation()
    if len(equationList) == 0:
      QTW.QMessageBox.warning(self.getController().getDesktop(), "warning", "No equation")
      return None
    combo = IFLX.XyzQComboBox(parent)
    intlist = []
    for i in range(len(equationList)):
      intlist.append(str(i))
    combo.addItems(intlist)
    combo.setCurrentIndex(0)
    return combo
###############################################################


class CaseSvl(_XyzConstrainBase):
  """
  general informations about case and for launch solverlab
  """
  _attributesList = [  # list, not a dict because sequential order list is used in files
    ("launchInBackground", "BoolXyz"),
    ("Equation", "SelectEquation"),
    ("NumberOfProcessors", "IntPosXyz"),
 ]
  _icon = "casesvl"

  _helpDict = {
    "launchInBackground": (u"launch solverlab calculus in background", u""),
    "Equation": ("Used equation for calculus"),
    "NumberOfProcessors": ("Number of processors for MPI/Multithread run calculus"),
    # "Target": (u"Define target materials and geometry parameters", u""),
    # "Simulation": (u"Define general simulation parameters", u""),
    # "Variables": (u"""(right-click to add/modify a variable)
# Define simulation parameters as python coding variables""", u""),
  }

  _defaultVersion = "1.0.0"

  def __init__(self):
    super(CaseSvl, self).__init__()
    self.setIsCast(True)
    self._setAllAttributesList()

  def setDefaultValues(self):
    self.NumberOfProcessors = 1
    return
    '''self.IonBeam.setDefaultValues()
    self.Simulation.setDefaultValues()
    self.Target.setDefaultValues()'''

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    return CFGSVL.isHidden(self, nameAttr)



###############################################################
class SimulationSvl(_XyzConstrainBase):
  """
  general informations about simulation solverlab
  """
  _attributesList = [  # list, not a dict because sequential order list is used in files
    ("launch_in_background", "BoolFalseSvl"),
  ]
  _icon = "simulationsvl"


  _helpDict = {
    "launch_in_background": (u"Solverlab calculus in background", u""),
  }

  def __init__(self):
    super(SimulationSvl, self).__init__()
    self.setIsCast(True)
    self._setAllAttributesList()

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""

    res = False
    return res


# factory pattern using xyzpy.utilsXyz._dictOfXyzClass

CLFX.appendAllXyzClasses([
  SelectEquation,
  BoolFalseSvl, BoolTrueSvl,
  SimulationSvl,
  CaseSvl,
])
