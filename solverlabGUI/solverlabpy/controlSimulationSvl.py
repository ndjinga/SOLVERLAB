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


import os
import sys
import subprocess as SP
import platform

import xyzpy.classFactoryXyz as CLFX
import xyzpy.loggingXyz as LOG
logger = LOG.getLogger()

#set classes Xyz in factory xyzpy.utilsXyz
from xyzpy.baseXyz import _XyzConstrainBase
import xyzpy.intFloatListXyz as IFLX

import solverlabpy.utilsSvl as UTS
import debogpy.debug as DBG

import solverlabpy.configSvl as CFGSVL


verbose = False


###############################################################
class UserModeSvl(IFLX.StrInListXyz):

  _allowedList = CFGSVL.getExistingModes()
  _defaultValue = _allowedList[0] # static
  pass


########################################################################################
class ControlSimulationSvl(_XyzConstrainBase):
  """
  general informations about
  """
  _attributesList = [ #list, not a dict because sequential order list is used in files
    ("hostName", "StrNoEditionXyz"),
    ("system", "StrNoEditionXyz"),
    ("platform", "StrNoEditionXyz"),
    # TODO synchronize with ModelSvl.userMode(self)  ("userMode", "UserModeSvl"),
    # TODO("settings", "SettingsXyz"),
  ]
  _icon = "controlnmd"
  _helpDict = {
    "hostName": (u"current hostname", ""),
    "system": (u"current system", ""),
    "platform": (u"Operating System details", ""),
    "userMode": (u"user mode of solverlab input data display", ""),
    #TODO? "settings": (u"Settings from uranieGui.ini", ""),
  }
  
  def __init__(self):
    super(ControlSimulationSvl, self).__init__()
    self._defautNameAsRoot = "Controls"
    self.setIsCast(True)
    self._setAllAttributesList()

  def setDefaultValues(self):
    super(ControlSimulationSvl, self).setDefaultValues()
    
    self.hostName = platform.node()
    self.system = platform.system()
    self.platform = platform.platform()
    # DBG.write("dir(userMode)", dir(self.userMode), True)

    """
    # no uranie, no root, may be later...
    cmd = "root-config --version"
    res = UTS.Popen(cmd, logger=logger)
    DBG.write("root-config --version", res, True)
    try:
      stdout = res.getValue()
      tmp = stdout.split('.')[1] #exception if not '.'
      self.rootVersion = stdout.strip() 
    except:
      self.rootVersion == "unknown"
         
    cmd = 'grep "#define URANIE_RELEASE " ${URANIESYS}/include/UVersion.h'
    res = UTS.Popen(cmd, logger=logger)
    DBG.write("URANIE_RELEASE", res, True)
    try:
      stdout = res.getValue()
      tmp = stdout.split('.')[1] #exception if not '.'
      self.uranieVersion = stdout.split('"')[1] #exception if not '"'
    except:
      self.uranieVersion == "unknown"
    """


#factory pattern using xyzpy.utilsXyz._dictOfXyzClass
CLFX.appendAllXyzClasses( [UserModeSvl, ControlSimulationSvl] )


