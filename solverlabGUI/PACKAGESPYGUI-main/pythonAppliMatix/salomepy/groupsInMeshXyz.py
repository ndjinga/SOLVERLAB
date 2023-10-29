#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os

from xyzpy.baseXyz import _XyzConstrainBase, BaseFreeXyz
from xyzpy.intFloatListXyz import StrNoEditionXyz, FileXyz
from salomepy import utilsMedLoader as UML
import xyzpy.classFactoryXyz as CLFX

from PyQt5 import QtGui, QtWidgets

import xyzpy.loggingXyz as LOG
logger = LOG.getLogger()

verbose = False

###############################################################
class SourceFileCis(FileXyz):
  _envvars = "HOME WORKDIR4CASSIS".split()
  _title = 'Select source file'
  _filter = 'source files for mesh and groups (*.med)\nany files(*)'
  
  def getActionsContextMenu(self):
    actions = []
    actions.append( self._createAction('Browse', None, 'Browse', self.browseDialog, 'browsefile') )
    return actions

###############################################################
class GroupsInMeshXyz(_XyzConstrainBase):
  """
  class to evaluate groups from mesh med file
  all groups are in attributes as string str(value)
  """

  _attributesList = [ #list, not a dict because sequential order list is used in files Cis
    ("sourceFile", "SourceFileCis"),
    ("groups", "BaseFreeXyz"),#TODO : check synonymes for warning ?
  ]
  
  _helpDict = {
    "Simulations": ("list of Cassis cases data parameters", ""),
  }
  _icon = "cassispy.resources.maillage"

  
  def __init__(self):
    super(GroupsInMeshXyz, self).__init__()
    self._defautNameAsRoot = "GroupsInMesh"
    self.setIsCast(True)
    self._currentGroups = []
    self._currentMesh = None
    self._setAllAttributesList()

  def setDefaultValues(self):
    self.sourceFile = "" #copy... strange unique name attribute
    self.groups = BaseFreeXyz()
    self.setCurrentGroups()

  def __setattr__(self, name, value):
    super(GroupsInMeshXyz, self).__setattr__(name, value)
    if name == "sourceFile": 
      if verbose: print("sourceFile '%s'" % value)
      self.setCurrentGroups(sourceFile=value)    

  def setCurrentGroups(self, sourceFile=None):
    """
    set all name of groups variables as attributes from existing MedFileName
    """
    self.clearGroups() #clear previous
    if str(sourceFile) == "": return
    theMesh = UML.getMeshFromFile(sourceFile)
    self._currentMesh = theMesh # could serve later...
    if theMesh == None:
      #in case of bad local path from dataInformation do not display groups names
      logger.warning("file mesh and groups not found in %s" % sourceFile)
      return
    self._currentGroups = UML.getGroupsFromMesh(theMesh)
    for level, groups in self._currentGroups:
      lev = BaseFreeXyz()
      for g in groups:
        if g != "":
          lev.__setattr__(g, StrNoEditionXyz(""))
      self.groups.__setattr__(level, lev)

  def clearGroups(self):
    self.groups = BaseFreeXyz()

  def getSourceFile(self):
    return str(self.sourceFile) #copy for precaution

  def getGroupNames(self, typeGroup):
    """
    returns list of current groups of typeGroup
    example typeGroup ["Nodes", "Beams", "Shells", "Volumes"]
    """
    try:
      group = getattr(self.groups, typeGroup)
      return group.getAttributes()
    except Exception as e:
      logger.error( "problem getGroupNames: %s" % e )
      return []
  
  def getAllGroupNames(self, extend=True):
    res = []
    for typeGroup in self.groups.getAttributes():
      if extend: 
        res += self.getGroupNames(typeGroup)
      else:
        res.append(self.getGroupNames(typeGroup))
    return res
    
  def getActionsContextMenu(self):
    return self.sourceFile.getActionsContextMenu()

#factory pattern using xyzpy.utilsXyz._dictOfXyzClass
CLFX.appendAllXyzClasses( [SourceFileCis, GroupsInMeshXyz] )
