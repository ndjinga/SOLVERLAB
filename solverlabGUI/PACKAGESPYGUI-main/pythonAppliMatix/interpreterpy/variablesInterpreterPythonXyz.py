#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import traceback
import pprint as PP

from interpreterpy.interpreterPython import InterpreterPython as IPY
from xyzpy.baseXyz import _XyzConstrainBase, BaseFreeXyz, BaseXyz
from xyzpy.intFloatListXyz import StrNoEditionXyz
import xyzpy.classFactoryXyz as CLFX
from xyzpy.guiXyz.dialogXmlXyz import DialogXyz, QTextEditXyz
import xyzpy.globalVariables as GV
import xyzpy.loggingXyz as LOG

from PyQt5 import QtGui
from PyQt5 import QtCore, QtGui, QtWidgets as QTW

logger = LOG.getLogger()
verbose = False


###############################################################
class VariablePythonXyz(StrNoEditionXyz):
  pass

_defaultSourceGlobal = r'''
##############################################
# initial default example code source python #
# you are free to modify at your convenience #
##############################################

"""
This is helpful to set your comments.
"""

# ever floating mode for division
from __future__ import division

MyFirstVar = 1/3   # is an example...

#... etc ...

'''

_defaultSourceElementary = r'''
##############################################
# initial default example code source python #
# you are free to modify at your convenience #
#                                            #
# use '_GV.etc' to access                    #
# of unique global variables namespace       #
##############################################

"""
This is helpful to set your comments.
"""

# ever floating mode for division
from __future__ import division

MyFirstVar = 1/3   # is an example...

#... etc ...

'''

###############################################################
class VariablesInterpreterPythonXyz(BaseFreeXyz):
  """
  class to evaluate variables from python code
  all variables are in attributes as string str(value) of all python types
  """

  _icon = "python"
  _inexisting = "Inexisting"
  _global = False

  _defaultSource = _defaultSourceElementary

  def __init__(self):
    super(VariablesInterpreterPythonXyz, self).__init__()
    self._defautNameAsRoot = "Variables"
    self.setSource(self._defaultSource) #strange unique name attribute
    if verbose: print("VariablesInterpreterPythonXyz.__init__ setCurrentVariables(_defaultSource)")

  def __setattr__(self, name, value):
    super(VariablesInterpreterPythonXyz, self).__setattr__(name, value)
    if name == "X__source__":
      self.setCurrentVariables()

  def setCurrentVariables(self, source=None):
    """
    set all namespace variables as attributes from existing or new source python
    """
    self.setSource(source)
    if self._global:
      aDict = {} # forbid loop auto reference
    else:
      aDict = {"_GV": GV.getGlobalVariables(".GlobalVariables")}
    ip = IPY(aDict)
    ip.runcode(self.getSource())
    ipvars = ip.getVars()
    if self._global:
      _ipvars = BaseXyz()
      _ipvars.setAttributesFromDict(ipvars)
      GV.setGlobalVariables(".GlobalVariables", _ipvars)

    # print("*** variables %s on treePyName '%s'" % (list(ipvars.keys()), self.getTreePyName()))
    # logger.warning("setCurrentVariables:\n%s" % PP.pformat(ipvars))

    self.clearCurrentVariables() #clear previous
    
    if verbose: 
      print("setCurrentVariables %s" % list(ipvars.keys()))
      for i in traceback.format_stack(): 
        print("  ### stack: %s" % i)
        break
    
    for k, value in list(ipvars.items()):
      self.__setattr__(k, VariablePythonXyz(str(value)))

    self.__tooltip__ = ip.get__doc__()

  def getAppendingToolTip(self):
    """
    used for append in tooltip of parent in tree
    TODO not used for the moment 1904
    """
    return self.__tooltip__

  def clearCurrentVariables(self):
    toDel = self.getAttributes()
    for k in toDel:
      if k != "X__source__": delattr(self, k)
    
  def setSource(self, source):
    if source != None:
      self.X__source__ = VariablePythonXyz(str(source))
    else:
      pass #do nothing as choice to keep previous self.X__source__
    return

  def getSource(self):
    return str(self.X__source__) #copy for precaution

  def getVariablesNames(self):
    """
    returns list of current attributes from self.X__source__, sorted 
    (as no way from namespace dict)
    """
    self.setCurrentVariables()
    return sorted([str(i) for i in self.getAttributes() if i != "X__source__"]) #copy for precaution
    
  def getFloatValueByName(self, name):
    try:
      return float(getattr(self, name))
    except Exception as e:
      logger.error( "WARNING : Type of variable {} is not float : {}".format(getattr(self, name), e))
      return 0

  def getActionsContextMenu(self):
    """no browse"""
    if verbose: print("%s.getContextMenu" % self._className)
    actions = []
    actions.append( self._createAction('Edit variables', None, 'Edit source python for variables', self.createEditorSource, 'editor') )
    return actions

  """
    if True:
      pph = pparent.geometry().height()
      ppx = pparent.geometry().x()
      ppy = pparent.geometry().y()
      dw = aDialog.width()
      dh = aDialog.height()
      parentPos = pparent.mapToGlobal(pparent.pos())
      px = 100 #parentPos.x() #parentPos.x() + pparent.width()/2 - dw/2
      py = 200 #parentPos.y() #parentPos.y() + pparent.height()/2 - dh/2
      dx = dw
      dy = dh
      #aDialog.setGeometry(px, py, dx, dy) 
    while True:
      print "ppparent", aDialog.parent(), px, py, dx, dy
      #aDialog.setGeometry( px+100, py+100, dw, dh )
      #aDialog.move(self.getDesktop().rect().bottomLeft())
      #aDialog.setGeometry(px, py, dx, dy)
      aDialog.exec_()
  """

  def createEditorSource(self, parent=None):
    controller = self.getController()
    aDialog = DialogXyz(parent=self.getDesktop()) # open in center desktop
    aWidget = QTextEditXyz()
    aWidgetErr = QTextEditXyz()
    aDialog.setUpWidgetLayout([aWidget, aWidgetErr])
    aWidgetErr.hide()
    aWidget.setValue(self.getSource())
    aDialog.setMinimumSize(500, 450)
    aDialog.setWindowTitle(self.getTreePyName())
    while True:
      aDialog.exec_()
      if aDialog.choice == "Cancel": return True

      newSource = aWidget.getValue()
      if verbose: 
        print("VariablesInterpreterPythonXyz.createEditor: new source:'\n%s'" % newSource)
      
      # test newSource is correct
      if self._global:
        aDict = {} # forbid loop auto reference
      else:
        aDict = {"_GV": GV.getGlobalVariables(".GlobalVariables")}
      ip = IPY(aDict)
      ip.runcode(newSource)
      if ip.isOk(): break
      #print "Ooops:\n",ip.getStdErr()
      aWidgetErr.setValue(ip.getStdErr())
      aWidgetErr.show()
      # TODO fix why??? if not QMessageBox, next aDialog.exec_() is NOT in center desktop
      # may be qt loop event flush ?
      QTW.QMessageBox.warning(self.getDesktop(), "error in source", ip.getStdErr())

    if controller == None:
      print("ERROR: no controller for action Edit variables on %s" % self._className)
    else:
      controller.setModelItemValueSignalList.emit( ["%s.setCurrentVariables(args[1])" % self.getTreePyName(), newSource ] )
      expanded = [self.getTreePyName()]
      controller.ExpandSignal.emit(expanded)
    return True

  """# no direct edition
  def createEditor(self, parent):
    controller = self.getControllerInParents(parent)
    aDialog = XyzQTextEdit(parent=parent)
    aDialog.setValue(self.getSource())
    aDialog.setMinimumSize(200, 300)
    return aDialog
  """

  def createEditor(self, parent):
    self.createEditorSource()
    return None

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    if nameAttr == "X__source__": return True
    return False

###############################################################
class VariablesGlobalInterpreterPythonXyz(VariablesInterpreterPythonXyz):
  _global = True
  _defaultSource = _defaultSourceGlobal
  pass


#factory pattern using xyzpy.utilsXyz._dictOfXyzClass
CLFX.appendAllXyzClasses( [
  VariablePythonXyz,
  VariablesInterpreterPythonXyz,
  VariablesGlobalInterpreterPythonXyz
] )
