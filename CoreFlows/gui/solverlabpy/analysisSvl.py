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
import time
import platform
import random
from PyQt5 import QtGui, QtCore, QtWidgets as QTW
import subprocess as SP
import pprint as PP
import debogpy.debug as DBG

from xyzpy.baseXyz import _XyzConstrainBase, ListOfBaseXyz

import xyzpy.utilsXyz as UXYZ
import xyzpy.classFactoryXyz as CLFX
import xyzpy.loggingXyz as LOG
logger = LOG.getLogger()

#set classes Xyz in factory xyzpy.utilsXyz
import xyzpy.intFloatListXyz as IFLX
import interpreterpy.variablesInterpreterPythonXyz
import solverlabpy.solverlabSettings as USET
import solverlabpy.solverlabFilePatterns as UFPA
import solverlabpy.caseSolverlab

import solverlabpy.configSvl as CFGSVL

verbose = False


###############################################################
# common methods
###############################################################
def drawUnivariate(self):
  """common method draw univariate for classes as AttributeSvl or ExpressionSvl"""
  dataFile = self.getFileInParent() 
  expressions = self.getExpressionsInParent() 
  realDat = os.path.expandvars(dataFile.name) #absolute name    
  pars = [self.getAttributeName()]
  expr = expressions.getExpressions()
  fils = [""]
  opts = [""]
  res = UFPA.visualize_TDS_univariate(realDat, pars, EXPR=expr, DIVIDE=None, FILS=fils, OPTS=opts)
  self.getController().showTabByName("Canvas")
  return

###############################################################
def drawGraphic(self):
  """common method draw histogram for classes as AttributeSvl or ExpressionSvl"""
  dataFile = self.getFileInParent() 
  expressions = self.getExpressionsInParent() 
  realDat = os.path.expandvars(dataFile.name) #absolute name    
  pars = [self.getAttributeName()]
  expr = expressions.getExpressions()
  fils = [""]
  opts = [""]
  print("drawGraphic EXPR %s" % expr)
  res = UFPA.visualize_TDS_graphic(realDat, pars, EXPR=expr, DIVIDE=None, FILS=fils, OPTS=opts)
  self.getController().showTabByName("Canvas")
  return

###############################################################
def getDefaultRowColumn(nb):
  """return (row, column)"""
  if nb <= 1: return (1,1)
  if nb == 2: return (1,2)
  nbrow = int(nb**.5)
  if nbrow*nbrow < nb: 
    nbrow += 1
  return (nbrow , nbrow)

###############################################################
def getCurrentRowColumn(indiceCurrent, nbmaxGrid, rowBefore=True):
  """
  returns current (row, column) in nbmax elements in grid
  for indice current (0 to nbmaxGrid)
  """
  nbr, nbc = getDefaultRowColumn(nbmaxGrid)
  #print "getCurrentRowColumn", nbmaxGrid, nbr, nbc
  if rowBefore:
    row = indiceCurrent%nbr
    column = indiceCurrent/nbr
  else:
    column = indiceCurrent%nbc
    row = indiceCurrent/nbc  
  return (row, column)

###############################################################
def _test():
  for nb in range(0, 25):
    print("_test nb %s" % getDefaultRowColumn(nb))
  for nb in range(0, 17):
    print("_test nb %s %s" % (nb, getDefaultRowColumn(nb)))
    for i in range(0, 25):
      print("_test i %s %s" % (i, getCurrentRowColumn(i, nb)))


###############################################################
# classes
###############################################################
class DataInformationsSvl(_XyzConstrainBase):
  """
  general informations about
  """
  _attributesList = [ #list, not a dict because sequential order list is used in files
    ("name", "StrXyz"),
    ("directory", "DirectoryXyz"),
    ("release", "ReleaseXyz"),
  ]
  _icon = "datainformation"

  _helpDict = {
    "name": (u"As name of case directory", u""),
    "directory": (u"As parent directory of case directory", u""),
    "release": (u"Release of solverlab files (major.minor.patch)", u""),
  }

  _defaultVersion = "1.0.0"

  def __init__(self):
    super(DataInformationsSvl, self).__init__()
    self.setIsCast(True)

  def setDefaultValues(self):
    aDir = USET.getVar("_SOLVERLABGUI_WORKDIR")
    aName = "study_" + UXYZ.getDateTimeNow()
    #aDirName = os.path.join(aDirw, aName)
    #aDir, aName = os.path.split(os.path.realpath(aDirName))
    self.name = aName
    self.directory = aDir
    self.release = self._defaultVersion
    #self.dateFile = "" #no file saved time.ctime(os.path.getmtime(fileName))

  def getVersion(self):
    return ("%i.%i.%i" %(self.major, self.minor, self.patch))

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    return CFGSVL.isHidden(self, nameAttr)

  def getEtudeWorkdir(self):
    """return as os.path.join(directory,name)"""
    return os.path.join(str(self.directory), str(self.name))

  def getEtudeWorkdirBrut(self):
    """return as ${SOLVERLABGUI_WORKDIR}"""
    r1 = os.path.join(str(self.directory), str(self.name))
    r2 = os.path.join("${SOLVERLABGUI_WORKDIR}", str(self.name))
    r1Exp = os.path.expandvars(r1)
    r2Exp = os.path.expandvars(r2)
    if r1Exp != r2Exp:
      logger.warning("incoherency '%s'\n  different of '%s'" % (r2Exp, r1Exp)) 
    return os.path.join(str(self.directory), str(self.name))


###############################################################
class AttributeDataFrameSvl(IFLX.StrXyz):
  _icon = "solverlabpy.resources.leaf_s"

  def getActionsContextMenu(self):
    #as contents of file no edition etc...
    #actions = super(AttributeDataFrameSvl, self).getActionsContextMenu()
    actions = [] #no modify menu
    actions.append( self._createAction('draw univariate', None, "Execute univariate draw for Uranie TDataServer", self.drawUnivariate, "draw") )
    actions.append( self._createAction('draw histogram', None, "Execute graphic draw for Uranie TDataServer", self.drawGraphic, "draw") )
    return actions

  def getAttributeName(self):
    param = str(self)
    res = "NoName"
    tmp = param.split('#COLUMN:')[1]
    tmp = tmp.split('|')[0]
    res = tmp.strip()
    return res

  def getFileInParent(self):
    return  self.parentAsAttribute().parentAsAttribute()

  def getExpressionsInParent(self):
    return  self.parentAsAttribute().parentAsAttribute().expressions

  def drawUnivariate(self):
    return drawUnivariate(self)

  def drawGraphic(self):
    return drawGraphic(self)

  """
  def drawCppDataServer(self):
    aFile = UFPA.execReplaces("visualize_TDS_univariate_@DATETIME@.C")
    dataFile = self.getFileInParent()
    allpars = dataFile.getAllAttributesName()
    
    realDat = os.path.expandvars(dataFile.name) #absolute name
    realDirDat, baseNameDat = os.path.split(realDat)
    
    realWorkdir = self.getRoot().getEtudeWorkdirExpanded()

    if realDirDat == os.path.join(realWorkdir, "data"):
      #use relative path name macro in "../macros", data in "../data"
      usedDat = os.path.join("..", "data", baseNameDat)
    else:
      #use absolute path name
      usedDat = realDat #use absolute path name

    if True: #try:
      par1 = self.getAttributeName()
      fil1 = "(%s > 0)" % par1
      opt1 = ""
    else: #except:
      par1 = "??"
      fil1 = "(?? > 0)"
      opt1 = "??"

    replaces = [
      ("@FILEC@", aFile),
      ("@FILEDAT@", usedDat), #os.path.expandvars(self.name)),
      ("@ALLATTRIBUTES@", ":".join(allpars)),
      ("@ATT1@", par1),
      ("@FIL1@", fil1),
      ("@OPT1@", opt1),
    ]

    #print TODO set generate macro or execute without is choice... on QDialog
    res = UFPA.visualize_TDS_univariate(realDat, par1, DIVIDE=None, FILS=fil1, OPTS=opt1)
    return

    code = UFPA.getFilePatterns("visualize_TDS_univariate.C", replaces)
    aBrutFile = os.path.join(self.getRoot().getEtudeWorkdirBrut(), "macros", aFile)
    aRealFile = os.path.expandvars(aBrutFile)

    controller = self.getController()
    message = "\ndrawDataServer create file '%s':\n%s" % (aRealFile ,code)
    controller.centralLogView.insertText(message)

    cmd = ".Analysis.macroManager.macros.append( args[1] )"
    arg1 = aBrutFile
    
    with open(aRealFile, "w") as f:
      f.write(code)
    controller.ExecRootProcessFileSignal.emit(aRealFile)
    controller.setModelItemValueSignalList.emit( [ cmd, arg1 ] )
    controller.RefreshSolverlabModelSignal.emit(None)
    return
    
  def drawPythonDataServer(self):
    aFile = UFPA.execReplaces("visualize_TDS_univariate_@DATETIME@.py")
    dataFile = self.getFileInParent()
    allpars = dataFile.getAllAttributesName()
    
    realDat = os.path.expandvars(dataFile.name) #absolute name
    realDirDat, baseNameDat = os.path.split(realDat)
    
    realWorkdir = self.getRoot().getEtudeWorkdirExpanded()

    if realDirDat == os.path.join(realWorkdir, "data"):
      #use relative path name macro in "../macros", data in "../data"
      usedDat = os.path.join("..", "data", baseNameDat)
    else:
      #use absolute path name
      usedDat = realDat #use absolute path name

    if True: #try:
      par1 = self.getAttributeName()
      fil1 = "(%s > 0)" % par1
      opt1 = ""
    else: #except:
      par1 = "??:??"
      fil1 = "(?? > 0) && (?? > 0)"
      opt1 = "??"

    replaces = [
      ("@FILEC@", aFile),
      ("@FILEDAT@", usedDat), #os.path.expandvars(self.name)),
      ("@ALLATTRIBUTES@", ":".join(allpars)),
      ("@ATT1@", par1),
      ("@FIL1@", fil1),
      ("@OPT1@", opt1),
    ]
    code = UFPA.getFilePatterns("visualize_TDS_univariate.py", replaces)
    aBrutFile = os.path.join(self.getRoot().getEtudeWorkdirBrut(), "macros", aFile)
    aRealFile = os.path.expandvars(aBrutFile)

    controller = self.getController()
    message = "\ndrawDataServer create file '%s':\n%s" % (aRealFile ,code)
    controller.centralLogView.insertText(message)

    cmd = ".Analysis.macroManager.macros.append( args[1] )"
    arg1 = aBrutFile
    
    with open(aRealFile, "w") as f:
      f.write(code)

    controller.ExecRootProcessFileSignal.emit(aRealFile)
    controller.setModelItemValueSignalList.emit( [ cmd, arg1 ] )
    controller.RefreshSolverlabModelSignal.emit(None)
    return
  """


###############################################################
class ListOfAttributeSvl(ListOfBaseXyz):
  _allowedClasses = [AttributeDataFrameSvl]
  _icon = "solverlabpy.resources.tree_s"

  def getAllAttributesName(self):
    res = [i.getAttributeName() for i in self]
    return res  


###############################################################
class ExpressionSvl(IFLX.ExpressionXyz):
  _icon = "solverlabpy.resources.leaf_s"

  def getActionsContextMenu(self):
    #as contents of file no edition etc...
    actions = super(ExpressionSvl, self).getActionsContextMenu() #get modify
    actions.append( self._createAction('draw univariate', None, "Execute univariate draw for Uranie TDataServer", self.drawUnivariate, "draw") )
    actions.append( self._createAction('draw histogram', None, "Execute graphic draw for Uranie TDataServer", self.drawGraphic, "draw") )
    return actions

  def getFileInParent(self):
    return  self.parentAsAttribute().parentAsAttribute()

  def getExpressionsInParent(self):
    return  self.parentAsAttribute().parentAsAttribute().expressions

  def drawUnivariate(self):
    return drawUnivariate(self)

  def drawGraphic(self):
    return drawGraphic(self)

  #assume correct syntax
  def createEditorData(self, parent=None):
    if verbose: print("%s.createEditorData parent %s" % (self._className, parent))
    
    controller =  self.getController()
    aDialog = self.getEditDialog()
    aDialog.setModal(True)
    aDialog.closeOnApply = True

    for i in range(6): #no more 6 isNameUnique error
      aDialog.exec_()
      #print "Test Expression Svltest validity", aDialog.choice
      if aDialog.choice == "Cancel": return True
      if aDialog.choice == "Apply": 
        newValue = aDialog.getLocalModel().getValue()
        if verbose: print("ExpressionSvl newValue %s" % newValue)
        name, value = newValue.split("=")
        if self.getName() == name: #supposed overriding myself, so name stay unique
          ok = True
        else:
          ok = self.isNameUnique(name)
        if not ok:
          QTW.QMessageBox.warning(parent, "warning", "left-hand '%s' of expression is used yet." % name)
          continue
        ok = self.isValidExpression(value)
        if not ok:
          QTW.QMessageBox.warning(parent, "warning", "right-hand '%s' of expression is incorrect." % value)
          continue
        break
    
    if not ok: return True

    if verbose: 
      print("%s.createEditorData: new value: '%s'" % (self._className, newValue))

    if controller == None:
      logger.error("%s.createEditorData: no controller for action View/Edit" % self._className)
    else:
      #print "******** ExpressionXyz: %s = '%s'" % (self.getTreePyName(), newValue)
      #controller.setModelItemValueSignal.emit( "%s = '%s'" % (self.getTreePyName(), newValue) )
      controller.setModelItemValueSignalList.emit( ["%s = args[1]" % self.getTreePyName(), self.__class__(newValue) ] )
    return True

  def isNameUnique(self, name):
    """
    test unicity of name in DataSvl
    search in lists DataSvl.attributes and DataSvl.expressions
    """
    allpars = self.getAllAttributesNameInParents()
    if name in allpars: return False
    return True

  def getAllAttributesNameInParents(self):
    #dataManager = self.getFirstInParents("dataManager")
    expressions = self.parentAsAttribute()
    dataSvl = expressions.parentAsAttribute()
    allpars = dataSvl.attributes.getAllAttributesName()
    allpars.extend(expressions.getAllAttributesName())
    return allpars

  def getName(self):
    try:
      res, _ = self.split("=")
    except:
      res = ""
    res = res.strip() 
    return res

  def getAttributeName(self):
    """synonym for coherency with AttributeSvl"""
    return self.getName()

  def isValidExpression(self, value):
    """
    test if expression is valid in a uranie TDataServer
    """

    """
    // uranie C++ example
    // Create a TDataServer
    TDataServer * tds = new TDataServer();

    tds->addAttribute("tu");
    tds->addAttribute("tl");
    tds->addAttribute("hu");
    tds->addAttribute("hl");

    // Create Tuple
    tds->createTuple();
    TString sFormula = "tu+tl+5*hu*hl";
   
    // tds->addAttribute("toto", sFormula.Data());
    TTreeFormula * inform = new TTreeFormula("toto", sFormula.Data(), tds->getTuple());
    cout << " *** GetNdim[" << sFormula << "][" << inform->GetNdim() << "]" << endl;
    if (inform->GetNdim() != 0) {
      cout << " **** OK[" << sFormula << "]" << endl;
    } else {
      cout <<  " ************************** PROBLEM[" << sFormula << "]" << endl;
    }
    return();
    """

    import URANIE
    # Create a TDataServer
    tds = URANIE.DataServer.TDataServer()
    allpars = self.getAllAttributesNameInParents()
    for i in allpars:
      tds.addAttribute(str(i), 1., 2.)
    tds.createTuple()
    inform = URANIE.ROOT.TTreeFormula("tttesttt", value, tds.getTuple())
    #print "inform.GetNdim()", inform.GetNdim()
    if inform.GetNdim() == 0: 
      return False
    return True


###############################################################
class ListOfExpressionSvl(ListOfBaseXyz):
  _allowedClasses = [ExpressionSvl]
  _icon = "solverlabpy.resources.tree_s"
  pass

  def getExpressions(self):
    return [str(i) for i in self]

  def getAllAttributesName(self):
    res = [i.getName() for i in self]
    return res  


'''
###############################################################
class ComboBoxSvl(QTW.QComboBox):

  appendSignal = QtCore.pyqtSignal(str, str)
  
  def __init__(self, *args, **kwargs):
    super(ComboBoxSvl, self).__init__(*args, **kwargs)
    self._editText = ":"
    self.currentIndexChanged.connect(self.on_indexChanged)
    self.editTextChanged.connect(self.on_editTextChanged)
    self.appendSignal.connect(self.on_append)
    self.setCurrentIndex(-1)
    self.clearEditText()

  def on_indexChanged(self, index):
    if (index == -1): return
    ct = str(self._editText) #self.lineEdit().text())
    print("\non_indexChanged %s // %s" % (index, ct) )
    self.appendSignal.emit(ct, str(self.itemText(index)))
    self.setCurrentIndex(-1)

  def on_append(self, currentText, newItem):
    item = "%s:" % newItem
    newText = currentText + item
    print("on_append %s + %s -> %s" % (currentText, newItem, newText))
    self._editText = newText
    self.lineEdit().setText(newText)

  def on_editTextChanged(self, text):
    if text == "": self.lineEdit().setText(self._editText)
    let = str(self.lineEdit().text())
    print("\non_editTextChanged %s // %s // %s" % (text, let, self._editText))
'''


###############################################################
class DialogAttributesSvl(QTW.QDialog):

  def __init__(self, *args, **kwargs):
    super(DialogAttributesSvl, self).__init__(*args, **kwargs)
    bbox = QTW.QDialogButtonBox
    self.bbox = bbox(bbox.Reset | bbox.Apply | bbox.Cancel)
    self.bbox.button(bbox.Reset).clicked.connect(self.on_reset)
    self.bbox.button(bbox.Apply).clicked.connect(self.on_apply)
    self.bbox.button(bbox.Cancel).clicked.connect(self.on_cancel)
    self.gbox, self.vl = self._createGbox("Attributes")
    #self.gbox.setFlat(true);
    layout = QTW.QVBoxLayout()
    layout.addWidget(self.gbox)
    layout.addWidget(self.bbox)
    self.setLayout(layout)
    self.setWindowTitle('Choose draw attribute(s)')
    self.resize(300,200)

  def _createGbox(self, name):
    """create gbox with vertical layout into"""
    gbox = QTW.QGroupBox(name)
    layout = QTW.QGridLayout() #QVBoxLayout()
    gbox.setLayout(layout)
    return (gbox, layout)

  def addItems(self, data):
    self.checks = []
    gbox, vl = self._createGbox("Current attributes selection")
    self.value = QTW.QLabel("")
    vl.addWidget(self.value)
    nbmax = len(data.attributes) + len(data.expressions)
    nbr, nbc = getDefaultRowColumn(nbmax)
    self.vl.addWidget(gbox, 0, 0, 1, nbc) #extend all column(s)
    ii = 0
    for i in data.attributes:
      c = QTW.QCheckBox(i.getAttributeName())
      c.setCheckable(True)
      c.setToolTip(str(i))
      c.setObjectName(str(i))
      #warning first param status is overriden by emit
      c.stateChanged.connect(lambda status=None, sender=c: self.on_checked(status, sender))
      c.setChecked(False)
      self.checks.append(c)
      irow, icol = getCurrentRowColumn(ii, nbmax)
      #print i.getAttributeName(), ii, nbmax, (irow+1, icol)
      self.vl.addWidget(c, irow+1, icol) #skip first row where gbox
      ii += 1
    for i in data.expressions:
      c = QTW.QCheckBox(i.getAttributeName())
      c.setCheckable(True)
      c.setToolTip(str(i))
      c.setObjectName(str(i))
      #warning first param status is overriden by emit
      c.stateChanged.connect(lambda status=None, sender=c: self.on_checked(status, sender))
      c.setChecked(False)
      self.checks.append(c)
      irow, icol = getCurrentRowColumn(ii, nbmax)
      #print i.getAttributeName(), ii, nbmax, (irow+1, icol)
      self.vl.addWidget(c, irow+1, icol) #skip first row where gbox
      ii += 1

  def on_checked(self, status, sender):
    if verbose: print("on_checked %s %s" % (status, sender.objectName()))
    item = ": %s :" % sender.text()
    currentText = self.value.text()
    if sender.isChecked():
      currentText += item
    else:
      currentText = currentText.replace(item, "")
    self.value.setText(currentText)

  def on_reset(self):
    #print("TODO %s.on_reset" % self.__class__.__name__)
    self.value.setText("")
    for i in self.checks: i.setChecked(False)
  
  def on_apply(self):
    #print("TODO %s.on_apply" % self.__class__.__name__)
    #self.getValue()
    self.accept()
    self.close()
  
  def on_cancel(self):
    #print("TODO %s.on_cancel" % self.__class__.__name__)
    self.close()

  def getValue(self):
    res = str(self.value.text())
    res = [i.strip() for i in res.split(":") if i != ""]
    return res

  #def setValue(self, value):
  #  return self.setValue(value)



###############################################################
class ComboBoxSvl(QTW.QComboBox):
  
  def __init__(self, *args, **kwargs):
    super(ComboBoxSvl, self).__init__(*args, **kwargs)
    self.setCurrentIndex(-1)
    self.clearEditText()

  def contextMenuEvent(self, event):
    print("ComboBoxSvl.contextMenuEvent")
 
  def getValue(self):
    res = str(self.lineEdit().text())
    return res


###############################################################
class ComboBoxAttributesSvl(QTW.QComboBox):
 
  def __init__(self, *args, **kwargs):
    super(ComboBoxAttributesSvl, self).__init__(*args, **kwargs)
    self.setCurrentIndex(-1)
    self.clearEditText()

  def contextMenuEvent(self, event):
    if verbose: print("ComboBoxAttributesSvl.contextMenuEvent")
    dial = DialogAttributesSvl(self)
    dial.addItems(self.data)
    res = dial.exec_()
    if res == dial.Accepted:
      values = dial.getValue()
      newText = " : ".join(values)
      if verbose: print("ComboBoxAttributesSvl result: '%s'" % newText)
      self.lineEdit().setText(newText)
      
  def getValue(self):
    res = str(self.lineEdit().text())
    res = res.replace(" ", "") #strip values as do not like blanks
    return res

     
    
###############################################################
class DrawDialogDataSvl(QTW.QDialog):
  def __init__(self, *args, **kwargs):
    super(DrawDialogDataSvl, self).__init__(*args, **kwargs)
    bbox = QTW.QDialogButtonBox
    self.bbox = bbox(bbox.Reset | bbox.Apply | bbox.Cancel | bbox.Help)
    self.bbox.button(bbox.Reset).clicked.connect(self.on_reset)
    self.bbox.button(bbox.Apply).clicked.connect(self.on_apply)
    self.bbox.button(bbox.Cancel).clicked.connect(self.on_cancel)
    self.bbox.button(bbox.Help).clicked.connect(self.on_help)
    self._helpHtml = "https://root.cern.ch/root/html534/THistPainter.html"

  def _setCombo(self, name, aClass, layout):
    combo = aClass()
    combo.setEditable(True)
    lab = QTW.QLabel(name)
    layout.addWidget(lab)
    layout.addWidget(combo)
    combo.data = self.data
    return combo

  def getDefaultFilters(self):
    res = """
( xxx > 0 )
( xxx < 0 ) && ( yyy > 0 )
( xxx < 0 ) || ( yyy > 0 )"""
    return res.split("\n")

  def getDefaultOptions(self):
    res = """
LP
PARA
CONT1
CONT5
SURF1
SURF5"""
    return res.split("\n")

  def _createCheck(self, name, tooltip, layout):
    c = QTW.QCheckBox(name)
    c.setCheckable(True)
    c.setToolTip(tooltip)
    c.setObjectName(name)
    c.setChecked(False)
    layout.addWidget(c)
    return c

  def _createWidget(self, aLayout):
    w = QTW.QWidget()
    w.setLayout(aLayout)
    return w

  def _initLayout(self, kwargs):
    self.data = kwargs["data"]
    if verbose: print("DrawDialogDataSvl._initLayout %s" % PP.pformat(kwargs))
    layout = QTW.QVBoxLayout()
    self.typeOfVisualize = self._setCombo("type of visualisation :", ComboBoxSvl, layout)
    self.typeOfVisualize.addItems(UFPA.getTypesOfVisualize())
    self.attributes = self._setCombo("attribute(s) to draw:", ComboBoxAttributesSvl, layout)
    allpars = self.data.getAllAttributesName()
    self.attributes.setToolTip("right-click to get smart dialog")
    self.attributes.addItems(allpars)
    self.filters = self._setCombo("filters(s) on attribute(s):", ComboBoxSvl, layout)
    self.filters.addItems(self.getDefaultFilters())
    self.options = self._setCombo("options(s) for draw:", ComboBoxSvl, layout)
    self.options.addItems(self.getDefaultOptions())
    
    hlayout = QTW.QHBoxLayout()
    self.checkCreateMacroC = self._createCheck("Create macro .C", "append corresponding editable macro file '.C' in macro directory", hlayout)
    self.checkCreateMacroPython = self._createCheck("Create macro .py", "append corresponding editable macro file '.py' in macro directory", hlayout)
    self.options.addItems(self.getDefaultOptions())
    
    layout.addWidget(self._createWidget(hlayout))
    layout.addWidget(self.bbox)
    self.setLayout(layout)
    self.setWindowTitle('Choose draw options')
    self.resize(400,200)
    
  def on_reset(self):
    #if verbose: print("%s.on_reset" % self.__class__.__name__)
    self.attributes.clearEditText()
    self.filters.clearEditText()
    self.options.clearEditText()
  
  def on_apply(self):
    #if verbose: print("%s.on_apply" % self.__class__.__name__)
    self.accept()
    self.close()
  
  def on_cancel(self):
    #if verbose: print("%s.on_cancel" % self.__class__.__name__)
    self.close()

  def on_help(self):
    nameBrowser = UXYZ.getBrowser()
    cmd = "%s %s &" % (nameBrowser, self._helpHtml)
    logger.info("on_help %s" % cmd)
    os.system(cmd)

  def getValue(self):
    typ = self.typeOfVisualize.getValue()
    att = self.attributes.getValue()
    fil = self.filters.getValue()
    opt = self.options.getValue()
    res = (typ, att, fil, opt)
    if verbose: 
      print("%s.getValue %s" % (self.__class__.__name__, PP.pformat(res)))
    return res

  def getValueAsDict(self):
    res = {}
    res["type"] = self.typeOfVisualize.getValue()
    res["attributes"] = self.attributes.getValue()
    res["filters"] = self.filters.getValue()
    res["options"] = self.options.getValue()
    return res

  #def setValue(self, value):
  #  print "TODO DrawDialogDataSvl setValue"
  #  return



###############################################################
class DataSvl(_XyzConstrainBase):
  """DataFrame"""
  _typesFiles = "*.dat *.js".split()
  _directory = USET.getVar("_SOLVERLABCODE_ROOT_DIR") #env var unix syntax

  _attributesList = [ #list, not a dict because sequential order list is used in files
    ("name", "FileViewerXyz"),
    ("drawFilter", "StrXyz"),
    ("attributes", "ListOfAttributeSvl"),
    ("expressions", "ListOfExpressionSvl"),
  ]
  _helpDict = {
    "name": (u"name of data frame file", u""),
    "drawFilter": (u"general filter for draw data frame", u""),
    "attributes": (u"attributes of data frame", u""),
    "expressions": (u"user attributes from expressions of attributes of data frame", u""),
  }
  _icon = "solverlabpy.resources.rootdb_s"

  def __init__(self):
    super(DataSvl, self).__init__()
    self.setIsCast(True)
    self.name = self._getClassFromAttributeName("name")("")
    self.drawFilter = self._getClassFromAttributeName("drawFilter")()
    self.attributes = self._getClassFromAttributeName("attributes")()
    self.expressions = self._getClassFromAttributeName("expressions")()
    self._dial = None

  def __setattr__(self, name, value):
    super(DataSvl, self).__setattr__(name, value)
    if name == "name": self.setAttributes()
    return

  def setAttributes(self):
    """
    append attributes from '#COLUMN_NAMES' from file .dat, 
    ordered as useful in uranie
    """
    params = self._getClassFromAttributeName("attributes")()
    if str(self.name) != "": #non defined file name is not error
      header = UFPA.getHeaderContentsFileDat(self.name)
      #ordered as useful in uranie
      headerUranie = UFPA.filterColumsNamesForUranie(header)
      for i in headerUranie:
        pari = params._allowedClasses[0](i)
        params.append(pari)
    self.attributes = params

  def getActionsContextMenu(self):
    actions = super(DataSvl, self).getActionsContextMenu()[0:2] #avoid useless insert etc
    actions.append( self._createAction('copy file in data directory', None, "Copy file in current directory data", self.copyFileInData, "copy") )
    #actions.append( self._createAction('draw graphic C', None, "Create and execute C graphic macro for Uranie TDataServer", self.drawCppDataServer, "draw") )
    #actions.append( self._createAction('draw graphic Python', None, "Create and execute Python graphic macro for Uranie TDataServer", self.drawPythonDataServer, "draw") )
    actions.append( self._createAction('draw something', None, "Execute dialog for draw Uranie TDataServer", self.drawDialogDataServer, "draw") )
    return actions
  
  def copyFileInData(self):
    controller = self.getController()    
    realDat = os.path.expandvars(self.name) #absolute name
    realDirDat, baseNameDat = os.path.split(realDat)
    realWorkdir = self.getRoot().getEtudeWorkdirExpanded()
    if realDirDat == os.path.join(realWorkdir, "data"):
      desktop = self.getDesktop()
      QTW.QMessageBox.warning(desktop, "warning", "file is already present in data directory:\n'%s' " % realDat)
      return False
    
    controller.UpdateEtudeDataSignal.emit(str(self.name))
    return True

  def drawDialogDataServer(self):
    desktop = self.getDesktop()    
    if self._dial == None:
      self._dial = DrawDialogDataSvl(desktop)
      self._dial._initLayout({"title":"TODO layout", "data": self})
    res = self._dial.exec_()
    if res != self._dial.Accepted: return
    value = self._dial.getValue()
    if verbose: print("DrawDialog process result: %s" % (PP.pformat(value)))
    TYPE, ATT, FIL, OPT = value
    FILEDAT = os.path.expandvars(self.name)
    aVisualize = UFPA.getVisualizeMethod(TYPE)
    expr = self.expressions.getExpressions()
    tds = aVisualize(FILEDAT, [ATT], EXPR=expr ,DIVIDE=None, FILS=[FIL], OPTS=[OPT], Verbose=True)
    controller = self.getController()
    controller.centralLogView.showTabByName("ROOTCanvasWidget")
    
    #create files macro if asked
    replaces = self._getReplacesForVisualizeMacro()
    if self._dial.checkCreateMacroC.isChecked():
      filePattern = "visualize_%s.C" % TYPE
      code = UFPA.getFilePatterns(filePattern)
      self._createMacroFile(filePattern, code, replaces)
      
    if self._dial.checkCreateMacroPython.isChecked():
      filePattern = "visualize_%s.py" % TYPE
      code = UFPA.getFilePatterns(filePattern)
      self._createMacroFile(filePattern, code, replaces)

  def _getReplacesForVisualizeMacro(self):
    realDat = os.path.expandvars(self.name) #absolute name
    realDirDat, baseNameDat = os.path.split(realDat)
    
    realWorkdir = self.getRoot().getEtudeWorkdirExpanded()

    if realDirDat == os.path.join(realWorkdir, "data"):
      #use relative path name macro in "../macros", data in "../data"
      usedDat = os.path.join("..", "data", baseNameDat)
    else:
      #use absolute path name
      usedDat = realDat #use absolute path name

    replaces = UFPA.getStdReplaces()
    allpars = self.getAllAttributesName()
    value = self._dial.getValue()
    TYPE, ATT, FIL, OPT = value
    EXPR = self._getExpressions(usedDat)

    replaces.extend( [
      ("@FILEDAT@", usedDat), #os.path.expandvars(self.name)),
      ("@ALLEXPRESSIONS@", EXPR),
      ("@ALLATTRIBUTES@", ":".join(allpars)),
      ("@ATT1@", ATT),
      ("@FIL1@", FIL),
      ("@OPT1@", OPT),
    ] )
    return replaces

  def _getExpressions(self, usedDat):
    exprs = self.expressions.getExpressions()
    if len(exprs) == 0:
      return '// tds->addAttribute("toto", "titi+tata");'
    else:
      res = ""
      for e in exprs:
        name, value = e.split('=')
        res += '  tds->addAttribute("%s", "%s");\n  ' % (name, value)
    return res[2:] #no 2 first whitespaces for indentation

  def _createMacroFile(self, filePattern, code, replaces):
    fileName, ext = os.path.splitext(filePattern)
    aFile = UFPA.execReplaces("%s_@DATETIME@%s" % (fileName, ext), replaces)
    replacesNew = [("@FILE@", aFile)]
    replacesNew.extend(replaces)
    aCode = UFPA.execReplaces(code, replacesNew)

    aBrutFile = os.path.join(self.getRoot().getEtudeWorkdirBrut(), "macros", aFile)
    aRealFile = os.path.expandvars(aBrutFile)

    controller = self.getController()
    message = "\ncreate macro file '%s':\n%s" % (aRealFile, aCode)
    controller.centralLogView.insertText(message)

    cmd = ".Analysis.macroManager.macros.append( args[1] )"
    arg1 = aBrutFile
    
    if verbose: print("_createMacroFile '%s' with replaces:\n%s" % (aRealFile, PP.pformat(replacesNew)))
    with open(aRealFile, "w") as f:
      f.write(aCode)

    #avoid execution... controller.ExecRootProcessFileSignal.emit(aRealFile)
    controller.setModelItemValueSignalList.emit( [ cmd, arg1 ] )
    controller.RefreshSolverlabModelSignal.emit(None)
    return

  def getAllAttributesName(self):
    res = []
    for p in self.attributes:
      res.append(p.getAttributeName())
    for p in self.expressions:
      res.append(p.getAttributeName())
    return res

  def getNameExpanded(self):
    return self.name.getNameExpanded()

###############################################################
class FileSvl(IFLX.FileViewerXyz):
  """initialize type of files extension for interest .C etc."""
  _typesFiles = "*.H *.C *.h *.c *.hxx *.cxx *.hpp *.cpp *.py *.bash *.so *.dat *.js".split()
  _directory = USET.getVar("_SOLVERLABCODE_ROOT_DIR") #env var unix syntax
  #_icon = #no icon
  pass

###############################################################
class UserFileSvl(IFLX.FileViewerXyz):
  """initialize type of files extension for user interest"""
  _typesFiles = "*".split()
  #_directory = USET.getVar("_SOLVERLABCODE_ROOT_DIR") #env var unix syntax
  #_icon = #no icon
  pass

###############################################################
class FunctionSvl(FileSvl):
  """initialize type of files extension for function files"""
  _typesFiles = "*".split()
  #_directory = USET.getVar("_ROOTSYS")+"/macros" #env var unix syntax
  _directory = USET.getVar("_SOLVERLABCODE_ROOT_DIR") + "/macros"  # env var unix syntax
  _icon = "cpp"
  pass

###############################################################
class LibrarySvl(FileSvl):
  """initialize type of files extension for libraries files .C .so"""
  _typesFiles = "*.C *.so".split()
  _directory = USET.getVar("_SOLVERLABCODE_ROOT_DIR")+"/lib" #env var unix syntax
  _icon = "solverlabpy.resources.branch_folder_s"
  pass


###############################################################
class MacroSvl(FileSvl):
  """initialize type of files extension for macro .C .py"""
  _typesFiles = "*.C *.py".split()
  _directory = USET.getSolverlabSysMacrosDir()
  _icon = "run" #"uraniepy.resources.folder"
  pass

###############################################################
class ListOfFileViewerXyz(ListOfBaseXyz):
  """
  base class, used only inheritage, modify _allowedClasses etc
  """
  _directory = "$HOME"
  _targetDirectory = "" #directory more than EtudeWorkdir as "macros" 
  _typesFiles = "*".split()
  _allowedClasses = [IFLX.StrXyz] #by default, modified by inheritage 
  _icon = "solverlabpy.resources.folder"
  _browseMessage = "Select files"

  def getDirectory(self):
    """could be dynamic in inheritage"""
    return self._directory
  
  def getTargetDirectory(self):
    """could be dynamic in inheritage"""
    return self._targetDirectory

  def getActionsContextMenu(self):
    actions = super(ListOfFileViewerXyz, self).getActionsContextMenu()
    for iClass in self._allowedClasses:
      aClass = iClass() #instancier pour avoir le nom str(iClass) durdur!
      name = aClass._defautNameAsRoot
      aText = 'Append item(s) %s' % name
      aIcon = aClass._icon
      actions.append( self._createAction(aText, None, aText, lambda status, aClass=iClass: self.addItemsSlot(status, aClass), aIcon) )
    return actions

  def addItemsSlot(self, status, aClass=None):
    """new pyqt5 state override in lambda"""
    return self.addItems(aClass)

  def addItems(self, aClass=None):
    selfClassName = self.__class__.__name__
    if aClass==None:
      logger.error( "%s.addItem %s not in %s" % (selfClassName, str(None), str(self._allowedClasses)) )
      return
    item = aClass()
    itemClassName = item.__class__.__name__
    logger.debug("%s.addItem standard by default with %s" % (selfClassName, itemClassName))
    if aClass not in self._allowedClasses:
      logger.error( "%s.addItem %s not in %s" % (selfClassName, itemClassName, str(self._allowedClasses)) )
      return
    self._itemBrowsing = item #to get type for browseViewerExecOnApply
    self.browseViewerDialog()
    
  def browseViewerDialog(self):
    desktop = self.getDesktop()
    controller = self.getController()    
    if controller == None:
      QTW.QMessageBox.warning(desktop, "warning", "%s: needs a controller for browseViewer" % self._title)
      return False

    #print "browse directory:", self._browseMessage, self._typesFiles, self.getDirectory()
    aDialog = controller.getExploreDir()
    strForWhat = self._browseMessage
    aDialog.lockSetDirRootPath(self.getDirectory(), filters=self._typesFiles)
    controller.showExploreDir()
    aDialog.lockForChangeModelController(strForWhat, self.browseViewerExecOnApply)
    return True

  def browseViewerExecOnApply(self, selectedFiles):
    """part of Apply on browseViewerDialog"""
    #print "%s.browseViewerExecOnApply" % self._className, selectedFiles
    if len(selectedFiles) == 0: 
      logger.warning("no file selected: as 'Cancel'")
      return False
    controller = self.getController()
    if controller == None:
      logger.error("no controller for %s action 'Apply'" % self._className)
      return False

    for name in selectedFiles:
      nameFile = str(name)
      if nameFile == "": continue #cancel
      cmd = "%s.append( args[1] )" % self.getTreePyName()
      try: #type immutable 
        arg1 = self._itemBrowsing.__class__(nameFile) #stuff set type (no casting) if appended in ListOf...
      except: #other types 
        arg1 = self._itemBrowsing.__class__()
        arg1.name = nameFile
      controller.setModelItemValueSignalList.emit( [ cmd, arg1 ] )
    controller.UpdateEtudeSignal.emit('all') #copy in local directory
    return True

  def addItem(self, aClass=None):
    """override method ListOfBaseXyz"""
    #print "create new file", self._patternName
    desktop = self.getDesktop()
    controller = self.getController()

    dial = QTW.QInputDialog
    adir = self.getRoot().getEtudeWorkdirBrut()
    valueInit = self._patternName
    basename, ok = dial.getText(
                     desktop, "choose new file name", 
                     "in %s" % adir, QTW.QLineEdit.Normal, valueInit)
    if not ok: return
    basename = str(basename)
    if basename == "": return
    name = os.path.join(self.getRoot().getEtudeWorkdirBrut(), self.getTargetDirectory(), basename)
    nameExp = os.path.expandvars(name)
    selfClassName = self.__class__.__name__
    if aClass==None:
      logger.error( "%s.addItem %s not in %s" % (selfClassName, str(None), str(self._allowedClasses)) )
      return
    try:
      item = aClass(name) #immutable as IFLX.FileViewerXyz
    except:
      item = aClass() #as DataSvl
      item.name = name
    itemClassName = item.__class__.__name__
    logger.debug("%s.addItem standard by default with %s" % (selfClassName, itemClassName))
    if aClass not in self._allowedClasses:
      logger.error( "%s.addItem %s not in %s" % (selfClassName, itemClassName, str(self._allowedClasses)) )
      return

    #try:
    #  item.setDefaultValues()
    #except:
    #  logger.warning("%s.setDefaultValues problem" % itemClassName)
    #self.append(item)

    if controller != None:
      if verbose: print("there is controller and refresh views")
      controller.createFilePattern(nameExp, self._patternName)
      try: 
        item.name = name #set attributes file existing...
      except:
        pass
      self.append(item)
      controller.refreshModelViews()
      controller.UpdateEtudeSignal.emit('all') #copy in local directory
    return

  def getNamesExpanded(self):
    res = [i.getNameExpanded() for i in self]
    return res

  def getNoLocal(self):
    files = self.getNamesExpanded()
    etudeDirExp = self.getRoot().getEtudeWorkdirExpanded()
    res = []
    for name in files:
      if etudeDirExp not in name: res.append(name)
    return res
    

###############################################################
class DataManagerSvl(ListOfFileViewerXyz):
  _directory = DataSvl._directory
  _targetDirectory = "data" 
  _typesFiles = DataSvl._typesFiles
  _allowedClasses = [DataSvl]
  _icon = "solverlabpy.resources.rootdb_s"
  _browseMessage = "Select Data files"
  _patternName = "aUserData.dat"

  def getActionsContextMenu(self):
    actions = super(DataManagerSvl, self).getActionsContextMenu()
    actions.append( self._createAction('copy all files in data directory', None, "Copy all file in current directory data", self.copyAllFileInData, "copy") )   
    return actions

  def copyAllFileInData(self):
    controller = self.getController()
    controller.UpdateEtudeDataSignal.emit("all")
    return True

###############################################################
class ListOfFunctionSvl(ListOfFileViewerXyz):
  _directory = FunctionSvl._directory
  _targetDirectory = "macros" 
  _typesFiles = FunctionSvl._typesFiles
  _allowedClasses = [FunctionSvl]
  _icon = "cpp" #"solverlabpy.resources.c_src_t"
  _browseMessage = "Select function files"
  _patternName = "aUserFunctions.C"
  pass

###############################################################
class ListOfLibrarySvl(ListOfFileViewerXyz):
  _directory = LibrarySvl._directory
  _targetDirectory = "macros" 
  _typesFiles = LibrarySvl._typesFiles
  _allowedClasses = [LibrarySvl]
  _icon = "solverlabpy.resources.branch_folder_s"
  _browseMessage = "Select library files"
  _patternName = "aUserLibrary.C"
  pass

###############################################################
class ListOfMacroSvl(ListOfFileViewerXyz):
  _directory = MacroSvl._directory
  _targetDirectory = "macros" 
  _typesFiles = MacroSvl._typesFiles
  _allowedClasses = [MacroSvl]
  _icon = "run" #"solverlabpy.resources.folder"
  _browseMessage = "Select macro files"
  _patternName = "aLaunching.C"
  pass

###############################################################
class ListOfUserFileSvl(ListOfFileViewerXyz):
  _directory = UserFileSvl._directory
  _targetDirectory = "" #as root of etude directory
  _typesFiles = UserFileSvl._typesFiles
  _allowedClasses = [UserFileSvl]
  _icon = "user" 
  _browseMessage = "Select your files"
  _patternName = "aUserFile"

  def getDirectory(self):
    """dynamic override _directory"""
    return "/tmp"
  
  def getTargetDirectory(self):
    """dynamic override _targetDirectory"""
    return self._targetDirectory

###############################################################
class MacroManagerSvl(_XyzConstrainBase):
  _attributesList = [ #list, not a dict because sequential order list is used in files Xyz
    ("functions","ListOfFunctionSvl"),
    ("libraries", "ListOfLibrarySvl"),
    ("macros", "ListOfMacroSvl"),
  ]
  _icon = "solverlabpy.resources.macromanager"
  _helpDict = {
    "functions": ("appends from directory %s" % ListOfFunctionSvl._directory, ""),
    "libraries": ("appends from directory %s" % ListOfLibrarySvl._directory, ""),
    "macros": ("appends from directory %s" % ListOfMacroSvl._directory, ""),
  }
  
  def __init__(self):
    super(MacroManagerSvl, self).__init__()
    self.setIsCast(True)
    self._defautNameAsRoot = "MacroManagerSolverlab"
    self._setAllAttributesList()

  def setDefaultValues(self):
    self.functions.setDefaultValues()
    self.libraries.setDefaultValues()
    self.macros.setDefaultValues()

  def getNamesExpanded(self):
    res = self.functions.getNamesExpanded()
    res.extend(self.libraries.getNamesExpanded())
    res.extend(self.macros.getNamesExpanded())
    return res

  def getNoLocal(self):
    res = self.functions.getNoLocal()
    res.extend(self.libraries.getNoLocal())
    res.extend(self.macros.getNoLocal())
    return res

###############################################################
class HistoryFileManagerXyz(_XyzConstrainBase):
  """
  store in string for xml save
  history on some files, with hashing control
  as cp, mv etc... actions from gui/model actions
  one line by action
  """
  _attributesList = [ #list, not a dict because sequential order list is used in files Xyz
    ("history","StrXyz"),
  ]
  _icon = "sablier"

  def __init__(self):
    super(HistoryFileManagerXyz, self).__init__()
    self.setIsCast(True)
    self._defautNameAsRoot = "historyFileManager"
    self._setAllAttributesList()
    #self.example()

  def setDefaultValues(self):
    self.history = ""
  
  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    res = False
    #TODO if nameAttr in ["history"]:
    #  res = True
    return res

  def clearHistory(self):
    self.history = ""

  """def example(self):
    print("HistoryFileManagerXyz.example")
    self.appendHistoryCopyOf('ini.py','end.py')"""
    
  def getFileHash(self, aFile):
    import hashlib as HASH
    filename = os.path.expandvars(aFile)
    with open(filename, "rb") as f: 
      hashvalue = HASH.sha1(f.read()).hexdigest()
    return hashvalue[0:20] #pas trop nen faut

  def getCompleteFileName(self, aFile):
    filename = os.path.expandvars(aFile)
    res = "%s:%s" % (platform.node(), os.path.realpath(filename))
    return res

  def getIdentFile(self, aFile):
    res = "%s, %s, %s" % (aFile, self.getCompleteFileName(aFile), self.getFileHash(aFile))
    return res
    
  def appendHistoryCopyOf(self, originFile, newFile):
    new = self.getIdentFile(newFile)
    origin = self.getIdentFile(originFile)
    currentdate = UXYZ.getDateTimeNow()
    #action description comprehension as (almost) current language
    action = "(%s)IsCopyOf(%s)At(%s)" % (new, origin, currentdate)
    self.appendHistoryAction(action)

  def appendHistoryAction(self, action):
    #TODO print self.history
    #print("HistoryFileManagerXyz\n****** history:\n%s\n****** action:\n%s" % (self.history, action))
    self.history = str(self.history) + "\n" + str(action)



###############################################################
class AnalysisSvl(_XyzConstrainBase):

  _attributesList = [ #list, not a dict because sequential order list is used in files Xyz
    ("dataInformations", "DataInformationsSvl"),
    ("caseSolverlab", "CaseSvl"),
    #("dataManager", "DataManagerSvl"),
    #("macroManager", "MacroManagerSvl"),
    #("userFileManager", "ListOfUserFileSvl"),
    #("historyFileManager", "HistoryFileManagerXyz"),
    #("pythonManager", "VariablesInterpreterPythonXyz"),
  ]
  _icon = "run"
  #_helpDict = {
  #  "": ("", ""),
  #}
  
  def __init__(self):
    super(AnalysisSvl, self).__init__()
    self.setIsCast(True)
    self._defautNameAsRoot = "Analysis"
    self._setAllAttributesList()

  def setDefaultValues(self):
    self.dataInformations.setDefaultValues()
  
  def _checkAttribute(self, attribute):
    """if elif checks on possibles attributes in inherited classes"""
    #if attribute == "BinaryOutput": return self.check_IntO1Xyz(attribute)
    #non permissif tout doit etre teste
    if attribute not in self.getAttributes():
      return CheckBaseXyz(False, self.__class__.__name__+": unexpected '"+attribute+"' attribute to check")
    #permissif ok comme dans la virtual method, mais sans logger warning
    #logger.warning(self.__class__.__name__+": ok but no check for attribute: '"+attribute+"'")
    return CheckBaseXyz(True, "")

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    return CFGSVL.isHidden(self, nameAttr)

  def appendHistoryFileManager(self, action):
    self.historyFileManager.appendHistory(action)
    
  def getActionsContextMenu(self):
    actions = super(AnalysisSvl, self).getActionsContextMenu()
    # actions.append( self._createAction('versioning commit', None, 'git commit', self.gitCommit, "solverlabpy.resources.gitcommit") )
    # actions.append( self._createAction('gitk', None, 'gitk', self.gitkLaunch, "solverlabpy.resources.gitcommit") )
    # actions.append( self._createAction('create doc', None, 'create documentation', self.createDocLaunch, "doc") )
    # actions.append( self._createAction('update rootlogon', None, 'update rootlogon.C file', self.updateRootlogonLaunch, "run") )
    # actions.append( self._createAction('packaging', None, 'packaging', self.packageLaunch, "tgz") )
    #actions.append( self._createAction('run python code', None, 'run python code', self.runPythonCode, "run") )
    #actions.append( self._createAction('versioning commit', None, 'git commit', self.gitCommit, "solverlabpy.resources.gitcommit") )
    # actions.append( self._createAction('search URANIE method', None, 'search URANIE method by name', self.searchURANIEMethod, "search") )
    # actions.append( self._createAction('print env ROOT', None, 'print ROOT context', self.printROOTContext, "search") )
    # actions.append(self._createAction('post treatments', None, 'post treatments plots etc', self.postTreatments, "search"))
    return actions

  def gitCommit(self):
    #print("DataInformationsSvl.gitCommit")
    controller = self.getController()
    controller.GitCommitEtudeSignal.emit("all") 

  def gitkLaunch(self):
    cmd = "gitk &"
    aDir = self.parentAsAttribute().getEtudeWorkdirExpanded()
    if aDir != None and os.path.isdir(aDir):
      self._proc = SP.Popen(cmd, shell=True, cwd=aDir)
    else:
      desktop = self.getDesktop()
      QTW.QMessageBox.warning(desktop, "warning", "directory not existing yet:\n%s" % aDir)
      return False

  def packageLaunch(self):
    controller = self.getController()
    controller.CreateCpackEtudeSignal.emit(None) 

  def createDocLaunch(self):
    controller = self.getController()
    controller.CreateDocSignal.emit(None) 

  def updateRootlogonLaunch(self):
    controller = self.getController()
    controller.UpdateEtudeRootlogonSignal.emit(None) 

  def runPythonCode(self):
    print("runPythonCode TODO")
      
  def searchURANIEMethod(self):
    controller = self.getController()
    controller.ExecSearchMethodInURANIESignal.emit(None) 

  def printROOTContext(self):
    controller = self.getController()
    controller.ExecPrintROOTContextSignal.emit(None) 

  def postTreatments(self):
    controller = self.getController()
    controller.ExecPostTreatmentsSignal.emit(None)

  def toFileSvl(self):
    # DBG.push_debug(True)
    # 4 input files .in for solverlab launch
    aDir = self.parentAsAttribute().getEtudeWorkdirExpanded() # "/tmp"
    join = os.path.join
    case = self.caseSolverlab
    name = self.dataInformations.name

    with open(join(aDir, "Configuration.in"), "w") as f:
      toFileConfigurationIn(case, f, name)

    with open(join(aDir, "Materials.in"), "w") as f:
      toFileMaterialIn(case, f, name)

    with open(join(aDir, "Structure.in"), "w") as f:
      toFileStructureIn(case, f, name)

    with open(join(aDir, "Composition.in"), "w") as f:
      toFileCompositionIn(case, f, name)

    # DBG.pop_debug()
    return


######################################################
# write files
######################################################
def toValue(aStr, name, value):
  # In Python, strings are immutable, so override aStr outside function
  nameEg = name + "="
  if hasattr(value, "strCfg"):
    strValue = str(value.strCfg())
  elif type(value) is list:
    # str(value)[1:-1].replace(" ", "") # '[1, 2, 3]' -> '1,2,3'
    strValue = ""
    for i in value:
      if type(i) == float:
        strValue += "%.3g" % i + ","  # '1,2,3'
      else:
        strValue += str(i) + "," # '1,2,3'
        # logger.warning("toValue verify type(i) %s %s" % (type(i), i))
    strValue = strValue[0:-1]
  else:
    strValue = str(value)
  res = aStr.replace(nameEg, nameEg + strValue)
  return res


def join(*v):
  """
  as os.path.join but set antislash as slash,
  even for windows, keep windows 'c:'
  """
  # https://stackoverflow.com/questions/12086224/why-not-os-path-join-use-os-path-sep-or-os-sep
  res = os.path.join(*v).replace("\\", "/")
  return res

def toFileConfigurationIn(case, aStream, name=""):
  res = """\
# Configuration file for solverlab %s

[IonBeam]
ionZ=
ionM=
ionE0=
ion_vx=
ion_vy=
ion_vz=
ion_distribution=
enter_x=
enter_y=
enter_z=
beam_spread=

[Simulation]
max_no_ions=
ion_dose=
display_interval=
storage_interval=
status_update_interval=
store_transmitted_ions=
store_exiting_recoils=
store_exiting_limit=
store_energy_deposit=
store_ion_paths=
store_recoil_cascades=
store_path_limit=
store_path_limit_recoils=
simulation_type=
flight_length_type=
flight_length_constant=
detailed_sputtering=
min_energy=
seed1=
seed2=
OutputFileBaseName=
normalize_output=
dpa_output=
transport_type=
multiple_collisions=
scattering_calculation=
do_not_store_damage=
no_headers_in_files=
store_range3d=
store_info_file=

[Target]
straggling_model=
MaterialsFileName=
TargetstructureFileName=

""" % name
  v = case.IonBeam
  res = toValue(res, "ionZ", v.ion.atomicNumber) #=14
  res = toValue(res, "ionM", v.ion.atomicWeight) #=28.0
  res = toValue(res, "ionE0", v.ionE0*1000.) #=50000 keV to eV
  res = toValue(res, "ion_vx", v.ion_angle_xy.toVx()) #=1
  res = toValue(res, "ion_vy", v.ion_angle_xy.toVy()) #=0
  res = toValue(res, "ion_vz", 0.) #v.ion_vz) #=0 not used
  res = toValue(res, "ion_distribution", v.ion_distribution) #=0
  res = toValue(res, "enter_x", v.enter_x) #=0.
  res = toValue(res, "enter_y", v.enter_y) #=20.6
  res = toValue(res, "enter_z", v.enter_z) #=20
  res = toValue(res, "beam_spread", v.beam_spread) #=1.5

  v = case.Simulation
  res = toValue(res, "max_no_ions", v.max_no_ions)  #=2000
  res = toValue(res, "ion_dose", v.ion_dose)
  res = toValue(res, "display_interval", v.display_interval)  #=100
  res = toValue(res, "storage_interval", v.storage_interval)  #=1000
  res = toValue(res, "status_update_interval", v.status_update_interval)  #=1000
  res = toValue(res, "store_transmitted_ions", v.store_transmitted_ions)  #=1
  res = toValue(res, "store_exiting_recoils", v.store_exiting_recoils)  #=0
  res = toValue(res, "store_exiting_limit", v.store_exiting_limit)  #=100
  res = toValue(res, "store_energy_deposit", v.store_energy_deposit)  #=1
  res = toValue(res, "store_ion_paths", v.store_ion_paths)  #=0
  res = toValue(res, "store_recoil_cascades", v.store_recoil_cascades)  #=0
  res = toValue(res, "store_path_limit", v.store_path_limit)  #=100
  res = toValue(res, "store_path_limit_recoils", v.store_path_limit_recoils)  #=100
  res = toValue(res, "simulation_type", v.simulation_type)  #=0
  res = toValue(res, "flight_length_type", v.flight_length_type)  #=0
  res = toValue(res, "flight_length_constant", v.flight_length_constant)  #=0.3
  res = toValue(res, "detailed_sputtering", v.detailed_sputtering)  #=1
  res = toValue(res, "min_energy", v.min_energy)  #=5
  res = toValue(res, "seed1", v.seed1)  #=39419293
  res = toValue(res, "seed2", v.seed2)  #=93145294
  res = toValue(res, "OutputFileBaseName", join("output", "ira"))  #as files .../output/ira.*
  res = toValue(res, "normalize_output", v.normalize_output)  #=1
  res = toValue(res, "dpa_output", v.dpa_output)
  res = toValue(res, "transport_type", v.transport_type)  #=0
  res = toValue(res, "multiple_collisions", v.multiple_collisions)  #=0
  res = toValue(res, "scattering_calculation", v.scattering_calculation)  #=0
  res = toValue(res, "do_not_store_damage", v.do_not_store_damage)  #=0
  res = toValue(res, "no_headers_in_files", v.no_headers_in_files)  #=0
  res = toValue(res, "store_range3d", v.store_range3d)  #=0
  res = toValue(res, "store_info_file", v.store_info_file)  #=0

  res = toValue(res, "straggling_model", case.Target.straggling_model) #=3
  res = toValue(res, "MaterialsFileName", join(".", "Materials.in")) #=./Materials.in
  res = toValue(res, "TargetstructureFileName", join(".", "Structure.in")) #=./Structure.in
  aStream.write(res)

def toFileMaterialIn(case, aStream, name=""):
  join = os.path.join

  mat = """\
[%s]
IsVacuum=
ElementCount=
Density=
ElementsZ=
ElementsM=
ElementsConc=
ElementsDispEnergy=
ElementsLattEnergy=
ElementsSurfEnergy=
ElementsReplEnergy=
IonSurfEnergy=

"""
  matVacuum = """\
[Vacuum]
IsVacuum=1
ElementCount=0
Density=1
ElementsZ=1
ElementsM=1
ElementsConc=1
ElementsDispEnergy=1
ElementsLattEnergy=1
ElementsSurfEnergy=1
ElementsReplEnergy=-1

"""

  res = """\
# Material file for solverlab %s

""" % name
  aStream.write(res)

  for v in case.Target.Materials:
    vacuum = bool(v.IsVacuum)
    if vacuum:
     aStream.write(matVacuum)
    else:
      tmp = mat % v.name  #=GaAs
      # DBG.write("v.IsVacuum", bool(v.IsVacuum))
      tmp = toValue(tmp, "IsVacuum", v.IsVacuum)  #=0
      tmp = toValue(tmp, "ElementCount", v.getElementCount())  #=2
      tmp = toValue(tmp, "Density", v.Density)  #=4.43e22
      tmp = toValue(tmp, "ElementsSymbol", v.getElementsSymbol())  #=Ga_As
      tmp = toValue(tmp, "ElementsZ", v.getElementsZ())  #=31,33
      tmp = toValue(tmp, "ElementsM", v.getElementsM())  #=69.72,74.92
      tmp = toValue(tmp, "ElementsConc", v.getElementsConc())  #=0.5,0.5
      tmp = toValue(tmp, "ElementsDispEnergy", v.getElementsDispEnergy())  #=20.0,25.0
      tmp = toValue(tmp, "ElementsLattEnergy", v.getElementsLattEnergy())  #=3.0,3.0
      tmp = toValue(tmp, "ElementsSurfEnergy", v.getElementsSurfEnergy())  #=2.0,1.2
      tmp = toValue(tmp, "ElementsReplEnergy", v.getElementsReplEnergy())  #=2.0,1.2
      tmp = toValue(tmp, "IonSurfEnergy", v.IonSurfEnergy)  #=2.0
      aStream.write(tmp)


def toFileStructureIn(case, aStream, name=""):
  join = os.path.join
  res = """\
# Structure file for solverlab %s

[Target]
cell_count_x=
cell_count_y=
cell_count_z=
cell_size_x=
cell_size_y=
cell_size_z=
periodic_boundary_x=
periodic_boundary_y=
periodic_boundary_z=
CompositionFileType=
CompositionFileName=
UseDensityMultiplicator=
DensityMultiplicatorFileName=
special_geometry=

""" % name
  v = case.Target.Structure

  res = toValue(res, "cell_count_x", v.get_cell_count_x())  #=50
  res = toValue(res, "cell_count_y", v.cell_count_y)  #=50
  res = toValue(res, "cell_count_z", v.cell_count_z)  #=4
  res = toValue(res, "cell_size_x", v.getCellSizeX())  #=1
  res = toValue(res, "cell_size_y", v.cell_size_y)  #=1
  res = toValue(res, "cell_size_z", v.cell_size_z)  #=10
  res = toValue(res, "periodic_boundary_x", v.periodic_boundary_x)  #=0
  res = toValue(res, "periodic_boundary_y", v.periodic_boundary_y)  #=0
  res = toValue(res, "periodic_boundary_z", v.periodic_boundary_y)  #=0
  res = toValue(res, "CompositionFileType", v.CompositionFileType)  #=0
  res = toValue(res, "CompositionFileName", join(".", "Composition.in"))  #=./Composition.in
  res = toValue(res, "UseDensityMultiplicator", v.UseDensityMultiplicator)  #=0
  res = toValue(res, "DensityMultiplicatorFileName", "none")  #=none
  res = toValue(res, "special_geometry", v.special_geometry)  #=0
  aStream.write(res)

def get_value_random_multiple_materials(case):
  """
  return random value for integer indice in Materials
  using getRandomConcMaterial
  """
  concs = [float(mat.TargetConcentration) for mat in case.Target.Materials]
  concs = normalize(concs)
  # example as concs = [.5, .4, .1] for respective concentrations of 3 materials (sum as 1.)

  # cumul of concs
  irandoms = [concs[0]]
  for r in concs[1:]: # cumul
    irandoms.append(irandoms[-1] + r)
  # example as irandoms = [.5, .9, 1.] for respective concentrations of 3 materials (last as 1.)

  res = getRandomConcMaterial(irandoms)
  print("irandoms conc", concs, irandoms, res)
  return res

def getRandomConcMaterial(irandoms):
  """
  returns random indice in Materials (using TargetConcentration)
  irandoms is list as cumul of TargetConcentration
  """
  r = random.random()
  for res, ir in enumerate(irandoms):
    if r <= ir: return res
  return res

def fn_homogenous_one_material(case, aStream, options={}):
  v = case.Target.Structure
  nx = v.get_cell_count_x()
  ny = v.cell_count_y
  nz = v.cell_count_z
  for ix in range(nx):
    for iy in range(ny):
      for iz in range(nz):
        val = 0  # as first of Materials
        aStream.write("%i %i %i %i\n" % (ix, iy, iz, val))

def fn_multiLayer_multiple_materials(case, aStream, options={}):
  v = case.Target.Structure
  nx = v.cell_multiLayer_x.getCount()
  ny = v.cell_count_y
  nz = v.cell_count_z
  errors = ""
  for ix in range(nx):
    val = v.getMultiLayerMaterial(ix)
    if val >= len(case.Target.Materials):
      errors += "\nUnknown material of layer %i (as Materials[%i], fixed as Materials[0])" % (ix, val)
      val = 0
    for iy in range(ny):
      for iz in range(nz):
        aStream.write("%i %i %i %i\n" % (ix, iy, iz, val))
  if len(errors) > 0:
    logger.error("User have to fix unknown materials for '%s':%s" % (v.cell_multiLayer_x, errors))


def fn_heterogenous_random_multiple_materials(case, aStream):
  v = case.Target.Structure
  nx = v.get_cell_count_x()
  ny = v.cell_count_y
  nz = v.cell_count_z
  seed = v.seedCompositon
  #try:
  #  seed = options["seed"]
  #except:
  #  # None or no argument seeds from current time or
  #  # from an operating system specific randomness source if available
  #  seed = None
  random.seed(seed)
  for ix in range(nx):
    for iy in range(ny):
      for iz in range(nz):
        # as random in Materials (with concentrations)
        val = get_value_random_multiple_materials(case)  # as random in Materials (with concentrations)
        aStream.write("%i %i %i %i\n" % (ix, iy, iz, val))


def toFileCompositionIn(case, aStream, name):
  """see 'Creating new composition file' line 1222 in utils.C"""
  typ = str(case.Target.Structure.typeComposition)
  res = """
# Composition file for solverlab %s
# typeComposition as '%s'
""" % (name, typ)
  if typ == "homogenous one material":
    fn_homogenous_one_material(case, aStream)
    # generate error reading file ! even at last
    # aStream.write(res)
    return
  if typ == "multiLayer materials":
    fn_multiLayer_multiple_materials(case, aStream)
    # generate error reading file ! even at last
    # aStream.write(res)
    return
  """
  if typ == "heterogenous random multiple materials":
    fn_heterogenous_random_multiple_materials(case, aStream)
    # generate error reading file ! even at last
    # aStream.write(res)
    return
  """
  msg = "unexpected Target.Structure.typeComposition '%s'" % typ
  logger.error(msg)
  aStream.write("# %s\n" % msg)
  # generate error reading file ! even at last
  # aStream.write(res)
  return

def normalize(aList):
  total = float(sum(aList))
  if total == 0.:
    logger.warning("Problem total=0 in normalize from %s" % aList)
    return [float(i) for i in aList] # do nothing
  return [i/total for i in aList]

#factory pattern using xyzpy.utilsXyz._dictOfXyzClass
CLFX.appendAllXyzClasses( [
  DataInformationsSvl,
  AttributeDataFrameSvl, ListOfAttributeSvl,
  ExpressionSvl, ListOfExpressionSvl,
  DataSvl, DataManagerSvl,
  FileSvl, ListOfFileViewerXyz,
  FunctionSvl, ListOfFunctionSvl,
  LibrarySvl, ListOfLibrarySvl,
  MacroSvl, ListOfMacroSvl, MacroManagerSvl,
  UserFileSvl, ListOfUserFileSvl,
  HistoryFileManagerXyz,
  AnalysisSvl,
] )
