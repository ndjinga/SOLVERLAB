#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""
Elementaries immutables useful classes for xyzpy
"""

import os
import sys
import pprint as PP

import xyzpy.utilsXyz as UXYZ #common procedures
from PyQt5 import QtCore, QtGui, QtWidgets
import xml.etree.ElementTree as ET
from xyzpy.guiXyz.sciQDoubleSpinBox import SciQDoubleSpinBox
import xyzpy.classFactoryXyz as CLFX
from salomepy.strEvent import *
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()
_messDone = [] # to get particulars log warning message only one time

# python3 int max value
_AMinInt = -2147483648  # -99999999 old value
_AMaxInt = +2147483647  # +99999999
_AMinFloat = -1e100
_AMaxFloat = +1e100

verbose = False

###############################################################
# utilities
###############################################################

###############################################################
def getControllerInParents(aXyzItem):
  """iterate search in parents if there is a controller, else returns None"""
  parent = aXyzItem
  for nb in range(30): #precaution if infinite loop in parent()
    try:
      controller = parent.getController()
      return controller
    except:
      try:
        parent = parent.parent()
      except:
        return None
  logger.error("infinite loop calling parent()")
  return None


###############################################################
class XyzQComboBox(QtWidgets.QComboBox):
  def setValue(self, value):
    if verbose: logger.info("XyzQComboBox.setValue %s %s" % (str(value), type(value)))
    #index = self.findData(value)
    index = self.findText(str(value))
    self.setCurrentIndex(index)
    #print "setCurrentIndex",index
    
  def getValue(self):
    """
    method implemented to avoid have to use equivalent methods with other names:
    
    - QLineEdit.text()
    - QComboBox.currentText()
    - QSpinBox/QDoubleSpinBox.value()
    - etc...
    """
    return str(self.currentText())
    
###############################################################
class XyzQLineEdit(QtWidgets.QLineEdit):
  def setValue(self, value):
    if verbose: logger.info("XyzQLineEdit.setValue %s" % value)
    self.setText(value)

  def getValue(self):
    """
    method implemented to avoid have to use equivalent methods with other names:
    
    - QLineEdit.text()
    - QComboBox.currentText()
    - QSpinBox/QDoubleSpinBox.value()
    - etc...
    """
    return str(self.text())

###############################################################
class XyzQLineEditBrowseFile(QtWidgets.QLineEdit):
  """used for for asynchronous contextMenuEvent actions"""
  def setValue(self, value):
    if verbose: logger.info("XyzQLineEditBrowseFile.setValue %s" % value)
    self.setText(str(value)) #all casting

  def getValue(self):
    """
    method implemented to avoid have to use equivalent methods with other names:
    
    - QLineEdit.text()
    - QComboBox.currentText()
    - QSpinBox/QDoubleSpinBox.value()
    - etc...
    """
    return str(self.text())

  def xxevent(self, event):
    print('INFO: XyzQLineEditBrowseFile: TODO no editing event:', strEvent(event))
    return super(XyzQLineEditBrowseFile, self).event(event)
    
  def xxxsetModelItem(self, modelItem):
    """
    | store TreePyName path from root because with leaf modelItem immutables.
    | modelItem could be modified an so non permanent in menu.exec\_ actions
    | use modelItem = self._rootItem.getValueByTreePyName self._modelItemTreePyName
    """
    self._rootItem = modelItem.getRoot()
    self._modelItemTreePyName = modelItem.getTreePyName()
    #DO NOT USE, for information...
    #self._modelItem = modelItem #DO NOT USE, for information...
  
  def getModelItem(self):
    return self._rootItem.getValueByTreePyName(self._modelItemTreePyName)

  def contextMenuEvent(self, event):
    print('INFO: %s.contextMenuEvent %s' % (self.__class__.__name__, strEvent(event)))
    if hasattr(self, "_modelItemTreePyName"):
      print("there is _modelItemTreePyName: %s" % self._modelItemTreePyName)
      value = self.getModelItem()
      try: 
        actions = value.getActionsContextMenu()
      except:
        print("WARNING: There is no context menu for %s" % value._className)
        return
      menu = QtWidgets.QMenu("ContextMenu%s" % value._className, self)
      for a in actions:
        action = self._createAction(a.Name, a.Shortcut, a.ToolTip, a.Call, a.Icon, a.Enable)
        menu.addAction(action)
      menu.exec_(self.mapToGlobal(event.pos()))
    return

  def _createAction(self, Name, Shortcut, ToolTip, Call, Icon=None, Enable=True):
    # print("%s._createAction" % (self.__class__.__name__, Name))
    action = QtWidgets.QAction(Name, self)
    action.setToolTip(ToolTip)
    action.setEnabled(Enable)
    action.triggered.connect(Call)
    return action

###############################################################
class XyzQLineEditAsQLabel(QtWidgets.QLineEdit):
  """used when no editor actions needed"""
  def setValue(self, value):
    if verbose: logger.info("XyzQLineEditAsQLabel.setValue %s " % value)
    self.setText(str(value)) #all castingqeditor
    self.setToolTip("invalid editor: all modifications will be lost")

  def getValue(self):
    """
    method implemented to avoid have to use equivalent methods with other names:
    
    - QLineEdit.text()
    - QComboBox.currentText()
    - QSpinBox/QDoubleSpinBox.value()
    - etc...
    """
    return str(self.text())

  def xxevent(self, event):
    print('INFO: XyzQLineEditAsQLabel: TODO no editing event:', strEvent(event))
    return super(XyzQLineEditAsQLabel, self).event(event)
    
  def xxxsetModelItem(self, modelItem):
    """
    | store TreePyName path from root because with leaf modelItem immutables
    | modelItem could be modified an so non permanent in menu.exec\_ actions
    | use modelItem = self.getModelItem()
    """
    self._rootItem = modelItem.getRoot()
    self._modelItemTreePyName = modelItem.getTreePyName()
    #self._modelItem = modelItem #DO NOT USE, for information...
  
  def getModelItem(self):
    return self._rootItem.getValueByTreePyName(self._modelItemTreePyName)

###############################################################
class XyzQSpinBox(QtWidgets.QSpinBox):
  def setValue(self, value):
    if verbose: logger.info("XyzSpinBox.setValue %s" % value)
    super(XyzQSpinBox, self).setValue(int(value))

  def getValue(self):
    """
    method implemented to avoid have to use equivalent methods with other names:
    
    - QLineEdit.text()
    - QComboBox.currentText()
    - QSpinBox/QDoubleSpinBox.value()
    - etc...
    """
    return str(self.value())


###############################################################
class _XyzImmBase(object):
  """
  | base of previous immutables elementary class of Xyz package
  | begin with _, because have not be instancied directly
  """
  def __init__(self, *args, **kwargs):
    self._nameObject = self.__class__.__name__
    self._controller = None
    self._parentAsAttribute = None
    self._isSet= True # every times set because free attributes names
    self._nameAsAttribute = None # attribute name from parentAsAttribute
    self._defautNameAsRoot = self.__class__.__name__
    self._className = str(self.__class__.__name__)
    self._listKnownStylesXml = "withIndex,withTypeClass" #etc...
    self._styleXml = "withIndex,withTypeClass"
    self._icon = None

  def isSet(self):
    return self._isSet # theorically always True
  
  def setIsSet(self, aBool):
    """used in inherided class"""
    if not issubclass(aBool.__class__, bool): 
      raise Exception(self._className + ".setIsSet: is not a boolean")
    if aBool != True:
      raise Exception(self._className + ".setIsSet: immutable is always True")
    # self._isSet = aBool
  
  def getNameObject(self):
    return self._nameObject
  
  def getNameAsAttributeAsRoot(self):
    return self._nameObject
  
  def setNameObject(self, name):
    self._nameObject = name
  
  def sameType(self, value):
    return self.__class__(value)

  def duplicate(self):
    return self.__class__(self)

  def getValue(self):
    return self.duplicate()

  def getRoot(self):
    ii = 0
    res = self
    while True:
      ii += 1
      if ii>25: #avoid loop
        raise Exception(self._className + ".getRoot: there is a loop in tree")
      parent = res.parentAsAttribute()
      if parent == None: 
        return res
      else:
        res = parent
    pass

  def getFirstInParents(self, nameAttr):
    ii = 0
    current = self
    while True:
      ii += 1
      if ii>25: #avoid loop
        raise Exception("%s.getFirstInParents: there is a loop in tree" % self._className)
      parent = current.parentAsAttribute()
      if parent == None: return None
      try:
        res = parent.__getattribute__(nameAttr)
        return res
      except: #continue
        current = parent #and continue
      continue

  def getController(self):
    return self.getRoot()._controller
  
  def getDesktop(self):
    try:
      return self.getController().getDesktop()
    except:
      return None 

  def getControllerInParents(self, parent):
    return getControllerInParents(parent)

  def controllerRefreshViewsSignalEmit(self):
    controller = self.getController()
    if controller != None:
      controller.refreshViewsSignal.emit()
    return

  def getTreePyName(self):
    """
    returns a string "in tree path name" as "python class attribute", 
    may be used in light python meta programmation
    (do not useful use because quickly generates bugs).
    
    | example of meta programmation:
    |  >>> aTreePath = anXyzModelLeaf.getTreePyName()
    |  >>> print aTreePath
    |  >>> ".tutu.titi.toto"
    |  >>> val = "b = root.%s" % aTreePath
    |  >>> print val
    |  >>> "b = root.tutu.titi.toto"
    |  >>> res = {"root": root}
    |  >>> res = exec(val, res) # python 2-3
    |  >>> print res["b"]
    """
    ii = 0
    res = self
    resStr = ""
    while True:
      ii += 1
      if ii>25: #avoid loop
        raise Exception(self._className + ".getTreePyName: there is a loop in tree:\n" + resStr)
      parent = res.parentAsAttribute()
      if parent == None: 
        return resStr
      else:
        resStr =  res.getPyName() + resStr
        res = parent
    pass
  
  def setValueByControllerSignal(self, newValue, treePyName=None):
    """use controller signal. if there is controller."""
    tree = treePyName
    if tree == None:
      tree = self.getTreePyName()
    controller = self.getController()
    cmd = "%s = '%s'" % (tree, newValue)
    if controller != None:
      controller.setModelItemValueSignal.emit( cmd )
    else:
      print('%s.setValueByControllerSignal("%s") not done, no controller.' % (self._className, cmd))
  
  def getPyName(self):
    tmp = self.getNameAsAttribute()
    if tmp != None:
      return "." + tmp
    if self.parentAsAttribute() == None:
      return ""
    return "[%i]" % self.parentAsAttribute().index(self)

  def parentAsAttribute(self):
    return self._parentAsAttribute
  
  def getNameAsAttribute(self):
    return self._nameAsAttribute

  def toStrXml(self):
    xml = self.toXml()
    return self._prettyPrintET(xml)
    
  def isHidden(self):
    """
    | isHidden with no parameter (for tree leaf only).
    | Ask for parent tree to have context answer
    | call parent.isHidden(self.getNameAsAttribute())
    """
    parent = self._parentAsAttribute #in model
    if parent == None:
      res = False
    else:
      res = parent.isHidden(self.getNameAsAttribute())
    #print "_XyzImmBase.isHidden parent",parent._className,res
    return res

  def toXml(self, **kwargs):
    """kwarg are for optional future option of added details in xml tree"""
    #print "_XyzImmBase.toXml(%s)" % kwargs 
    res = ET.Element(self._className)
    res.text = self._toStr()
    nameAsAttribute = self.getNameAsAttribute()
    attribs = {}
    hidden = self.isHidden()
    if hidden == True: #not hidden by default
      attribs["hidden"] = str(hidden)
    
    try: #kwargs could be {} as empty
      styleXml = kwargs["styleXml"]
    except:
      styleXml = "withTypeClass,withTreePyName"
    
    if "withoutTypeClass" not in styleXml:
      attribs["typeClass"] = self.__class__.__name__  #not empty as default
    
    if nameAsAttribute != None:
      #attribs["typeClass"] = self.__class__.__name__
      res.tag = nameAsAttribute
    
    try:
      if "withoutTreePyName" not in styleXml:
        attribs["treePyName"] = self.getTreePyName()
    except:
      #logger.warning("Cannot set treePyName '%s'" % res.tag)
      print("Cannot set treePyName '%s'" % res.tag)
      attribs["treePyName"] = "?"
      pass

    res.attrib = attribs
    return res

  def _toStr(self):
    res = UXYZ.toStrForXml(self)  # formatting tagXml.text
    return res

  def createEditor(self, parent):
    """
    virtual method

    | a QComboBox example for inherited:
    | >>> combo = QtWidgets.QComboBox(parent)
    | >>> combo.addItems(["AnExemple 0", "AnExemple 1", "AnExemple etc..."])
    | >>> combo.setCurrentIndex(1)
    | >>> return combo
    """
    return None

  def setDefaultValues(self):
    """
    | virtual method which do NOTHING for immutables.
    | default values is implemented in __new__.
    | it is for convenience in tree/branch of recursive BaseXyz.setDefaultValues().
    """
    #print "%s.setDefaultValue not implemented, NEVER, for immutables" % self.__class__.__name__
    return

  def getActionsContextMenu(self):
    if verbose: print("_XyzImmBase %s.getActionsContextMenu" % str(self.__class__.__name__))
    actions = []
    parent = self.parentAsAttribute()
    isList = hasattr(parent, 'insert')
    if isList:  # parent is type list, add ContextMenu insert up/down on self which is element of list
      iindex = parent.index(self)
      tmp = self._createAction(
        'Delete', None, 'Delete me', 
         lambda status=None, i=iindex: parent.delItemSlot(status, i), 'deleteitem')
      actions.append(tmp)
      # print("isList parent ImmBase", isList, iindex)
      for iClass in parent._allowedClasses:
        name = iClass()._defautNameAsRoot #instancier pour avoir le nom str(iClass) durdur!
        tmp = self._createAction(
          'Insert up '+name, None, 'Insert item above'+name, 
          lambda status=None, i=int(iindex), c=iClass: parent.insertItemSlot(status, i, c), 'insertitemup')
        actions.append(tmp)
        tmp = self._createAction(
          'Insert down '+name, None, 'Insert item below'+name, 
          lambda status=None, i=int(iindex+1), c=iClass: parent.insertItemSlot(status, i, c), 'insertitembelow')
        actions.append(tmp)
    return actions

  def _createAction(self, Name, Shortcut, ToolTip, Call, Icon=None, Enable=True):
    # print("_XyzImmBase %s._createAction '%s'" % (self._className, Name))
    import xyzpy.baseXyz as BXYZ
    action = BXYZ.BaseXyz()
    action.Name = Name
    action.ClassName = str(self.__class__.__name__)
    action.TreePath = self.getTreePyName()
    action.Shortcut = Shortcut
    action.ToolTip = ToolTip
    action.Icon = Icon
    action.Call = Call  # Slot
    action.Enable = Enable
    return action


###############################################################
class BoolWithIndeterminatedXyz(str, _XyzImmBase):
  """
  almost immutable boolean as string that could be "True", "False" and "indeterminated"
  
  note:
  
  - used in namelist flux.cnf when QuadUse not set.
  - 0 is int(False)
  - 1 is int(True)
  - True is bool("True")
  - True is bool("False")
  - True is bool("indeterminated")
  - so use getValue() to get [True, False, None]
  - see http://stackoverflow.com/questions/2172189/why-i-cant-extend-bool-in-python
  """
  _defaultValue = "Indeterminated"
  
  def __new__(cls, value="defaultValue"):
    strValue = str(value).strip() #.replace(" ", "")
    if strValue == "defaultValue":
      val = cls._defaultValue
    elif strValue == "Indeterminated" or strValue == "None":
      val = "Indeterminated"
    elif strValue == "True" or strValue == "False":
      val = strValue
    elif strValue == ".true." or strValue == ".TRUE.": #fortran thanks
      val = "True"
    elif strValue == ".false." or strValue == ".FALSE.": #fortran thanks
      val = "False"
    else:
      raise Exception("BoolWithIndeterminatedXyz: '%s' is not in ['True', 'False', 'Indeterminated']"  % strValue)
    return str.__new__(cls, val)

  def __repr__(self):
    #__repr__ goal is to be unambiguous
    cl = self.__class__
    return cl.__name__ + '(' + self.getValue().__repr__() + ')'

  def __bool__(self): #only Python 3.x
    # print("BoolWithIndeterminatedXyz.__bool__", self)
    if self == "True": return True
    if self == "False": return False
    raise Exception("BoolWithIndeterminatedXyz: forbidden test with Indeterminated, use getValue() and assume None")

  def __nonzero__(self): #only Python 2.x
    if self == "True": return True
    if self == "False": return False
    raise Exception("BoolWithIndeterminatedXyz: forbidden test with Indeterminated, use getValue() and assume None")

  def getValue(self):
    #could be risky? wait and see
    if self == "True": return True
    if self == "False": return False
    return None

  def createEditor(self, parent):
    combo = XyzQComboBox(parent)
    combo.addItems(["True", "False", "Indeterminated"])
    combo.setValue("True")
    return combo


###############################################################
class BoolXyz(str, _XyzImmBase):
  """
  almost immutable boolean as string that could be "True", "False"
  
  note:
  
  - 0 is int(False)
  - 1 is int(True)
  - True is bool("True")
  - True is bool("False")
  - so use getValue() to get [True, False]
  - see http://stackoverflow.com/questions/2172189/why-i-cant-extend-bool-in-python
  """
  _defaultValue = "False"
  
  def __new__(cls, value=None):
    strValue = str(value).strip()# replace(" ", "")
    if value == None:
      val = cls._defaultValue
    elif strValue == "True" or strValue == "False":
      val = strValue
    elif strValue == ".true." or strValue == ".TRUE.": #fortran thanks
      val = "True"
    elif strValue == ".false." or strValue == ".FALSE.": #fortran thanks
      val = "False"
    else:
      raise Exception("BoolXyz: '%s' is not in ['True', 'False']"  % strValue)
    return str.__new__(cls, val)

  def __repr__(self):
    #__repr__ goal is to be unambiguous
    cl = self.__class__
    return cl.__name__ + '(' + self.getValue().__repr__() + ')'

  def __bool__(self): #only Python 3.x
    if self == "True": return True
    if self == "False": return False
    raise Exception("BoolXyz: forbidden value: %s" % self)

  def __nonzero__(self): #only Python 2.x
    if self == "True": return True
    if self == "False": return False
    raise Exception("BoolXyz: forbidden value: %s" % self)

  def getValue(self):
    #could be risky? wait and see
    if self == "True": return True
    if self == "False": return False
    raise Exception("BoolXyz: Problem: forbidden value: %s" % self)

  def createEditor(self, parent):
    combo = XyzQComboBox(parent)
    combo.addItems(["True", "False"])
    combo.setCurrentIndex(0)
    return combo


###############################################################
class BoolNoEditionXyz(BoolXyz):
  def createEditor(self, parent):
    return None
  

###############################################################
class NoneXyz(str, _XyzImmBase):
  """
  | almost immutable None as string "None", be only "None"
  | use getValue() to get None
  """
  _defaultValue = "None"
  
  def __new__(cls, value=None):
    val = cls._defaultValue
    return str.__new__(cls, val)

  def __repr__(self):
    #__repr__ goal is to be unambiguous
    cl = self.__class__
    return cl.__name__ + '(' + self.getValue().__repr__() + ')'

  def getValue(self):
    #could be risky? wait and see
    return None

  def createEditor(self, parent):
    combo = XyzQComboBox(parent)
    combo.addItems([sel._defaultValue])
    combo.setCurrentIndex(0)
    return combo

###############################################################
class IntXyz(int, _XyzImmBase):
  """
  int as all integers
  
  - int is immutable so you can't modify it after they are created
  - Since __init__ is called after the object is constructed,
    it is too late to modify the value for immutable types.
  - see http://stackoverflow.com/questions/2673651/inheritance-from-str-or-int
  - see https://docs.python.org/2/reference/datamodel.html
  """
  _defaultValue = 0
  
  def __new__(cls, value=None):
    if value == None:
      val = cls._defaultValue
    else:
      val = int(value)
    return int.__new__(cls, val)

  def __repr__(self):
    #__repr__ goal is to be unambiguous
    cl = self.__class__
    return cl.__name__ + '(' + int.__repr__(self) + ')'

  def __str__(self):
    res = int.__repr__(self)
    return res

  """no use
  def __str__(self):
    #__str__ goal is to be readable
    return super(self.__class__, self).__repr__()
  """
  def createEditor(self, parent):
    spin = XyzQSpinBox(parent)
    spin.setRange(_AMinInt, _AMaxInt)
    spin.setSingleStep(1)
    spin.setValue(self._defaultValue)
    return spin


###############################################################
class IntRangeXyz(IntXyz):
  """
  immutable int only in range [-99999999, +99999999]

  - restricted values only in liste define in allowedRange = [str1, str2, ...]
  - this class could to be EZ derived for modifing allowedRange
  """
  _allowedRange = [_AMinInt, _AMaxInt]

  def __new__(cls, value=None):
    miMa = cls._allowedRange
    if value == None:
      if 0 >= miMa[0] and 0 <= miMa[1]:
        val = 0
      else:
        val = miMa[0]
    else:
      val = int(value)
    cls._parentAsAttribute = None
    if val >= miMa[0] and val <= miMa[1]:
      obj = int.__new__(cls, val)
      return obj
    else:
      raise Exception("%s: '%s' is not in %s" % (cls.__name__, str(value).strip(), str(miMa)))

  def __repr__(self):
    # __repr__ goal is to be unambiguous
    cl = self.__class__
    miMa = self._allowedRange
    return "%s(%s, %s)" % (cl.__name__, int.__repr__(self), str(miMa))

  def createEditor(self, parent):
    spin = XyzQSpinBox(parent)
    miMa = self._allowedRange
    spin.setRange(miMa[0], miMa[1])
    spin.setSingleStep(1)
    spin.setValue(int(self))
    return spin


###############################################################
class Int0Xyz(IntRangeXyz):
  """ int only in [0]"""
  _allowedRange = [0, 0]

  def createEditor(self, parent):
    combo = XyzQComboBox(parent)
    combo.addItems(["0"])
    combo.setCurrentIndex(0)
    return combo

###############################################################
class Int3Xyz(IntRangeXyz):
  """int only in [3]"""
  _allowedRange = [3, 3]

  def createEditor(self, parent):
    combo = XyzQComboBox(parent)
    combo.addItems(["3"])
    combo.setCurrentIndex(0)
    return combo

###############################################################
class Int01Xyz(IntRangeXyz):
  """int only in [0,1]"""
  _allowedRange = [0, 1]
  _defaultValue = 0

  def createEditor(self, parent):
    combo = XyzQComboBox(parent)
    combo.addItems(["0", "1"])
    combo.setCurrentIndex(0)
    return combo


###############################################################
class Int02Xyz(IntRangeXyz):
  """int only in [0,1,2]"""
  _allowedRange = [0, 2]
  _defaultValue = 0

  def createEditor(self, parent):
    combo = XyzQComboBox(parent)
    combo.addItems(["0", "1", "2"])
    combo.setCurrentIndex(0)
    return combo


###############################################################
class Int12Xyz(IntRangeXyz):
  """int only in [1,2]"""
  _allowedRange = [1, 2]
  _defaultValue = 1

  def createEditor(self, parent):
    combo = XyzQComboBox(parent)
    combo.addItems(["1", "2"])
    combo.setCurrentIndex(0)
    return combo

###############################################################
class Intm11Xyz(IntRangeXyz):
  """int only in [-1, 0, 1]"""
  _allowedRange = [-1, 1]
  _defaultValue = 0

  def createEditor(self, parent):
    combo = XyzQComboBox(parent)
    combo.addItems(["-1", "0", "1"])
    combo.setCurrentIndex(0)
    return combo

###############################################################
class Int03Xyz(IntRangeXyz):
  """int only 0 to 3"""
  _allowedRange = [0, 3]
  _defaultValue = 0

  def createEditor(self, parent):
    combo = XyzQComboBox(parent)
    combo.addItems(["0", "1", "2", "3"])
    combo.setCurrentIndex(0)
    return combo

###############################################################
class Int13Xyz(IntRangeXyz):
  """int only 1 to 3"""
  _allowedRange = [1, 3]
  _defaultValue = 1

  def createEditor(self, parent):
    combo = XyzQComboBox(parent)
    combo.addItems(["1", "2", "3"])
    combo.setCurrentIndex(0)
    return combo

###############################################################
class Int05Xyz(IntXyz):
  """int only in range [0, 5]"""
  _allowedRange = [0, 5]
  _defaultValue = 0

  def sameType(self, value):
    return self.__class__(value, self.minMax)

  def createEditor(self, parent):
    combo = XyzQComboBox(parent)
    combo.addItems(["0", "1", "2", "3", "4", "5"])
    combo.setCurrentIndex(0)
    return combo

###############################################################
class IntPosXyz(IntRangeXyz):
  """int only positive"""
  _allowedRange = [0, _AMaxInt]
  _defaultValue = 0

  def createEditor(self, parent):
    spin = XyzQSpinBox(parent)
    spin.setRange(0, _AMaxInt)
    spin.setSingleStep(1)
    spin.setValue(0)
    return spin

###############################################################
class IntSupEq1Xyz(IntRangeXyz):
  """int only positive >= 1"""
  _allowedRange = [1, _AMaxInt]
  _defaultValue = 1

  def createEditor(self, parent):
    spin = XyzQSpinBox(parent)
    spin.setRange(1, _AMaxInt)
    spin.setSingleStep(1)
    spin.setValue(1)
    return spin


###############################################################
class IntSupEq2Xyz(IntRangeXyz):
  """int only positive >= 2"""
  _allowedRange = [2, _AMaxInt]
  _defaultValue = 2

  def createEditor(self, parent):
    spin = XyzQSpinBox(parent)
    spin.setRange(2, _AMaxInt)
    spin.setSingleStep(1)
    spin.setValue(2)
    return spin


###############################################################
class IntNegXyz(IntRangeXyz):
  """int only negative"""
  _allowedRange = [_AMinInt, -1]
  _defaultValue = -1

  def createEditor(self, parent):
    spin = XyzQSpinBox(parent)
    spin.setRange(_AMinInt,  -1)
    spin.setSingleStep(1)
    spin.setValue(-1)
    return spin

###############################################################
class FloatXyz(float, _XyzImmBase):
  """float"""
  _defaultValue = 0.

  def __new__(cls, value = None):
    if value == None:
      val = cls._defaultValue
    else:
      val = float(value)
    cls._parentAsAttribute = None
    return float.__new__(cls, val)

  def __repr__(self):
    #__repr__ goal is to be unambiguous
    cl = self.__class__
    return cl.__name__ + '(' + float.__repr__(self) +')'

  def __str__(self):
    return float.__repr__(self)

  def createEditor(self, parent):
    spin = SciQDoubleSpinBox(parent)
    spin.setRange(_AMinFloat, _AMaxFloat)
    spin.setDecimals(12)
    spin.setSingleStep(1)
    spin.setValue(self)
    return spin


###############################################################
class FloatRangeXyz(FloatXyz):
  """
  immutable float only in range [-1e100, +1e100]

  - restricted values only in liste define in _allowedRange = [str1, str2]
  - this class could to be EZ derived for modifing allowedRange
  """
  _defaultValue = 0.
  _allowedRange = [_AMinFloat, _AMaxFloat]

  def __new__(cls, value=None, minMax=None):
    if value == None:
      val = cls._defaultValue
    else:
      val = float(value)
    if minMax == None:
      miMa = cls._allowedRange
    else:
      miMa = [float(minMax[0]), float(minMax[1])]
    cls._parentAsAttribute = None
    if val >= miMa[0] and val <= miMa[1]:
      obj = float.__new__(cls, val)
      obj.minMax = miMa
      return obj
    else:
      raise Exception("FloatRangeXyz: '%s' is not in %s" % (str(value).strip(), str(miMa)))

  def __repr__(self):
    # __repr__ goal is to be unambiguous
    cl = self.__class__
    return cl.__name__ + '(' + float.__repr__(self) + ', ' + str(self.minMax) + ')'

  def sameType(self, value):
    return self.__class__(value, self.minMax)

  def createEditor(self, parent):
    spin = SciQDoubleSpinBox(parent)
    spin.setRange(self.minMax[0], self.minMax[1])
    spin.setSingleStep(1)
    spin.setValue(self)
    return spin


###############################################################
class FloatPosXyz(FloatRangeXyz):
  """float only positive"""
  _allowedRange = [0., _AMaxFloat]
  _defaultValue = 0.

  def createEditor(self, parent):
    spin = SciQDoubleSpinBox(parent)
    spin.setRange(0., _AMaxFloat)
    spin.setSingleStep(1)
    spin.setValue(self)
    return spin

###############################################################
class Float01Xyz(FloatRangeXyz):
  """float only in range [0., 1.]"""
  _allowedRange = [0., 1.]
  _defaultValue = 0.

  def sameType(self, value):
    return self.__class__(value)

  def createEditor(self, parent):
    spin = SciQDoubleSpinBox(parent)
    spin.setRange(self.minMax[0], self.minMax[1])
    spin.setSingleStep(0.1)
    spin.setValue(self)
    return spin


###############################################################
class StrXyz(str, _XyzImmBase):
  """
  this is immutable string, STRIPPED. # n.b. python 3 str is unicode
  """
  _defaultValue = ""
  _lenMax = None
  
  def __new__(cls, value = None):
    if value == None:
      val = cls._defaultValue
    else:
      try:
        val = str(value).strip()
      except: # from unicode to string
        val = value.encode("utf-8", "ignore") # may be latin-1 better
        val = val.strip()
    if cls._lenMax is not None:
      if cls._lenMax < len(val):
        msg = "string too long, need no more %i characters: '%s...etc'" % \
              (cls._lenMax, val[0:10])
        raise Exception(msg)
    cls._parentAsAttribute = None
    obj = str.__new__(cls, val)
    return obj
    
  def __repr__(self):
    #__repr__ goal is to be unambiguous
    cl = self.__class__
    return cl.__name__ + '(' + str.__repr__(self) + ')'

  def sameType(self, value):
    return self.__class__(value)

  def createEditor(self, parent):
    spin = XyzQLineEdit(parent)
    spin.setValue(self)
    return spin


###############################################################
class StrNoEditionXyz(StrXyz):
  def createEditor(self, parent):
    if verbose: logger.warning("%s.createEditor: no edition allowed" % self._className)
    return None


###############################################################
def getEnvvarized(value, envvars, verbose=False):
  """
  do inverse of os.path.expandvars() (on named list env vars only).
  search longest expandvars first, stop at first found replace

  | from  value '/home/johndoe/hello'
  |       envvars ["HOME", etc.]
  | returns '${HOME}/hello' (if you are johndoe of course)
  | and also for envvars environment variables
  """
  aStr = "" + value  # do the job 'cast to str' if value is StrEnvVarXyz for example
  # print("getEnvvarized", aStr, envvars)
  if "$" in aStr:  # do nothing, envvarized yet
    return aStr
  else:
    # creates [ ("/home/wambeke", "${HOME}"), (), ... ]
    envValues = [(os.path.expandvars("${%s}" % aEnv), "${%s}" % aEnv) for aEnv in envvars]

    # choose to search/replace longest expandvars first
    envValues.sort(key=lambda x: len(x[0]), reverse=True)
    # print("envValues %s" % PP.pformat(envValues))

    for expanded, aEnv in envValues:
      if expanded in aStr:
        res = aStr.replace(expanded, aEnv)
        if verbose: print("getEnvvarized %s\n  '%s' ->\n  '%s'" % (envvars, aStr, res))
        return res # stop at first found replace
  return aStr


###############################################################
class StrEnvVarXyz(StrXyz):
  """
  this is immutable string, STRIPPED. # n.b. python 3 str is unicode
  including envvarize & expandvar
  """
  _defaultValue = ""
  _lenMax = None
  _envvars = ['HOME']
  _envvarized = True # assume or not expandvar/envvarize

  def __new__(cls, value = None):
    if value == None:
      val = cls._defaultValue
    else:
      try:
        val = str(value).strip()
      except: # from unicode to string
        val = value.encode("utf-8", "ignore") # may be latin-1 better
        val = val.strip()
    if cls._lenMax is not None:
      if cls._lenMax < len(val):
        msg = "string too long, need no more %i characters: '%s...etc'" % \
              (cls._lenMax, val[0:10])
        raise Exception(msg)

    if cls._envvarized:
      valenv = getEnvvarized(val, cls._envvars)
      # print("StrEnvVarXyz %s '%s' -> '%s'" % (cls._envvars, val, valenv))
      val = valenv

    cls._parentAsAttribute = None
    obj = str.__new__(cls, val)
    return obj

  def __repr__(self):
    #__repr__ goal is to be unambiguous
    cl = self.__class__
    return cl.__name__ + '(' + str.__repr__(self) + ')'

  def __str__(self):
    if self._envvarized:
      return getEnvvarized(self, self._envvars)
    else:
      return self

  def getNameExpanded(self):
    return os.path.expandvars(str(self))

  def getNameEnvvarized(self):
    """
    from '/home/johndoe/etc' returns '${HOME}/etc' (if you are johndoe of course)
    and also for _envvars environ var
    """
    return getEnvvarized(aStr, self._envvars)

  def sameType(self, value):
    return self.__class__(value)

  def createEditor(self, parent):
    spin = XyzQLineEdit(parent)
    spin.setValue(self)
    return spin


###############################################################
class FileXyz(StrEnvVarXyz):
  """
  createEditor as qdialog get file name
  accept $ENVVAR in self._envvars []
  """
  _envvars = ['HOME']
  _title = 'Select file'
  _filter = 'any files(*)' #for classical QFileDialog browse
  #_filter example: _filter = 'element files (*.xml *.XML)\nany files(*)'
  _directory = None # default directory for QFileDialog opening
  _envvarized = True # assume or not expandvar/envvarize
   
  def getActionsContextMenu(self):
    # actions are list of BaseXyz(), see _XyzImmBase._createAction
    actions = super(FileXyz, self).getActionsContextMenu()
    actions.append( self._createAction("Browse file", None, "Browse", self.browseDialog, "browsefile") )
    name = str(self)
    _, ext = os.path.splitext(name)
    if ext in ".med".split():
      b = self._createAction("Help med API", None, "help() for MEDLoader", self.quickHelpMed, "editor")
      actions.append(b)
      a = self._createAction("Dump asci contents", None, "Dump asci contents with mdump", self.quickEdit, "editor")
      actions.append(a)
    else:
      a = self._createAction("Quick asci edit", None, "Edit asci content", self.quickEdit, "editor")
      actions.append(a)
    # disable/enable last action, as not pertinent
    # TODO mimetypes test ?
    a.Enable = os.path.isfile(name)
    return actions

  def quickEdit(self):
    aFile = str(self)
    controller = self.getController()
    if controller is None:
      logger.warning("no controller for action quickEdit on %" % self._className)
    else:
      try:
        controller.centralLogView.centralWidget().quickEditFiles([aFile])
      except Exception as e:  # sometimes (in test) ControllerXyz object has no attribute centralLogView
        logger.warning("problem for action quickEdit : %s" % e)

  def quickHelpMed(self):
    import helppy.helpPager as HP
    import MEDLoader
    # aFile = HP.helpToFile(MEDLoader, "/tmp/${USER}/helpPager.tmp")
    aFile = HP.helpToFile(MEDLoader.MEDFileData(), "/tmp/${USER}/helpPager.tmp")
    controller = self.getController()
    if controller is None:
      logger.warning("no controller for action quickHelpMed on %" % self._className)
    else:
      try:
        controller.centralLogView.centralWidget().quickEditFiles([aFile])
      except Exception as e:  # sometimes (in test) ControllerXyz object has no attribute centralLogView
        logger.warning("problem for action quickEdit : %s" % e)

  def browseDialog(self):
    
    #TODO decide... import salomepy.xsalomesession as XSS
    #if verbose: print "%s.browseDialog" % self._className
    #TODO decide... desktop = XSS.getDesktop() #TODO get sender and treeView...
    desktop = self.getDesktop()
    aDialog = QtWidgets.QFileDialog()
    urlsIni = [i for i in aDialog.sidebarUrls()]
    
    urls = [i for i in urlsIni]
    for var in self._envvars:
      aDir = os.getenv(var)
      if aDir != None:
        oneUrl = QtCore.QUrl.fromLocalFile(aDir)
        if oneUrl not in urls: urls.append(oneUrl)

    aFileCurrent = ""
    aFileOrDir = str(self) # user may be set an incorrect string.
    if aFileOrDir != "":
      aFileCurrent = aFileOrDir
      aFileCurrent = os.path.realpath(os.path.expandvars(aFileCurrent))
      oneUrl = QtCore.QUrl.fromLocalFile(os.path.dirname(aFileCurrent))
      if oneUrl not in urls: urls.append(oneUrl)
    
    aDialog.setSidebarUrls(urls)

    if aFileCurrent == "":
      aFileCurrent = self._directory

    # print("aFileCurrent", aFileCurrent)
    choose = aDialog.getOpenFileName(desktop, self._title, str(aFileCurrent), self._filter)
    # print("choose file", choose) # is tuple
    choose = str(choose[0])
    if choose == "": return True #cancel
    realPath = os.path.realpath(choose)
    if not os.path.isfile(realPath):
      QtWidgets.QMessageBox.warning(desktop, "warning", "Needs to be a file\n'%s'" % realPath)
      return False

    if self._envvarized:
      realPath = getEnvvarized(realPath, self._envvars)

    #do refresh views if controller
    controller = self.getController() #_controllerForActions
    if controller == None:
      print("ERROR: no controller for action browse on",self._className)
    else:

      #controller.setModelItemValueSignal.emit( "%s = '%s'" % (self._treePyNameForActions, realPath) )
      #controller.setModelItemValueSignal.emit( "%s = '%s'" % (self.getTreePyName(), realPath) )
      
      #stuff set type (no casting) if appended in ListOf...
      cmd = "%s = args[1]" % self.getTreePyName()
      arg1 = self.__class__(realPath)
      controller.setModelItemValueSignalList.emit( [ cmd, arg1 ] )
    return True

###############################################################
class FileOnlyBrowseXyz(FileXyz):
  """
  as FileXyz, only menu browse, without direct str edition
  """
  def createEditor(self, parent):
   return None

###############################################################
class FileViewerXyz(StrEnvVarXyz):
  """
  createEditor as qdialog get file name
  accept $ENVVAR in self._envvars []
  """
  _envvars = ['HOME']
  _title = 'Select file'
  _filter = 'any files(*)' #for classical QFileDialog browse
  #_filter example: _filter = 'element files (*.xml *.XML)\nany files(*)'
  _typesFiles = "*".split()
  _directory = "$HOME"
  _envvarized = False # assume or not expandvar/envvarize

  def createEditor(self, parent): #no direct edition: use 'Browse'
    return None
    
  def getActionsContextMenu(self):
    #print "+++++ %s.getContextMenu" % self._className
    actions = super(FileViewerXyz, self).getActionsContextMenu()
    actions.append( self._createAction('BrowseViewer', None, 'Browse with Viewer', self.browseViewerDialog, 'browsefile') )
    #TODO launch ROOT macro .C
    _, ext = os.path.splitext(self)
    if ext in [".C", ".py"]: 
      actions.append( self._createAction('Process', None, 'ROOT Process file', self.rootProcessFile, 'run') )
    return actions
  
  def browseViewerDialog(self):
    #print "%s.browseViewerDialog" % self._className
    desktop = self.getDesktop()
    controller = self.getController()
          
    if controller == None:
      QtWidgets.QMessageBox.warning(desktop, "warning", "%s: needs a controller for browseViewer" % self._title)
      return False

    aDialog = controller.getExploreDir()
    strForWhat ="Select ONE file"
    aDialog.lockSetDirRootPath(self._directory, filters=self._typesFiles)
    controller.showExploreDir()
    aDialog.lockForChangeModelController(strForWhat, self.browseViewerExecOnApply)
    return True

  def browseViewerExecOnApply(self, selectedFiles):
    """part of Apply on browseViewerDialog"""
    #print "%s.browseViewerExecOnApply" % self._className, selectedFiles
    if len(selectedFiles) > 1: logger.warning("get only first file of multiple selection files")
    if len(selectedFiles) == 0: 
      logger.warning("no file selected: as 'Cancel'")
      return False
    nameFile = str(selectedFiles[0])
    if nameFile == "": return True #cancel
    controller = self.getController()
    if controller == None:
      logger.error("no controller for %s action 'Apply'" % self._className)
    else:
      cmd = "%s = args[1]" % self.getTreePyName()
      arg1 = self.__class__(nameFile) #stuff set type (no casting) if appended in ListOf...
      controller.setModelItemValueSignalList.emit( [ cmd, arg1 ] )
    return True

  def rootProcessFile(self):
    controller = self.getController()
    controller.ExecRootProcessFileSignal.emit(os.path.expandvars(str(self)))
    return True  

###############################################################
class DirectoryXyz(StrEnvVarXyz):
  """
  createEditor as qdialog get file name
  accept $ENVVAR in self._envvars []
  """
  _envvars = ['HOME', 'DEFAULT_WORKDIR']
  _title = 'Select directory'
  # _filter = 'any directory(*)'  # for classical QFileDialog browse
  _directory = None  # default directory for QFileDialog opening
  _envvarized = True  # assume or not expandvar/envvarize

  def getActionsContextMenu(self):
    actions = super(DirectoryXyz, self).getActionsContextMenu()
    actions.append(self._createAction('Browse', None, 'Browse', self.browseDialog, 'browsefile'))
    return actions


  def browseDialog(self):

    # TODO decide... import salomepy.xsalomesession as XSS
    # if verbose: print "%s.browseDialog" % self._className
    # TODO decide... desktop = XSS.getDesktop() #TODO get sender and treeView...
    desktop = self.getDesktop()
    aDialog = QtWidgets.QFileDialog()
    aDialog.setFileMode(aDialog.DirectoryOnly)
    # aDialog.setOption(aDialog.ShowDirsOnly)
    # print("aDialog.FileMode", dir(aDialog.FileMode))
    urlsIni = [i for i in aDialog.sidebarUrls()]

    urls = [i for i in urlsIni]
    for var in self._envvars:
      aDir = os.getenv(var)
      if aDir != None:
        oneUrl = QtCore.QUrl.fromLocalFile(aDir)
        if oneUrl not in urls: urls.append(oneUrl)

    aDirCurrent = ""
    aFileOrDir = str(self) # user may be set an incorrect string.
    if aFileOrDir != "":
      if os.path.isdir(aFileOrDir):
        aDirCurrent = aFileOrDir
      else:
        aDirCurrent, _ = os.path.split(aFileOrDir)
      oneUrl = QtCore.QUrl.fromLocalFile(aDirCurrent)
      if oneUrl not in urls: urls.append(oneUrl)

    aDialog.setSidebarUrls(urls)

    if aDirCurrent == "":
      aDirCurrent = self._directory

    choose = aDialog.getExistingDirectory(desktop, self._title, str(aDirCurrent))
    # print("choose dir", choose)
    choose = str(choose)
    if choose == "": return True  # cancel
    realPath = os.path.realpath(choose)
    if not os.path.isdir(realPath):
      QtWidgets.QMessageBox.warning(desktop, "warning", "Needs to be a directory\n'%s'" % realPath)
      return False

    if self._envvarized:
      realPath = getEnvvarized(realPath, self._envvars)

    # do refresh views if controller
    controller = self.getController()  # _controllerForActions
    if controller == None:
      print("ERROR: no controller for action browse on", self._className)
    else:

      # controller.setModelItemValueSignal.emit( "%s = '%s'" % (self._treePyNameForActions, realPath) )
      # controller.setModelItemValueSignal.emit( "%s = '%s'" % (self.getTreePyName(), realPath) )

      # stuff set type (no casting) if appended in ListOf...
      cmd = "%s = args[1]" % self.getTreePyName()
      arg1 = self.__class__(realPath)
      controller.setModelItemValueSignalList.emit([cmd, arg1])
    return True


###############################################################
class DateXyz(StrXyz):
  """
  TODO better implementation.
  accept or not string with
  createEditor as qdialog get directory name,
  accept $ENVVAR?
  """
  pass


###############################################################
class StrInListXyz(str, _XyzImmBase):
  """
  immutable string, base class, NOT FOR DIRECT USE:
  
  | - restricted values only in liste define in \_allowedList = [str1, str2, ...]
  | - this class have to be derived for use
  |
  | example:
  | >>> class MyFruits(StrInListXyz):
  | >>>   _allowedList = ['apple', 'blueberry']
  | >>>   _defaultValue = 'blueberry' #by default is _allowedList[0]
  | >>>   pass
  """
  _allowedList = ["StrInListXyz class have to be derived for use"]
  _defaultValue = None
  
  def __new__(cls, value=None):
    allowList = cls._allowedList
    if value == None:
      if cls._defaultValue == None:
        val = allowList[0]
      else:
        val = cls._defaultValue
    else:
      val = str(value).strip()
    cls._parentAsAttribute = None
    if val in allowList :
      obj = str.__new__(cls, val)
      return obj
    else:
      raise Exception("%s: '%s' is not in:\n%s" % (cls.__name__, str(value).strip(), PP.pformat(allowList)))

  def __repr__(self):
    #__repr__ goal is to be unambiguous
    cl = self.__class__
    return "%s(%s, %s)" % ( cl.__name__, str.__repr__(self), str(self._allowedList) )

  def sameType(self, value):
    return self.__class__(value)

  def createEditor(self, parent):
    combo = XyzQComboBox(parent)
    combo.addItems(self._allowedList)
    combo.setCurrentIndex(0)
    return combo


###############################################################
class StrModelXyz(StrInListXyz):
   _allowedList = ["FeCu", "FeCu_p", "FeHe"]
   pass

###############################################################
class StrTypeProgXyz(StrInListXyz):
   _allowedList = ["g", "c"]
   pass

###############################################################
class NameInGeomShapeXyz(StrXyz):
  """
  immutable string, base class for inherited. 
  
  | - find in salome browser all existing geom shapes
  | - createEditor with a QComboBox of shape names
  |
  | example:
  | >>> class MyGrainShapes(NameInGeomShapeXyz):
  | >>>   _searchStr = 'grain_'
  | >>>   pass
  """
  _searchStr = "" #will search string in geom name shapes
  
  def __init__(self, *args, **kwargs):
    super(NameInGeomShapeXyz, self).__init__(*args, **kwargs)
    self._defautNameAsRoot = "NameInGeomShape"

  """find in salome browser all existing geom '*_searchStr*' shapes"""
  def createEditor(self, parent):
    """will find in salome browser all existing geom self._searchStr shapes"""
    return self.createEditorInGeom(parent)

  def createEditorInGeom(self, parent):
    import salomepy.xsalomesession as XSS
    res = XSS.getGeomSobj(verbose)
    XSS.getOrLoadSobj(res)
    items = XSS.getSobjPath(res)
    if verbose: logger.info("createEditorInGeom %s" % inStr)
    if items == None:
      mess = "No Active Salome Geometry, problem for ComboBox of Geom objets names"
      if mess not in _messDone:
        logger.error(mess)
        UXYZ.QMessageBoxWarning(parent, mess)
        _messDone.append(mess)
      #QtWidgets.QMessageBox.warning(parent, "warning", mess)
      return None
    if verbose: logger.info("items from geom %s" % items)
    inStr = self._searchStr
    if inStr != "":
      itemsBox = [item.replace("/Geometry/","") for item in items if inStr in item]
    else:
      itemsBox = [item.replace("/Geometry/","") for item in items]
    itemsBox = [item for item in itemsBox if item != ""] #erase for Ref.Entry objects... indirection
    if len(itemsBox) == 0:
      #create one with CreateShapeBox action for example
      mess = "No item name with search string '%s' in Geometry Browser" % inStr
      UXYZ.QMessageBoxWarning(parent, mess)
      return None
    combo = XyzQComboBox(parent)
    combo.addItems(itemsBox)
    combo.setValue(itemsBox[0])
    return combo


###############################################################
def checkRelease(value):
  def checkItem(val, value):
    try:
      res = int(val)
      if res < 0 or res > 99: raise Exception("incorrect release string '%s'" % value)
      return str(res)
    except:
      raise Exception("incorrect release string '%s'" % value)

  try:
    a = value.split(".")
    res = [checkItem(i, value) for i in a]
    if len(res) == 0: raise Exception("incorrect release string '%s'" % value)
    if len(res) > 3: raise Exception("incorrect release string '%s'" % value)
    if len(res) == 1: res = res[0]+".0.0"
    if len(res) == 2: res = res[0]+"."+res[1]+".0"
    if len(res) == 3: res = res[0]+"."+res[1]+"."+res[2]
    return res
  except:
    raise Exception("incorrect release string '%s'" % value)

###############################################################
class ReleaseXyz(StrNoEditionXyz):
  _defaultVersion = "0.0.0"
  
  def __new__(cls, value = None):
    if value == None:
      val = cls._defaultVersion
    else:
      val = checkRelease(value) #str(value).strip(".")
    cls._parentAsAttribute = None
    obj = str.__new__(cls, val)
    return obj

  def checkRelease(self, value):
    return checkRelease(value)

  def getActionsContextMenu(self):
    """no browse"""
    if verbose: logger.info("%s.getContextMenu" % self._className)
    actions = []
    actions.append( self._createAction('Modify', None, 'set other release', self.createEditorData, 'modify') )
    return actions

  def createEditor(self, parent):
    # if verbose: print("%s.createEditor: no direct edition allowed" % self._className)
    return None
  
  def createEditorData(self, parent=None):
    if verbose: logger.info("%s.createEditorData parent %s" % (self._className, parent))
    
    controller = self.getController()
    aDialog = self.getEditDialog()
    aDialog.setModal(True)
    aDialog.closeOnApply = True
    aDialog.exec_()
    if aDialog.choice == "Cancel": return True
    
    res = aDialog.getLocalModel()
    newValue = "%i.%i.%i" %(res.major, res.minor, res.patch)
    if verbose: 
      logger.info("%s.createEditorData: new value: '%s'" % (self._className, newValue))

    if controller == None:
      logger.error("%s.createEditorData: no controller for action View/Edit" % self._className)
    else:
      controller.setModelItemValueSignal.emit( "%s = '%s'" % (self.getTreePyName(), newValue) )
      #controller.setModelItemValueSignalList.emit( ["%s.setCurrentData(args[1])" % self.getTreePyName(), newValue ] )
    return True

  def getEditDialog(self):
    import xyzpy.guiXyz.dialogXmlXyz as DXYZ
    import xyzpy.majorMinorPatchXyz as MMP
    if verbose: logger.info("%s.editDialog" % self._className)

    a = MMP.MajorMinorPatchXyz()
    try:
      major, minor, patch = self.split(".")
    except:
      major, minor, patch = self._defaultVersion.split(".")
    a.major = major
    a.minor = minor
    a.patch = patch

    widDialog = DXYZ.DialogXmlXyz()
    widDialog.setFromXml(a.toStrXml(), modeView="allInTab")
    return widDialog


###############################################################
class ExpressionXyz(StrNoEditionXyz):
  _defaultValue = "NoName=NoValue"
  
  def getActionsContextMenu(self):
    """no browse"""
    if verbose: logger.info("%s.getContextMenu" % self._className)
    actions = [ self._createAction('Modify', None, 'set other release', self.createEditorData, 'modify') ]
    return actions

  def createEditor(self, parent):
    # if verbose: print("%s.createEditor: no direct edition allowed" % self._className)
    return None
  
  def createEditorData(self, parent=None):
    print("%s.createEditorData parent %s" % (self._className, parent))
    controller = self.getController()
    aDialog = self.getEditDialog()
    aDialog.setModal(True)
    aDialog.closeOnApply = True
    aDialog.exec_()
    if aDialog.choice == "Cancel": return True
    
    res = aDialog.getLocalModel()
    newValue = res.getValue()
    if verbose: 
      logger.info("%s.createEditorData: new value: '%s'" % (self._className, newValue))

    if controller == None:
      logger.error("%s.createEditorData: no controller for action View/Edit" % self._className)
    else:
      #print "******** ExpressionXyz: %s = '%s'" % (self.getTreePyName(), newValue)
      #controller.setModelItemValueSignal.emit( "%s = '%s'" % (self.getTreePyName(), newValue) )
      controller.setModelItemValueSignalList.emit( ["%s = args[1]" % self.getTreePyName(), self.__class__(newValue) ] )
    return True

  def getEditDialog(self):
    import xyzpy.guiXyz.dialogXmlXyz as DXYZ
    import xyzpy.nameValueXyz as NV
    if verbose: logger.info("%s.editDialog" % self._className)

    a = NV.NameValueXyz()
    try:
      name, value = self.split("=")
    except:
      name, value = self._defaultValue.split("=")
    a.name = name
    a.value = value

    widDialog = DXYZ.DialogXmlXyz()
    widDialog.setFromXml(a.toStrXml(), modeView="allInTab")
    return widDialog


#print "dir!",dir()
#print "locals!",locals()
#print "__builtins__", dir(__builtins__)
#listOfXyzLocals =[]
#listOfXyzLocals = [value for key, value in locals().items() if key[-3:] == "Xyz"]
#print listOfLocalsXyz
#keys = [key for key in locals().keys() if key[-3:] == "Xyz"]
#CLFX.appendAllXyzClasses( [locals()[key] for key in keys] )


#factory pattern using xyzpy.utilsXyz._dictOfXyzClass
CLFX.appendAllXyzClasses( locals() )
