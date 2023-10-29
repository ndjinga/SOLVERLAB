#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


__all__ = ['BaseXyz', 'BaseFreeXyz', 'BaseFreeAllXyz', 'ListOfBaseXyz', '_XyzConstrainBase']

import traceback
import json
import pprint as PP
import xml.etree.ElementTree as ET
from functools import partial

import xyzpy.loggingXyz as LOG
import xyzpy.utilsXyz as UXYZ
import xyzpy.classFactoryXyz as CLFX
# import debogpy.debug as DBG

logger = LOG.getLogger()

verbose = False
verboseEvent = verbose
debug = verbose

_messSome = []

########################################################################################
# json utils
########################################################################################
def dumper(obj):
    """goal is to json explore subclass object as dict"""
    return obj.__dict__

def dumperType(obj):
    """goal is to get a "_type" trace json subclass object, but ignore attributes begining with '_'"""
    typeatt = "_type"
    aDict = dict((k,v) for k, v in list(obj.__dict__.items()) if k[0] != "_" or k == typeatt)
    if typeatt not in aDict: aDict[typeatt] = obj.__class__.__name__
    return aDict

def jsonDumps(obj):
    """to get direct default jsonDumps method"""
    return json.dumps(obj, default=dumperType, sort_keys=True, indent=2)

########################################################################################
class BaseXyz(object):
  """
  goal is as a simple dictionary with keys
  a = BaseXyz()
  a.tintin = "reporter"
  a.milou = "dog"
  print "tintin:",a.tintin
  ...as...
  a = {}
  a["tintin"] = "reporter"
  a["milou"] = "dog"
  print "tintin:",a["tintin"]
  """

  def __repr__asList(self):
    #__repr__ goal is to be unambiguous
    aList = [] #an ordered list representation is better for test (and visualize) (in)equality
    for k in sorted(self.__dict__.keys()):
      if k[0] != '_':
        aList.append( [k, self.__dict__[k]]  )
    return self.__class__.__name__ + " = " + aList.__repr__()

  def setAttributesFromDict(self, aDict):
    """no deep copy of values"""
    for key, value in aDict.items():
      self.__setattr__(key, value)
    # print("setAttributesFromDict", self)

  def __repr__(self):
    #__repr__ goal is to be unambiguous, easy human readeable
    return self._reprIndent()

  def _reprIndent(self, indent=""):
    res = ""
    newIndent = indent + "   "
    for k in sorted(self.__dict__.keys()):
      if k[0] != '_':
        kk = self.__dict__[k]
        if issubclass(BaseXyz, kk.__class__):
          res += "\n" + newIndent + "[%s, " % k + kk._reprIndent(newIndent) + "]"
        else:
          res += "\n" + newIndent + "%s" % [k, kk]
    return self.__class__.__name__ + " = " + res

  def jsonDumps(self):
    return jsonDumps(self)



########################################################################################
class _BaseCommonXyz(object):
  """
  common methods for BaseXyz ListOfBaseXyz
  """
  _lf = "\n"
  _icon = None

  def getAppendingToolTip(self):
    """return tooltip to append in tooltip of parent in tree"""
    return None

  def getNameAsAttribute(self):
    return self._nameAsAttribute

  def getNameAsAttributeAsRoot(self):
    return self._defautNameAsRoot

  def getNameObject(self):
    return self._nameObject

  def objectName(self):
    try:
      return self._nameObject
    except: #if _parentAsAttribute not set yet, (during initialisation)
      #print "_nameObject not set yet (during initialisation)?"
      return "Unknown (during initialisation)?"

  def setNameObject(self, name):
    self._nameObject = name

  def getAttributeType(self, nameAttribute):
    return self._attributesDict[nameAttribute]

  def sameType(self, value):
    return self.__class__(value)

  def parentAsAttribute(self):
    try:
      return self._parentAsAttribute
    except: #if _parentAsAttribute not set yet, (during initialisation)
      #print "_parentAsAttribute not set yet (during initialisation)?"
      return None

  def getRoot(self):
    ii = 0
    res = self
    while True:
      ii += 1
      if ii > 25:  # avoid loop
        raise Exception("%s.getRoot: there is a loop in tree" % self._className)
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
      if ii > 25:  # avoid loop
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

  def setController(self, controller):
    root = self.getRoot()
    if root._controller == None:
      root._controller = controller
    else:
      raise Exception("%s.setController : controller set yet for" % self._className)

  def controllerRefreshViewsSignalEmit(self):
    controller = self.getController()
    if controller != None:
      controller.refreshViewsSignal.emit()
    return

  def getTreePyName(self):
    """
    used in light python meta programmation

    | example:
    | >>> val = 'b=root.tutu.titi.toto'
    | >>> namespace = {"root": root}
    | >>> res = exec(val, namespace) # python 2-3
    | >>> b = res["b"]
    """
    ii = 0
    res = self
    resStr = ""
    while True:
      ii += 1
      if ii > 25:  # avoid loop
        raise Exception("%s.getTreePyName: there is a loop in tree:\n%s" % (self._className, resStr))
      parent = res.parentAsAttribute()
      if parent == None:
        return resStr
      else:
        resStr =  res.getPyName() + resStr
        res = parent
    pass

  def setValueByTreePyName(self, treePyName, newValue):
    if treePyName == '':
      logger.error("%s.setValueByTreePyName of root forbidden (as change self)" % self._className)
      return
    cmd = 'myself%s = newValue' % treePyName
    if verbose: logger.info("%s.setValueByTreePyName exec('%s')" % (self._className, val))
    namespace = {"myself": self, "newValue": newValue}
    namespace = UXYZ.exec2or3(cmd, namespace)  # TODO python 2-3
    return

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
      logger.error('%s.setValueByControllerSignal("%s") not done because no controller.' % (self._className, cmd))

  def getValueByTreePyName(self, treePyName):
    """
    you cannot change local variables in function scope in Python 3 using exec,
    even though it was possible in Python 2. Not even previously declared variables.
    In Python 2, using the exec statement meant the compiler knew to switch off the local scope optimizations
    With exec() being a function, that option is no longer available and function scopes are now always optimized.
    Moreover, in Python 2, the exec statement explicitly copies all variables found in locals()
    back to the function locals using PyFrame_LocalsToFast,
    but only if no globals and locals parameters were supplied.
    The proper work-around is to use a new namespace (a dictionary) for your exec()
    """
    namespace = {"myself": self}
    cmd = "res = myself%s" % treePyName
    namespace = UXYZ.exec2or3(cmd, namespace)  # TODO python 2-3
    res = namespace["res"]
    if verbose:
      logger.info("%s.getValueByTreePyName exec('%s') = %s" % (self._className, cmd, type(res)))
    return res

  def getPyName(self):
    tmp = self.getNameAsAttribute()
    if tmp != None:
      return "." + tmp
    if self.parentAsAttribute() == None:
      return ""
    #case parent is listofbase
    return "[%i]" % self.parentAsAttribute().index(self)

  def getDictOfTreePyNameStrValues(self, aDict=None):
    """
    to get easy representation of values in self
    """
    if aDict == None:
      res = {}
    else:
      res = aDict
    for key, value in list(self.items()):
      if issubclass(value.__class__, BaseFreeXyz) or \
         issubclass(value.__class__, ListOfBaseXyz):
        value.getDictOfTreePyNameStrValues(res)
      else:
        treePyName = value.getTreePyName()
        res[treePyName] = str(value)
    return res

  def setValueFromDictOfTreePyNameStrValues(self, aDict):
    for treePyName, strValue in list(aDict.items()):
      namespace = {"myself": self, "strValue": strValue}
      cmd = 'myself%s = myself%s.sameType(strValue)' % (treePyName, treePyName)
      try:
        namespace = UXYZ.exec2or3(cmd, namespace)  # TODO python 2-3
      except Exception as e:
        msg = "problem model setValue exec('%s'):\n%s\n" % (cmd, e)
        logger.critical(msg)
        raise Exception(msg)
    return

  def isSet(self):
    return self._isSet

  def setIsSet(self, aBool):
    """used in inherided class"""
    if not issubclass(aBool.__class__, bool):
      raise Exception("%s.setIsSet parameter is not a boolean" % self._className)
    self._isSet = aBool

  def isHidden(self, nameAttr):
    """virtual method to know if attribute is currently displayed in treeView and other dialog widget"""
    return False

  def toStrXml(self, **kwargs):
    xml = self.toXml(**kwargs)
    return UXYZ.prettyPrintET(xml)

  def toFileXml(self, fileName, **kwargs):
    strXml = self.toStrXml(**kwargs)
    file = open(fileName, "w")
    file.write(strXml)
    file.close()

  def setFromFileXml(self, fileName, verbose=False):
    raise Exception("%s.setFromFileXml: use 'aData1 = UXYZ.fromFileXml(testFile1)' instead" % self._className)

  def setDefaultStyleXml(self, kwargs):
    if "StyleXml" not in kwargs:
      kwargs["StyleXml"]=self._styleXml

  def toXml(self, **kwargs):
    """kwarg are for optional future option of added details in xml tree"""
    return UXYZ.toXml(self, **kwargs)

  def fromXml(self, tagXmlOrStrXml, **kwargs):
    """kwarg are for optional future option of added details in xml tree"""
    raise Exception("fromXm: virtual method, not set for class '%s'" % self.__class__.__name__)
    tagXml = None
    if issubclass(tagXmlOrStrXml.__class__, str):
      tagXml = ET.fromstring(tagXmlOrStrXml)
    else:
      tagXml = tagXmlOrStrXml
    if verbose:
      logger.info("%s.fromXml:\n'%s' '%s' '%s'" % (self._className, tagXml.tag, tagXml.text, str(tagXml.attrib)))
    #return

  def duplicate(self):
    """with xml format deep copy, cpu expensive copy, controller not set"""
    if verbose: logger.info("%s.duplicate" % self._className)
    return UXYZ.fromXml(self.toXml())

  def setDefaultValues(self):
    """virtual"""
    logger.warning("%s.setDefaultValues method is virtual, as nothing done." % self._className)

  def getActionsContextMenu(self):
    # print("_BaseCommonXyz %s.getActionsContextMenu" % self._className)
    actions = []
    try:
      import iradinapy.configIra as CFGIRA
      mode = CFGIRA.getCurrentMode()
    except:
      mode = "advanced"
      pass
    if mode == "advanced":
      tmp = self._createAction('Edit dialog', None, 'Edit dialog widget', self.editDialog, 'dialogwidget')
      actions.append(tmp)
    tmp = self._createAction('Delete', None, 'Delete me', self.delItemAsMe, 'deleteitem')
    actions.append(tmp)
    parent = self.parentAsAttribute()
    isList = hasattr(parent,'insert')
    if isList: # parent is type list, add ContextMenu insert up/down on self which is element of list
      iindex = parent.index(self)
      # print("isList parent BaseCommon", isList, iindex)
      for iClass in parent._allowedClasses:
        cClass = CLFX.toClass(iClass)
        name = CLFX.toString(cClass) # instancier pour avoir le nom str(iClass) durdur!
        tmp = self._createAction('Insert up '+name, None, 'Insert item above'+name, lambda status=None, i=int(iindex), c=cClass: parent.insertItemSlot(status, i, c), 'insertitemup')
        actions.append(tmp)
        tmp = self._createAction('Insert down '+name, None, 'Insert item below'+name, lambda status=None, i=int(iindex+1), c=cClass: parent.insertItemSlot(status, i, c), 'insertitembelow')
        actions.append(tmp)
    tmp = self._createAction('Delete child items', None, 'Delete all child items', self.delAllItems, 'deletechilditem')
    actions.append(tmp)
    return actions

  def _createAction(self, Name, Shortcut, ToolTip, Call, Icon=None, Enable=True):
    # print("_BaseCommonXyz._createAction %s" % Name)
    action = BaseXyz()
    action.Name = Name
    action.ClassName = str(self.__class__.__name__)
    action.TreePath = self.getTreePyName()
    action.Shortcut = Shortcut
    action.ToolTip = ToolTip
    action.Icon = Icon
    action.Call = Call  # Slot
    action.Enable = Enable
    return action

  def editDialog(self):
    import xyzpy.guiXyz.dialogXmlXyz as DXYZ
    logger.debug("%s.editDialog standard by default" % self._className)
    widDialog = DXYZ.DialogXmlXyz()
    #widDialog.resize(500, 500)
    #widDialog.setFromXml(self.toStrXml(modeView="allInTab")) #(self.xmlData)
    controller = self.getController()
    if controller == None:
      logger.warning("%s.editDialog with no controller" % self._className)
      self._widDialog = widDialog
      widDialog.setFromXml(self.toStrXml(), modeView="allInTab") #(self.xmlData)
      widDialog.show()
    else:
      treePyName = self.getTreePyName() #for correspondance in model of controller
      controller.assumeDialogBox(widDialog, treePyName)

      #widDialog.setFromXml(self.toStrXml(), modeView="allInTab", treePyNameRoot=treePyName) #(self.xmlData)

      #TODO send all data, not subtree
      if debug: logger.info("%s.editDialog send all data, not subtree" % self._className)
      dataStrXml = self.getRoot().toStrXml()
      widDialog.setFromXml(dataStrXml, modeView="allInTab", treePyNameRoot=treePyName) #(self.xmlData)
    return

  def delItemAsMe(self):
    # DBG.write("%s.delItemAsMe with no parameters" % self._className, " ")
    myTreePyName = self.getTreePyName()
    root = self.getRoot()
    controller = self.getController()
    if myTreePyName == "":
      controller.clearModel()
      return
    namespace = {"root": root}
    cmd = "del(root%s)" % myTreePyName
    res = UXYZ.exec2or3(cmd, namespace)  # TODO python 2-3
    # warning it is del(self)! and so get controller before, (TODO could be event...)
    if controller != None:
      controller.refreshModelViews()
    return

  def delItemSlot(self, status):
    # DBG.write("%s.delItem Slot with no parameters" % self._className, status, True)
    self.delItem()
    return

  def receiveRequest(self, strXmlRequest):
    """
    virtual method: all caseXyz could be a root as a model
    asynchronous treatment of a request from a Controller requestToModelSignal QtCore.pyqtSignal

    in controller
    - create signal and connect it
      >> aController.requestToModelSignal.connect(aModel.receiveRequest)

    in controller
    - have sendRequestToModel method where
      >> aController.requestToModelSignal.emit(aRequest)

    in view (for example)
    - have getController, AND DO NOT directly modify the model
      >> aController = aView.getController()
      >> aRequest = controller.getRequest() # for example
      >> # now define what request details exactly (modify model or else)
      >> aController.sendRequestToModel(aRequest)
    """
    if verbose:
      print("_BaseCommonXyz %s receiveRequest:\n%s" % (self.objectName(), aRequest))
    logger.error("_BaseCommonXyz %s have NOT to to be instancied" % self.objectName())
    return False

  def jsonDumps(self):
    return jsonDumps(self)


########################################################################################
class BaseFreeXyz(_BaseCommonXyz):
  """
  goal is a base class like-standard-class-python plus tree path, plus xml conversion.
  all dynamic new attributes without "_" prefix
  must have parentAsAttribute() method to reach root
  (through __setattr__, like parentAsAttribute() method)
  and so we get possibility to walk on tree:
  >>> root = BaseFreeXyz()
  >>> root.att1 = BaseFreeXyz()
  >>> root.att1.att2 = BaseFreeXyz()
  >>> root.att1.att2.att3 = BaseFreeXyz()
  >>> id( root.att1.att2.att3.parentAsAttribute() ) == id( root.att1.att2 )
  etc...
  this class must be used with immutables types class defined in intFloatListXyx.py, or more inherited.

  there is a choice important, do not change it:
  a == b only if ids are equals
  use instead a.equal(b) == True if contents are equals
  """

  def __init__(self, *args): #nothing to do with args for yet, compatibility with immutables init
    self._verbose = verbose #for debug
    self._nameObject = self.__class__.__name__
    self._isSet= True #every times set because free attributes names
    self._parentAsAttribute = None
    self._attributesAsDict = {}

    self._numeroAttribute = 0 #to reach an order of attributes by sequence of create attributes
    self._numeroAttributesAsDict = {}

    self._nameAsAttribute = None #attribute name from parentAsAttribute
    self._className = str(self.__class__.__name__)
    #if instance is root of tree, name "as it was attribute chid name" (in treeView for example)
    self._defautNameAsRoot = self.__class__.__name__
    self._listKnownStylesXml = "withIndex,withTypeClass" #etc...
    self._styleXml = "withIndex,withTypeClass"
    self._controller = None
    self._haveHidden = False
    """
    example:
    self._attributesAsDict = {("model": StrModelXyz), ("restart": Int01Xyz), etc... }
    """

  def __setattr__(self, name, value):
    if name[0] == "_": return object.__setattr__(self, name, value) #a classical attribute
    if self._verbose: logger.info("%s.__setattr__ %s %s" % (self._className, name, value))
    obj = value #default
    if self._verbose: logger.info("%s.__setattr__ try %s %s %s" % (self._className,  name, id(obj), id(value)))

    if hasattr(obj, "_parentAsAttribute"):
      if obj._parentAsAttribute != None: #avoid loop in tree
        raise Exception("%s.__setattr__: parent of value is set yet" % self._className)
    else: #value without _parentAsAttribute: integer for non mutable (for example)
      raise Exception("%s.__setattr__: forbidden type without _parentAsAttribute: %s" % \
                       (self._className, obj.__repr__()))

    if name in list(self.__dict__.keys()): #attribute existing yet
      self.__delattr__(name)

    object.__setattr__(self, name, obj)
    if self._verbose:
      logger.info("BaseFreeXyz.__setattr__ done %s %s\n%s" % \
                  (name, id(obj), PP.pformat(list(self.__dict__.keys()))))
    obj._parentAsAttribute = self #could be risky, wait and see
    done = [j for i, j in list(self._numeroAttributesAsDict.items())]
    if not name in done: # do not works ... in self._attributesAsDict:
      self._numeroAttributesAsDict[self._numeroAttribute] = name
      self._numeroAttribute += 1
    self._attributesAsDict[name] = obj.__class__
    obj._nameAsAttribute = name
    return

  def __delattr__(self,  name):
    if verbose: logger.info("%s.__delattr__ on %s " % (self._className, name))
    try:
      del self._attributesAsDict[name]
    except:
      if len(_messSome) <= 5:
        msg = "%s.__delattr__ (as BaseFreeXyz) with no '%s' in _attributesAsDict.keys()\n" % (self._className, name)
        for i in traceback.format_stack(): msg += i + "\n"
        msg += "__dict__:\n%s\n" % PP.pformat(list(self.__dict__.keys()))
        msg += "_attributesAsDict:\n%s\n" % PP.pformat(list(self._attributesAsDict.keys()))
        logger.error(msg)
        _messSome.append(True)
      pass
    try:
      tmp = super(BaseFreeXyz, self).__getattribute__(name)
    except:
      tmp = None
    if tmp != None:
      tmp._parentAsAttribute = None
      #print "reset _parentAsAttribute on __delattr__",tmp, id(tmp), tmp._parentAsAttribute
    return super(BaseFreeXyz, self).__delattr__(name)

  def __iter__(self):
    """iter as alphabeticat order"""
    for name in sorted(self._attributesAsDict.keys()):
      if name in list(self.__dict__.keys()):
        yield super(BaseFreeXyz, self).__getattribute__(name)
      else:
        logger.warning("%s.__iter__ (as BaseFreeXyz) attribute '%s' is undefined" % (self._className, name))
        continue

  def items(self):
    for name in sorted(self._attributesAsDict.keys()):
      if self._verbose: logger.info("items %s:\n%s" % (name, PP.pformat(list(self.__dict__.keys()))))
      if name in list(self.__dict__.keys()):
        yield name, super(BaseFreeXyz, self).__getattribute__(name)
      else:
        logger.warning("%s.items (as BaseFreeXyz) attribute '%s' is undefined" % (self._className, name))
        continue

  def itemsByNumero(self):
    if self._verbose:
      logger.info("%s.itemsByNumero (as BaseFreeXyz)\n%s" % \
                  (self._className, PP.pformat(sorted(self._numeroAttributesAsDict.keys()))))
    for no in sorted(self._numeroAttributesAsDict.keys()):
      name = self._numeroAttributesAsDict[no]
      # if self._verbose: print("BaseFreeXyz.itemsByNumero", name, list(self.__dict__.keys()))
      if name in list(self.__dict__.keys()):
        yield name, super(BaseFreeXyz, self).__getattribute__(name)
      else:
        raise Exception("%s.itemsByNumero attribute '%s' is undefined" % (self._className, name))

  def _getNextValidLine(self, stream, verbose=False):
    """
    to reach next line in stream
    each empty line or line starting with # or ! or & is ignored
    """
    line = ""
    while True:
      line = stream.readline()
      if not line:
        if verbose: logger.info("%s._getNextValidLine: eof reached" % self._className)
        raise Exception("End of File: " + stream.name)
      line = line.strip()
      if verbose: logger.info("_getNextValidLine: %s" % line)
      if "#" == line[0]: continue
      if "!" == line[0]: continue #namelist
      if "&" == line[0]: continue #namelist
      if "" == line: continue
      return line

  def __repr__(self):
    #__repr__ goal is to be unambiguous
    aList = [] #an ordered list representation is better for test (and visualize) (in)equality
    for k in sorted(self.__dict__.keys()):
      if k[0] != '_' or k == '_isSet':
        aList.append( [k, self.__dict__[k].__class__.__name__, self.__dict__[k]]  )
    return self._className + " = " + aList.__repr__()

  def equal(self, b):
    """
    as identity on __repr__
    """
    res = True
    if b == None: return False
    if not issubclass(b.__class__, BaseFreeXyz):
      raise Exception("%s instance is not agreed for test equality with %s" % (b.__class__.__name__, self._className))
    ii_cl = self.__class__.__name__
    jj_cl = b.__class__.__name__
    if ii_cl != jj_cl: return False
    if self.__repr__() != b.__repr__(): return False
    return True

  """
  do not do that! risky use a.equal(b)
  def __ne__(self, b):

  def __ne__(self, b):
    # https://docs.python.org/2/reference/datamodel.html
    # Accordingly, when defining __eq__(), one should also define __ne__()
    # so that the operators will behave as expected
    return not self.__eq__(b)
  """

  def check_equality(self, b):
    """idem equal but return first not equal as CheckBaseFreeXyz(ok, why)"""
    """
    as identity on __repr__
    """
    if b == None: return False
    if not issubclass(b.__class__, BaseFreeXyz):
      raise Exception("%s instance is not agreed for test equality with %s" % (b.__class__.__name__, self._className))
    ii_cl = self.__class__.__name__
    jj_cl = b.__class__.__name__
    if ii_cl != jj_cl:
      return CheckBaseXyz(False, "Are not same class '%s' != '%s'" % (ii_cl, jj_cl) )
    ii_repr = self.__repr__()
    jj_repr = b.__repr__()
    if ii_repr != jj_repr:
      #return CheckBaseXyz(False, "Are different")
      #return CheckBaseXyz(False, "Are different:\n%s\n%s" % (ii_repr, jj_repr) )
      why = ""
      strXml1 = self.toStrXml().split('\n')
      strXml2 = b.toStrXml().split('\n')
      for i in range(len(strXml1)):
        l1 = strXml1[i]
        l2 = strXml2[i]
        if l1 != l2:
          return CheckBaseXyz(False, "Are Different in toStrXml line: %i\n%s\n%s" % (i, l1, l2) )
      return CheckBaseXyz(False, "Are different in repr:\n%s\n%s" % (ii_repr, jj_repr) )
    return CheckBaseXyz(True, "Are equal")

  def getCurrentAttributes(self):
    """
    sorted current attributes for a may be incompletly setted inheritance class
    alphabetical ordered attributes output: this method could be overriden
    return ["firstAttrName", "secondAttrName", ...]
    """
    if self.isSet():
      attributes = [k for k in sorted(self._attributesAsDict.keys())]
      return sorted(attributes)
    else:
      raise Exception("%s : is not set" % self._className)

  def getAttributes(self):
    """
    sorted all attributes for a expected completly setted inheritance class
    alphabetical ordered attributes output: this method could be overriden
    return ["firstAttrName", "secondAttrName", ...]
    """
    return self.getCurrentAttributes()
    #this is for next BaseXyz
    if self.isSet():
      attributes = [k for k, _  in self._attributesAsList]
      return sorted(attributes)
    else:
      raise Exception("%s : is not set" % self._className)

  '''def yygetActionsContextMenu(self):
    if verbose: print("%s.getContextMenu" % str(self.__class__.__name__))
    actions = []
    actions.append( self._createAction('Edit dialog', None, 'Edit dialog widget', self.editDialog, None) )
    actions.append( self._createAction('Delete child items', None, 'Delete all child items', self.delAllItems, None) )
    actions.append( self._createAction('Delete', None, 'Delete me', self.delItem, None) )
    return actions'''

  def delAllItems(self):
    selfClassName = self.__class__.__name__
    for nameAttr in list(self._attributesAsDict.keys()):
      delattr(self, nameAttr)
    controller = self.getController()
    if controller != None:
      controller.refreshModelViews()
    return


  def readElementTreeContents(self, tag, listsByTag=[], leaf=None, verbose=False):
    """
    generic reading elementtree tag recursively
    to get model make with BaseFreeXyz(s) and ListofBaseXyz(s)
    usage:
    res = BaseFreeXyz()
    res.readContents(tag, listsByTag)
    listsByTag=['time', ...] (for example) is expected tag names
    with multiples xml brother occurrences,
    stored as named 'time_s'  ListofBaseXyz (for example)
    """
    if verbose: print("readElementTreeContents %s" % tag.tag)
    if leaf is None:
      typeleaf = CLFX.getXyzClassFromName("StrNoEditionXyz")
    else:
      typeleaf = leaf
    attrs = False
    for k, v in tag.items():  # attributes
      if verbose: print("  %s attribute %s=%s" % (tag.tag, k, v))
      self.__setattr__(k, typeleaf(v))
      attrs = True

    children = {}
    for t in tag.getchildren():
      child = BaseFreeXyz()
      child.readElementTreeContents(t, listsByTag)
      if t.tag not in listsByTag:
        self.__setattr__(t.tag, child)
      else:
        try:
          alist = children[t.tag]
        except:
          children[t.tag] = []
          alist = children[t.tag]
        # print("+++ %s append" % t.tag)
        alist.append(child)

    for k, tv in children.items():
      tmp = ListOfBaseXyz()
      for vv in tv:
        tmp.append(vv)
      self.__setattr__("%s_s" % k, tmp)

    try:
      text = tag.text.replace(" ", "").replace("\n", "").replace("\t", "")
    except:
      text = ""
    if len(text) == 0: return  # tag.text supposed useless as empty
    # tag.text supposed useful
    text = tag.text.replace("\n", " ").replace("\t", " ")
    # print("text not empty '%s' %s" % (text, len(text)))
    self.__setattr__('text', typeleaf(text))
    return



########################################################################################
class BaseFreeAllXyz(BaseFreeXyz):

  def __setattr__(self, name, value):
    import xyzpy.intFloatListXyz as IFLX

    if name[0] == "_": return object.__setattr__(self, name, value) #a classical attribute
    if self._verbose: logger.info("%s.__setattr__ %s %s" % (self._className, name, value))
    obj = value #default
    if self._verbose: logger.info("%s.__setattr__ try %s %s %s" % (self._className, name, id(obj), id(value)))

    if hasattr(obj, "_parentAsAttribute"):
      if obj._parentAsAttribute != None: #avoid loop in tree
        #better duplicate implicit, useful, but permissive
        if issubclass(obj.__class__, IFLX._XyzImmBase):
          obj = obj.duplicate()
          if obj._parentAsAttribute != None:
            raise Exception("%.__setattr__ parent of immutable value is set yet" % self._className)
        else: #do not know why
          raise Exception("%.__setattr__ parent of complex value is set yet" % self._className)
    else: #value without _parentAsAttribute: integer for non mutable (for example)
      obj = self._toXyz(obj) #convert to one acceptable


    if name in list(self.__dict__.keys()): #attribute existing yet
      self.__delattr__(name)

    object.__setattr__(self, name, obj)
    if self._verbose:
      logger.info("%s.__setattr__ done %s %s\n%s" % \
                  (self._className, name, id(obj), PP.pformat(sorted(list(self.__dict__.keys())))))
    try:
      obj._parentAsAttribute = self #could be risky, wait and see
    except: #no work for list or dict
      pass
    if name not in self._attributesAsDict:
      self._numeroAttributesAsDict[self._numeroAttribute] = name
      self._numeroAttribute += 1
    self._attributesAsDict[name] = obj.__class__
    try:
      obj._nameAsAttribute = name
    except: #no work for list or dict
      pass
    return

  def _toXyz(self, obj):
    import xyzpy.intFloatListXyz as IFLX

    if type(obj) == int: newObj = IFLX.IntXyz(obj)
    elif type(obj) == float: newObj = IFLX.FloatXyz(obj)
    elif type(obj) == str: newObj = IFLX.StrXyz(obj)
    elif type(obj) == bool: newObj = IFLX.BoolXyz(obj)
    elif obj == None: newObj = obj
    elif type(obj) == list:
      newObj = ListOfBaseXyz()
      for i in obj: newObj.append(self._toXyz(i))
    elif type(obj) == dict:
      #newObj = obj #obj mutable, dangerous, illogic, have to be replaced by BaseFreeAllXyz
      raise Exception("%s.__setattr__ forbidden type dict without _parentAsAttribute, fix as have to be replaced by BaseFreeAllXyz\n%s" % \
                      (self._className, obj.__repr__()))
    else:
      raise Exception("%s.__setattr__ forbidden type without _parentAsAttribute\n%s" % \
                      (self._className, obj.__repr__()))
    return newObj

########################################################################################
class ListOfBaseXyz(list, _BaseCommonXyz):
  """
  goal is a base class like-standard-list-python plus tree path, plus xml conversion.

  | All dynamic new element appended must have same type,
  | and must have parentAsAttribute() method to reach root
  | (through __setattr__, like parentAsAttribute() method)
  | and so we get possibility to walk on tree.
  |
  | examples:
  | >>> root = ListOfBaseXyz()
  | >>> a = BaseFreeXyz()
  | >>> root.append(a)
  | >>> id( root[0].parentAsAttribute() ) == id( root )
  | >>> id( root[0] ) == id( a )
  |
  | this class must be used with immutables types class defined in intFloatListXyx.py,
  | or more inherited.
  |
  | there is a choice important, do not change it:
  | a == b only if ids are equals
  | use instead a.equal(b) == True for test contents are equals
  """

  # list of allowed (sub)classes to be append/insert in this list
  # include previously defined only as Class
  # include now can be defined as Class or String
  # example : _allowedClasses = [ IntXyz , "FloatXyz" , ... ]
  _allowedClasses = [] # by default only class of first element appended allowed

  def __init__(self, nameObject=""):
    self._verbose = verbose #for debug
    self._nameObject = self.__class__.__name__ #"ListOfUnknownName" #self.__class__.__name__
    self._isSet= True #every times set because
    self._parentAsAttribute = None
    self._nameAsAttribute = None #attribute name from parentAsAttribute
    self._className = str(self.__class__.__name__)
    #if instance is root of tree, name "as it was attribute chid name" (in treeView for example)
    self._defautNameAsRoot = self.__class__.__name__
    self._listKnownStylesXml = "withIndex,withTypeClass" #etc...
    self._styleXml = "withIndex,withTypeClass"
    self._controller = None
    #if "[]" in nameObject: print plante

  def __setattr__(self, name, value):
    # existing by default in this inherided classes of list
    if name[0] == "_":
      return list.__setattr__(self, name, value) # a classical attribute
    raise Exception("%s.__setattr__ forbidden in this list inherided class" % self._className)

  def _testInAllowedClasses(self, obj):
    if len(self._allowedClasses) == 0:
      # if _allowedClasses not defined, first item is unconditionaly allowed, but others have to be same type.
      if len(self) != 0:
        if obj.__class__.__name__ != self[0].__class__.__name__:
          raise Exception("%s._testInAllowedClasses type %s unexpected, needs %s" % \
                          (self._className, obj.__class__.__name__, self[0].__class__.__name__))
      return True #ok
    for iClass in self._allowedClasses:
      aClass = CLFX.toClass(iClass) # iClass can be string or class 211021
      if issubclass(obj.__class__, aClass):
        return True # ok
    raise Exception("%s._testInAllowedClasses type %s unexpected, needs subclasses in:\n%s" % \
                    (self._className, obj.__class__.__name__, PP.pformat(self._allowedClasses)))

  def appendDefault(self):
    #aClassName = self._allowedClasses[0]
    aClass = self._allowedClasses[0] #CLFX.getXyzClassFromName(aClassName)
    value = aClass()
    self.append(value)

  def append(self, value):
    """only one type in list"""
    obj = value #default
    if type(value) == str: #automatic cast if only one allowed class
      if len(self._allowedClasses) == 1:
        obj = self._allowedClasses[0](value)
    ok = self._testInAllowedClasses(obj) #ko: raise

    if hasattr(obj, "_parentAsAttribute"):
      if obj._parentAsAttribute != None: #avoid loop in tree
        raise Exception("%s.append: parent of value is set yet" % self._className)
    else: #value without _parentAsAttribute: integer for non mutable (for example)
      raise Exception("%s.append: forbidden type without _parentAsAttribute: %s" % (self._className, str(obj)))

    res = super(ListOfBaseXyz, self).append(obj)
    if len(self._allowedClasses) == 0:
      #set instance empty self._allowedClasses as class of first element appended, definitif
      self._allowedClasses = [obj.__class__]
    if self._nameObject == "ListOfBaseXyz": #"ListOfUnknownName":
      self._nameObject = str("ListOf" + obj.__class__.__name__)

    if self._verbose: logger.info("%s.append done %s %s" % (self._className, len(self), type(obj)))
    obj._parentAsAttribute = self #could be risky, wait and see
    #obj._nameAsAttribute = "_IndexListXyz" #means that obj is accessible through list index
    return res

  def insert(self, i, value):
    #print "ListOfBaseXyz.insert", i, value
    obj = value #default
    ok = self._testInAllowedClasses(obj) #ko: raise

    if hasattr(obj, "_parentAsAttribute"):
      if obj._parentAsAttribute != None: #avoid loop in tree
        raise Exception("%s.insert: parent of value is set yet" % self._className)
    else: #value without _parentAsAttribute: integer for non mutable (for example)
      raise Exception("%s.insert: forbidden type %s without _parentAsAttribute method:\n%s" % \
                      (self._className, type(obj), PP.pformat(obj)))

    res = super(ListOfBaseXyz, self).insert(i, value)

    if self._nameObject == "ListOfBaseXyz": #"ListOfUnknownName":
      self._nameObject = str("ListOf" + obj.__class__.__name__)

    if self._verbose: logger.info("%s.insert done at %s %s %s" % (self._className, i, len(self), type(obj)))
    obj._parentAsAttribute = self #could be risky, wait and see
    #obj._nameAsAttribute = "_IndexListXyz" #means that obj is accessible through list index
    return res

  def __setitem__(self, i, value):
    #print "ListOfBaseXyz.__setitem__", i, value
    obj = value #default
    ok = self._testInAllowedClasses(obj) #ko: raise

    if hasattr(obj, "_parentAsAttribute"):
      if obj._parentAsAttribute != None: #avoid loop in tree
        raise Exception("%s.__setitem__: parent of value is set yet" % self._className)
    else: #value without _parentAsAttribute: integer for non mutable (for example)
      raise Exception("%s.__setitem__: forbidden type without _parentAsAttribute: %s" % (self._className, str(obj)))

    res = super(ListOfBaseXyz, self).__setitem__(i, value)
    if self._nameObject == "ListOfBaseXyz": #"ListOfUnknownName":
      self._nameObject = str("ListOf" + obj.__class__.__name__)

    if self._verbose: logger.info("%s.__setitem__ done %s %s" % (self._className, len(self), type(obj)))
    obj._parentAsAttribute = self #could be risky, wait and see
    #obj._nameAsAttribute = "_IndexListXyz" #means that obj is accessible through list index
    return res

  def itemsByNumero(self):
    if self._verbose: logger.info("%s.itemsByNumero len %s" % (self._className, len(self)))
    for no in range(len(self)):
      name = "_%i_" % no
      # if self._verbose: print("BaseFreeXyz.itemsByNumero", name)
      yield name, self[no]

  def items(self):
    for no in range(len(self)):
      name = "_%i_" % no
      if self._verbose: logger.info("%s.items %s" % (self._className, name))
      yield name, self[no]

  def __repr__(self):
    aList = [] #an ordered list representation is better for test (and visualize) (in)equality
    for i in self:
      aList.append( i )
    return self._nameObject + " = " + aList.__repr__()

  def equal(self, b):
    """
    as identity on __repr__
    """
    res = True
    if b == None: return False
    if not issubclass(b.__class__, ListOfBaseXyz):
      raise Exception("%s.equal: type %s is not agreed for test equality"  % \
                      (self._className, b.__class__.__name__))
    if len(self) != len(b): return False
    if self.__repr__() != b.__repr__(): return false
    return True

  """
  do not do that! risky use a.equal(b)
  def __ne__(self, b):

  def __ne__(self, b):
    # https://docs.python.org/2/reference/datamodel.html
    # Accordingly, when defining __eq__(), one should also define __ne__()
    # so that the operators will behave as expected
    return not self.__eq__(b)
  """

  def check_equality(self, b):
    """idem equal but return first not equal as CheckBaseFreeXyz(ok, why)"""
    """
    as identity on __repr__
    """
    res = True
    if b == None: return False
    if not issubclass(b.__class__, ListOfBaseXyz):
      raise Exception("%s.check_equality: type %s is not agreed for check equality"  % \
                      (self._className, b.__class__.__name__))
    if len(self) != len(b):
      return CheckBaseXyz(False, "have not same length")
    for i in range(len(self)):
      ii = self[i]
      jj = b[i]
      ii_repr = ii.__repr__()
      jj_repr = jj.__repr__()
      if ii_repr != jj_repr:
        return CheckBaseXyz(False, "have element %i different:\n%s\n%s" % (i, ii_repr, jj_repr))
    return CheckBaseXyz(True, "Are equal")

  def getActionsContextMenu(self):
    actions = _BaseCommonXyz.getActionsContextMenu(self)
    ii = 0
    for iClass in self._allowedClasses:
      aClass = CLFX.toClass(iClass)() # instancier pour avoir le nom str(iClass) durdur!
      name = aClass._defautNameAsRoot
      aText = 'Append item %s' % name
      aIcon = aClass._icon
      logger.debug("create lambda action for %s in %s" % (iClass, self.__class__.__name__))
      # new in PyQt5 state override in lambda
      # https://stackoverflow.com/questions/35819538/using-lambda-expression-to-connect-slots-in-pyqt
      #tmp = self._createAction(aText, None, aText, lambda status=None, lClass=aClass.__class__: self.addItemSlot(status, lClass), aIcon)
      tmp = self._createAction(aText, None, aText, partial(self.addItemSlot, None, aClass.__class__), aIcon)
      actions.append(tmp)
      ii += 1
    return actions

  def addItemSlot(self, status, aClass=None):
    """new in PyQt5 state override in lambda"""
    return self.addItem(aClass)

  def addItem(self, aClass=None):
    logger.debug("%s.addItem %s" % (self._className, aClass))
    if aClass == None:
      logger.error("%s.addItem '%s' not in allowed classes:\n%s" % \
                   (self._className, "None", PP.pformat(self._allowedClasses)))
      return
    try:
      item = aClass()
    except:
      raise Exception("%s.addItem problem adding '%s', allowed classes:\n%s" % \
                      (self._className, type(aClass), PP.pformat(self._allowedClasses)))
    itemClassName = item.__class__.__name__
    logger.debug("%s.addItem standard by default for '%s'" % (self._className, itemClassName))
    ok = False
    for c in self._allowedClasses:
      cClass = CLFX.toClass(c)
      if aClass == cClass:
          ok = True
          break
    if not ok:
      logger.error("%s.addItem '%s' not in allowed classes:\n%s" % \
                   (self._className, itemClassName, PP.pformat(self._allowedClasses)))
      return
    item.setDefaultValues() # if problem exception go to caller
    self.append(item)
    controller = self.getController()
    if controller != None:
      if verbose: logger.info("%s there is controller and refresh views" % self._className)
      controller.refreshModelViews()
      # listof itself and new element and childs
      new = self.getTreePyName() + r"[%i]" % (len(self) - 1)
      expanded = [self.getTreePyName(),
                  new,
                  new + ".",
                 ]
      logger.debug("expand %s" % expanded)
      controller.ExpandSignal.emit(expanded)
    return

  def insertItem(self, iindex, aClass=None):
    # DBG.write("%s.insertItem" % self._className, [iindex, aClass])
    if aClass==None:
      logger.error("%s.insertItem '%s' not in allowed classes:\n%s" % \
                   (self._className, "None", PP.pformat(self._allowedClasses)))
      return
    sClass = CLFX.toString(aClass)
    cClass = CLFX.toClass(aClass)
    itemClassName = sClass
    logger.debug("%s.insertItem standard by default for '%s'" % (self._className, itemClassName))
    if cClass not in CLFX.toClassList(self._allowedClasses):
      logger.error( "%s.insertItem '%s' not in allowed classes:\n%s" % \
                   (self._className, itemClassName, PP.pformat(self._allowedClasses)))
      return

    item = cClass()
    item.setDefaultValues()
    self.insert(iindex, item)
    controller = self.getController()
    if controller != None:
      if verbose: logger.info("%s.insertItem there is controller and refresh views" % self._className)
      controller.refreshModelViews()
      # listof itself and new element and childs
      new = self.getTreePyName() + r"[%i]" % (iindex)
      expanded = [self.getTreePyName(),
                  new,
                  new + ".",
                 ]
      logger.debug("expand %s" % expanded)
      controller.ExpandSignal.emit(expanded)
    return

  def insertItemSlot(self, status, iindex, aClass=None):
    """
    for pyqt4 to pyqt5 change status as boolean
    https://stackoverflow.com/questions/28319889/pyqt5-whats-this-boolean-passed-when-a-menu-action-is-triggered
    """

    # DBG.write("%s.insertItemSlot" % self._className, [status, iindex, aClass])
    self.insertItem(iindex, aClass)
    return

  def delAllItems(self):
    for i in reversed(list(range(len(self)))):
      del(self[i])
    controller = self.getController()
    if controller != None:
      controller.refreshModelViews()
    return

  def delItem(self, index=None):
    # DBG.write("%s.delItem index" % self._className, index, True)
    if index == None:
      raise Exception("%s.delItem needs index not None" % self._className)
    del(self[index])
    controller = self.getController()
    if controller != None:
      controller.refreshModelViews()
    return

  def delItemSlot(self, status, index=None):
    # DBG.write("%s.delItemSlot index" % self._className, index, True)
    self.delItem(index)
    return

  def setDefaultValues(self):
    """not virtual method, do nothing. Default value is initial empty list []"""
    pass


########################################################################################
class _XyzConstrainBase(BaseFreeXyz):
  """
  goal is a not-to-be-instancied base class like-standard-class-python plus tree path, plus xml conversion.
  all constrained typed (NOT dynamic) attributes without "_" prefix
  must have parentAsAttribute() method to reach root
  (through __setattr__, like parentAsAttribute() method)
  and so we get possibility to walk on tree:
    root = BaseFreeXyz()
    root.att1 = BaseFreeXyz()
    root.att1.att2 = BaseFreeXyz()
    root.att1.att2.att3 = BaseFreeXyz()
    id( root.att1.att2.att3.parentAsAttribute() ) == id( root.att1.att2 )
    etc...

  this class must be use only for inheritance of constrained user new classes.

  this class must be used with immutables types class defined in intFloatListXyx.py, or more inherited.

  there is a choice important, do not change it:
  a == b only if ids are equals
  use instead a.equal(b) == True if contents are equals
  """

  #list of allowed attributes, not a dict because often sequential order list is used in asci data files to read
  _attributesList = [] #empty in this not-to-be-instancied base class

  #dict of help for allowed attributes.
  _helpDict = {} #empty in this not-to-be-instancied base class
  """example:
  _attributesList = [
    ("model", "StrModelXyz"),
    ("restart", Int01Xyz), etc...
  ]
  _helpDict = {
    "model": ("a short blabla for Help", "a long blabla for Help"),
    "restart": ("a short blabla for Help", "a long blabla for Help"), etc...
  }
  """

  def __init__(self):
    super(_XyzConstrainBase, self).__init__()
    #never set because this class must be use only for inheritance of constrained user new classes
    self._isSet = False #never set because this class must be use only for inheritance of constrained user new classes
    self._isCast = False #if True try to cast attribute to expected type, else need strictly good type for it
    self._attributesDict=dict((key, value) for (key, value) in self._attributesList)

  def _getAllowedAttributeNames(self):
    return [i for i, j in self._attributesList]

  def _getStrAllowedAttributeNames(self):
    return PP.pformat(self._getAllowedAttributeNames())

  def __setattr__(self, name, value):
    PPFAN = self._getStrAllowedAttributeNames # shortcut to method
    if name[0] == "_": return object.__setattr__(self, name, value) #a classical attribute
    if self._verbose: logger.info("%s.__setattr__ %s %s" % (self._className, name, value))
    obj = value #default
    if self._verbose: logger.info("%s.__setattr__ try %s %s %s" % (self._className, name, id(obj), id(value)))

    if name not in list(self._attributesDict.keys()):
      raise Exception("%s.__setattr__  unexpected attribute name '%s' not in allowed attributes:\n%s" % \
                      (self._className, name, PPFAN()))

    valueClass = value.__class__.__name__
    expectedClass = self._attributesDict[name]

    if not self._isCast:
      if valueClass != expectedClass:
        raise Exception("%s.__setattr__ unexpected type class %s for attribute %s (expected %s)" % \
                        (self._className, valueClass, name, expectedClass))
    else:
      if valueClass != expectedClass:
        if True: #try: TODO for debug
          if verbose:
            logger.debug("%s.__setattr__ create and set attribute with new casted instance %s for value, not value %s herself." % \
                           (self._className, expectedClass, valueClass))
          obj = CLFX.getXyzInstanceClassFromNameAndValue(expectedClass, value)
        else: #except: TODO for debug
          raise Exception("%s.__setattr__ unexpected and uncastable class %s for attribute %s (expected %s)" % \
                          (self._className, valueClass, name, expectedClass))

    if hasattr(obj, "_parentAsAttribute"):
      if obj._parentAsAttribute != None: #avoid loop in tree
        raise Exception("%s.__setattr__  parent of value is set yet" % self._className)
    else: #value without _parentAsAttribute: integer for non mutable (for example)
      raise Exception("%s.__setattr__ forbidden type without _parentAsAttribute: " % \
                      (self._className, obj.__repr__()))

    if name in list(self.__dict__.keys()): #attribute existing yet
      if verbose: logger.info("%s.__setattr__ attribute '%s' existing yet, delete before setattr" % \
                              (self._className, name))
      self.__delattr__(name)

    object.__setattr__(self, name, obj)
    if self._verbose:
      logger.info("%s.__setattr__ done %s %s:\n%s" % \
                  (self._className, name, id(obj), PP.pformat(sorted(list(self.__dict__.keys())))))
    obj._parentAsAttribute = self #could be risky, wait and see
    self._attributesAsDict[name] = obj.__class__
    obj._nameAsAttribute = name
    return

  def setIsCast(self, aBool):
    """used in inherided class"""
    if not issubclass(aBool.__class__, bool):
      raise Exception("%s.setIsSet argument is not a boolean" % self._className)
    self._isCast = aBool

  def __iter__(self):
    """user defined order of existing attributes, not sorted"""
    for name, _ in self._attributesList:
      if name not in list(self._attributesAsDict.keys()): continue
      try:
        res = super(_XyzConstrainBase, self).__getattribute__(name)
      except:
        logger.warning("%s.__iter__ (as XyzConstrainBase) attribute '%s' is undefined" % (self._className, name))
        continue
      yield res

  def itemsByNumero(self):
    """order is same as iteritem"""
    if self._verbose: logger.info("%s.itemsByNumero\n%s" % (self ._className, PP.pformat(self._attributesList)))
    # if self._verbose: print("XyzConstrainBase.itemsByNumero as _attributesList order", self._className, self._attributesList)
    for name, _ in self._attributesList:
      if name not in list(self._attributesAsDict.keys()): continue
      try:
        res = (name, super(_XyzConstrainBase, self).__getattribute__(name))
      except:
        # if self._verbose: print("Warning: " + self._className + ": attribute '" + name + "' is undefined")
        continue
      yield res

  def _setFromList(self, attributeName, line, classNameAppend):
    """
    | example:
    | attributeName is ListOfBaseXyz
    | line is 'value1 value2 ...'
    | >>> self._setFromList(attributeName, line, "FloatPosXyz")
    """
    tmp = line.split()
    aClass = self._getClassFromAttributeName(attributeName)
    anAttribute = aClass()
    anAppendClass = CLFX.getXyzClassFromName(classNameAppend)
    for i in tmp:
      anAppendInstance = anAppendClass(i)
      anAttribute.append(anAppendInstance)
    setattr(self, attributeName, anAttribute)

  def _getClassFromAttributeName(self, attributeName):
    try:
      aClassName = self._attributesDict[attributeName]
    except:
      raise Exception("%s._getClassFromAttributeName attribute '%s' unknown" % (self._className, attributeName))
    aClass = CLFX.getXyzClassFromName(aClassName)
    if aClass == None:
      raise Exception("%s._getClassFromAttributeName attribute '%s' unknown class '%s' in factory" % \
                      (self._className, attributeName, aClassName))
    else:
      return aClass

  def _getAttributesInLine(self, attributes, line):
    tmpatt = attributes.split()
    tmpval = line.split()
    try:
      for i in range(len(tmpatt)):
        attributeName = tmpatt[i]
        value = tmpval[i]
        setattr(self, attributeName, value)
      return
    except:
      raise Exception("%s._getAttributesInLine problem for attributes '%s' with line\n'%s'" % \
                      (self._className, attributes, line))

  def _setAllAttributesList(self):
    for attributeName, _ in self._attributesList:
      aClass = self._getClassFromAttributeName(attributeName)
      setattr(self, attributeName, aClass() )

  def getActionsContextMenu(self):
    actions = super(_XyzConstrainBase, self).getActionsContextMenu()
    actions.append( self._createAction('Reset inexisting child items', None, 'Reset expected inexisting child item', self.resetExpectedItems, 'resetinexistingchilditem') )
    return actions

  def resetExpectedItems(self):
    controller = self.getController()
    for attributeName, _ in self._attributesList:
      #print "attributeName",attributeName
      if not hasattr(self, attributeName):
        new = self._getClassFromAttributeName(attributeName)()
        new.setDefaultValues()
        setattr(self, attributeName, new )
    if controller != None:
      controller.refreshModelViews()
    return

  def setDefaultValues(self):
    """not virtual, could be used"""
    if verbose: logger.info("%s.setDefaultValues (as ConstrainBase)" % self._className)
    for attributeName, _ in self._attributesList:
      #print "attributeName",attributeName
      new = self._getClassFromAttributeName(attributeName)()
      #new.setDefaultValues()
      setattr(self, attributeName, new )

  def getToolTips(self):
    """get tooltips for attributes as ordeded list as _attributesList"""
    res = []
    for name, _ in self._attributesList:
      try:
        tooltip, _ = self._helpDict[name]
      except:
        logger.warning("No tooltip for %s" % self._className)
        tooltip = "No tooltip"
      res.append(tooltip)
    return res

########################################################################################
class CheckBaseXyz(_XyzConstrainBase): #could be listed, and xml-streamed

  _attributesList = [ #list, not a dict because sequential order list is used in files cnf
    ("ok", "BoolXyz"),
    ("why", "StrXyz"),
    ("traceback", "StrXyz"),
  ]

  def __init__(self, ok=False, why="Uknown status"):
    super(CheckBaseXyz, self).__init__()
    self.setIsCast(True)
    self.setIsSet(True)
    if type(ok) is not bool:
      raise Exception("%s init with ok not bool: %s" % (self._className, ok))
      # logger.error("%s init with ok not bool: %s" % (self._className, ok))
      self.ok = False
    else:
      self.ok = ok
    self.why = why
    self._traceback = ""

  def reset(self):
    self.ok = False
    self.why = "Uknown status"
    self.traceback = ""
    self.setIsSet(True)

  def setDefaultOk(self, why="Default ok status"):
    self.ok = True
    self.why = StrXyz(why)
    self._traceback = StrXyz("")
    self.setIsSet(True)

  def getValue(self):
    return (self.ok.getValue(), self.why.getValue())



################################################################"
class HelpXyz(BaseFreeXyz):
  """
  general instance to work with help data from _helpDict to __commonHelps__
  """

  """ no needs... may be...
  _attributesList = [
    ("shortHelp", "StrXyz"),
    ("longHelp", "StrXyz"),
  ]

  _helpDict = {
    "shortHelp": ("a simple text help", "like a small tooltip"),
    "longHelp": ("a long text help", "like a big tooltip"),
   }
  """

  def __init__(self):
    super(HelpXyz, self).__init__()


#factory pattern using xyzpy.utilsXyz._dictOfXyzClass
CLFX.appendAllXyzClasses( [ BaseXyz, BaseFreeXyz, ListOfBaseXyz, HelpXyz ] )
