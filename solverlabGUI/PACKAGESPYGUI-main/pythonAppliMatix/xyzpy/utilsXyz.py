#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
utilities for output format text XML for xyzpy
"""

import os
import sys
import pprint as PP
import subprocess as SP
import traceback
import xml.etree.ElementTree as ET
import platform
from datetime import datetime

import xyzpy.loggingXyz as LOG
import returncodepy.returnCode as RCO
import xyzpy.stringIO as IOX


logger = LOG.getLogger()
aMutableForOneMessage = []

verbose = False
verboseExec = False # verbose for warning of meta programmation call exec()

myDir = os.path.split(os.path.realpath(__file__))[0]
originDir = myDir.split("PACKAGESPY")[0]

#default formats for treeViewXyz
FORMATS_TREEVIEW = {
  "FMT_FLOAT": "{0:13.5e}",  #"{0:13.5g}", for cassis...
  "FMT_INT": "{0:>13}",
}

def currentFuncName(n=0):
  """
  for current caller func name, specify 0 or no argument.
  for name of caller of caller current func, specify 1.
  etc.
  """
  return sys._getframe(n + 1).f_code.co_name

def getBrowser(basename=False):
  res = os.getenv("USER_WEBBROWSER")
  if verbose: print("$USER_WEBBROWSER", res)
  if res is None:
    if platform.system() == "Windows":
      res = r'"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe"'
      if not os.path.exists(res[1:-1]): # skip '"'
        res = r'"C:\Program Files\Internet Explorer\iexplore.exe"'
    else: # linux etc
      res = "firefox" # or "xdg-open"
  if verbose: print("WEB BROWSER", platform.system(), res)
  if basename:
    # from r'"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe"' to 'chrome'
    return os.path.splitext(os.path.basename(res))[0]
  else:
    return res

def getFileBrowser(basename=False):
  res = os.getenv("USER_FILEBROWSER")
  if verbose: print("$USER_FILEBROWSER", res)
  if res is None:
    if platform.system() == "Windows":
      res = r'"C:\WINDOWS\explorer.exe"'
      if not os.path.exists(res[1:-1]): # skip '"'
        res = r'"C:\Program Files\Internet Explorer\iexplore.exe"'
    else: # linux etc
      res = "caja" # or "xdg-open"
  if verbose: print("FILE BROWSER", platform.system(), res)
  if basename:
    # from r'"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe"' to 'chrome'
    return os.path.splitext(os.path.basename(res))[0]
  else:
    return res

def getEditor(basename=False):
  res = os.getenv("USER_EDITOR")
  if verbose: print("$USER_EDITOR", res)
  if res is None:
    if platform.system() =="Windows":
      res = r'"C:\Program Files (x86)\Notepad++\notepad++.exe"'
      if not os.path.exists(res[1:-1]): # skip '"'
        res = r'"C:\Windows\System32\notepad.exe"'
      os.environ["USER_EDITOR"] = res
    else: # linux etc
      res = getInstalledLinuxDefaultEditor()
    os.environ["USER_EDITOR"] = res
    logger.debug("default editor (set from $USER_EDITOR) %s is %s" % (platform.system(), res))
  if basename:
    # from r'"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe"' to 'chrome'
    return os.path.splitext(os.path.basename(res))[0]
  else:
    return res

def getInstalledLinuxDefaultEditor():
  editors = "pluma kate gedit emacs".split()
  for edi in editors:
    rc = Popen("which %s" % edi) # , logger=logger)
    if rc.isOk():
      # logger.debug("getInstalledLinuxDefaultEditor is %s" % rc.getValue())
      return edi
  return "$USER_EDITOR not set"

def toStrForXml(value, newLineEvery=10):
    """
    goal is to write an xml tag.text, as data formatted, almost human eyes readable:
    
    - exponent notation for float ' -1.00000e+00' 13 digits,
      i.e for int  etc...         '            1' 
      (allows aligned columns)
    - exception for type str: always a first blank if len > 13 (without truncature)
      best for python fortran free format read etc...
    - value type accepted are only elementary immutable or list of immutables,
      (and intFloatListXyz types, for example).
      unless raise Exception
    
    warning: python bool(str) is True if str != "" and so use bool with attention:
    
    - toStrForXml(aBoolValue) return "True" or "False" and so
      you have not to use bool(toStrForXml(aBoolValue)) which ever True
    - use intFloatListXyz.Bool(toStrForXml(aBoolValue)).getValue() instead (!sorry!)
    """
    
    #if type(value)==list:
    #if hasattr(value, "_toStr"):
    #  return value._toStr()
    if issubclass(value.__class__, list):
      res=""
      ii=0
      for i in value:
        res += toStrForXml(i)
        ii += 1
        if ii%newLineEvery == 0: # almost human eyes readable: newline every 5 items
          res += "\n"
          ii=0
      return res
    if issubclass(value.__class__, bool):
      if value:
        return "{0:>13}".format("True")
      else:
        return "{0:>13}".format("False")
    if issubclass(value.__class__, str):
      res = "{0:>13}".format(value)
      if res[0] != " ": res = " " + res # always at least a blanks
      return res
    if issubclass(value.__class__, int):
      return "{0:>13}".format(value)
    if issubclass(value.__class__, float):
      return "{0:13.5e}".format(value)
    if value is None:
      return "None"
    # have to implement _toStrForXml(self) in class value
    raise Exception("toStrForXml type not implemented (works only for immutables): " + str(type(value)))

def toStrForColumns(value):
  return toStrForXml(value, newLineEvery=20)

def toStrForTreeView(value, formats=FORMATS_TREEVIEW):
  """
  filter dislpay role ot column value
  value is string, tries as int, but float, but i.e toStrForXml
  """
  try:
    return formats["FMT_INT"].format(int(value))
  except:
    pass
  try:
    return formats["FMT_FLOAT"].format(float(value))
  except:
    pass
  return toStrForXml(value)

def decodeBytes(aStr):
  if type(aStr) == bytes:
    # print("\ndecodeBytes %s..." % aStr[0:20])
    return aStr.decode("utf-8") # as str
  else:
    return aStr

def prettyPrintETMageia4PyXML(aRootOrNode):
  """
  work fine but needs PyXML distribution because xml.dom.ext sometimes not present
  """
  from xml.dom.ext.reader import Sax2
  from xml.dom.ext import PrettyPrint
  reader = Sax2.Reader()
  docNode = reader.fromString(ET.tostring(aRootOrNode))
  tmpStream = IOX.StringIO()
  PrettyPrint(docNode, stream=tmpStream)
  res = tmpStream.getvalue()
  return decodeBytes(res)

def prettyPrintETSalome730(aRootOrNode):
  """
  | works only from Python 2.7 for salome 730
  | see https://mail.python.org/pipermail/tutor/2010-August/078020.html
  """
  from xml.dom import minidom
  if isinstance(aRootOrNode, ET.Element("").__class__):
    txt = ET.tostring(aRootOrNode)
    try:
      res = minidom.parseString(txt).toprettyxml(encoding='utf-8', indent="  ")
    except:
      traceback.print_exc() #better explicit verbose problem
      print("\nproblem minidom.parseString(txt):\n")
      """
      ii=1
      for i in txt.split('\n'): 
        print("%i//%s//" % (ii, i))
        ii += 1
      """
    return decodeBytes(res)
  if isinstance(aRootOrNode, ET.ElementTree().__class__):
    txt = ET.tostring(aRootOrNode.getroot())
    res = minidom.parseString(txt).toprettyxml(encoding='utf-8', indent="  ")
    return decodeBytes(res)
  raise Exception("prettyPrintETSalome730: problem type(aRootOrNode): " + str(type(aRootOrNode)) + "not implemented")

def prettyPrintET(aRootOrNode):
  try:
      from xml.dom.ext import PrettyPrint
      PyXML = True
  except:
      PyXML = False
  if PyXML: 
      return prettyPrintETMageia4PyXML(aRootOrNode) #if package set PyXML distribution
  else:
     if aMutableForOneMessage == []:
       # user = os.getenv("USERNAME")
       # if user == "wambeke": #only for me
       #   print "\nWARNING: prettyPrintET: unknown PyXML import: get prettyPrintET as salome730 centos for", platformDist[0:2]
       aMutableForOneMessage.append("done")
     return prettyPrintETSalome730(aRootOrNode)

def stripAll(aStr):
  """remove all space and let newlines tabs etc..."""
  return aStr.replace(" ","")

def setDefaultStyleXml(item, kwargs):
  """
  | used for toXml details of method toXml(item)
  | if item is simple class without attribute _styleXml,
  | set default key styleXml at "withIndex,withTypeClass"
  """
  if "styleXml" not in kwargs:
    try:
      kwargs["styleXml"] = item._styleXml
    except:
      kwargs["styleXml"] = "withIndex,withTypeClass" #default if simple class dict or list
  return
    
def toXml(item, **kwargs):
  """
  | create an ET.Elementtree from item
  | item is a baseXyz or subclasses or dict or list or int or float or str ...
  | kwarg are for optional future option of added details in xml tree
  | (see setDefaultStyleXml)
  """
  
  import xyzpy.baseXyz as BXYZ
  import xyzpy.intFloatListXyz as IFLX

  setDefaultStyleXml(item, kwargs)
  if issubclass(item.__class__, IFLX._XyzImmBase) or \
     issubclass(item.__class__, BXYZ.ListOfBaseXyz) or \
     issubclass(item.__class__, BXYZ.BaseFreeXyz) or \
     issubclass(item.__class__, BXYZ.BaseXyz):
    res = _toXmlFromXyz(item, **kwargs)
    return res
    
  if issubclass(item.__class__, list): #simple list
    res = _toXmlFromList(item, **kwargs)
    return res

  if issubclass(item.__class__, tuple): #as simple list
    res = _toXmlFromTuple(item, **kwargs)
    return res

  if issubclass(item.__class__, dict): #simple Dist
    res = _toXmlFromDict(item, **kwargs)
    return res

  if issubclass(item.__class__, int): #simple Dist
    res = toXml(IFLX.IntXyz(item), **kwargs)
    return res

  if issubclass(item.__class__, float): #simple Dist
    res = toXml(IFLX.FloatXyz(item), **kwargs)
    return res

  if issubclass(item.__class__, str): #simple Dist
    res = toXml(IFLX.StrXyz(item), **kwargs)
    return res
  
  raise Exception("toXml: problem type '%s' unexpected" % str(type(item)))

def _toXmlFromXyz(item, **kwargs):
  """
  | create an ET.Elementtree from item
  | item is a baseXyz or ListOfBaseXyz or_XyzImmBase and subclasses
  | kwarg are for optional future option of added details in xml tree
  | have to be called from method toXml_
  """
  
  import xyzpy.baseXyz as BXYZ
  import xyzpy.intFloatListXyz as IFLX

  if not (issubclass(item.__class__, IFLX._XyzImmBase) or \
          issubclass(item.__class__, BXYZ.ListOfBaseXyz) or \
          issubclass(item.__class__, BXYZ.BaseFreeXyz) or \
          issubclass(item.__class__, BXYZ.BaseXyz)):
    raise Exception("_toXmlFromXyz: problem type '%s' unexpected" % str(type(item)))
  
  setDefaultStyleXml(item, kwargs)
  nameTag = item.getNameAsAttribute()
  
  if nameTag == None:
    nameTag = item.getNameAsAttributeAsRoot() #getNameObject()
      
  res = ET.Element(nameTag)
  
  parent = item.parentAsAttribute()
  if parent != None:
    hidden = parent.isHidden(nameTag)
  else:
    hidden = False
  if hidden == True: #not hidden by default
    res.attrib["hidden"] = str(hidden)

  if "withoutTypeClass" not in kwargs["styleXml"]:
    res.attrib["typeClass"] = item._className
  #TODO if "withTreePyName" in kwargs["styleXml"]:
  #print "_toXmlFromXyz kwargs", kwargs["styleXml"]
  if "withoutTreePyName" not in kwargs["styleXml"]:
    res.attrib["treePyName"] = item.getTreePyName()
  index = 0
  
  if item._icon != None: res.attrib["icon"] = item._icon
  
  if issubclass(item.__class__, IFLX._XyzImmBase):
    res = item.toXml(**kwargs)
    return res
  
  withIndex = False
  if issubclass(item.__class__, BXYZ.ListOfBaseXyz):
    withIndex = True
    
  for i in item:
    ires = i.toXml(**kwargs)
    if withIndex: ires.attrib["index"] = str(index)
    res.append(ires)
    index += 1
  if withIndex: res.attrib["size"] = str(len(item))
  return res

def _toXmlFromList(item, **kwargs):
  """
  | create an ET.Elementtree from item:
  | item is a list orsubclasses
  | kwarg are for optional future option of added details in xml tree
  | have to be called from method toXml_
  """

  import xyzpy.baseXyz as BXYZ
  import xyzpy.intFloatListXyz as IFLX

  if issubclass(item.__class__, BXYZ.ListOfBaseXyz):
    raise Exception("_toXmlFromList: problem type '%s' better with _toXmlFromXyz" % str(type(item)))
  
  if not issubclass(item.__class__, list): #simple list
    raise Exception("_toXmlFromList: problem type '%s' unexpected" % str(type(item)))

  res = ET.Element("aList")
  
  if "withTypeClass" in kwargs["styleXml"]:
    res.attrib["typeClass"] = "list"
  if "withTreePyName" in kwargs["styleXml"]:
    res.attrib["treePyName"] = "?"
  
  index = 0
  for i in item:
    ires = toXml(i, **kwargs)
    ires.tag = "_%i_"%index
    ires.attrib["index"] = str(index)
    res.append(ires)
    index += 1
  res.attrib["size"] = str(len(item))
  return res

def _toXmlFromTuple(item, **kwargs):
  """
  | create an ET.Elementtree from item:
  | item is a tuple orsubclasses
  | kwarg are for optional future option of added details in xml tree
  | have to be called from method toXml_
  """

  import xyzpy.baseXyz as BXYZ
  import xyzpy.intFloatListXyz as IFLX

  if issubclass(item.__class__, BXYZ.ListOfBaseXyz):
    raise Exception("_toXmlFromList: problem type '%s' better with _toXmlFromXyz" % str(type(item)))
  
  if not issubclass(item.__class__, tuple): #simple list
    raise Exception("_toXmlFromList: problem type '%s' unexpected" % str(type(item)))

  res = ET.Element("aTuple")
  
  if "withTypeClass" in kwargs["styleXml"]:
    res.attrib["typeClass"] = "tuple"
  if "withTreePyName" in kwargs["styleXml"]:
    res.attrib["treePyName"] = "?"
  
  index = 0
  for i in item:
    ires = toXml(i, **kwargs)
    ires.tag = "_%i_"%index
    ires.attrib["index"] = str(index)
    res.append(ires)
    index += 1
  res.attrib["size"] = str(len(item))
  return res

def _toXmlFromDict(item, **kwargs):
  """
  | create an ET.Elementtree from item:
  | item is a dict orsubclasses
  | kwarg are for optional future option of added details in xml tree
  | have to be called from method toXml_
  """

  import xyzpy.baseXyz as BXYZ
  import xyzpy.intFloatListXyz as IFLX

  if issubclass(item.__class__, BXYZ.BaseFreeXyz):
    raise Exception("_toXmlFromDict: problem type '%' better with _toXmlFromXyz" % str(type(item)))
  
  if not issubclass(item.__class__, dict): #simple dict
    raise Exception("_toXmlFromDict: problem type '%' unexpected" % str(type(item)))

  res = ET.Element("aDict")
  
  if "withTypeClass" in kwargs["styleXml"]:
    res.attrib["typeClass"] = "dict"
  if "withTreePyName" in kwargs["styleXml"]:
    res.attrib["treePyName"] = "?"
  
  for name, value in list(item.items()):
    ires = toXml(value, **kwargs)
    ires.tag = name.replace(" ", "_")
    res.append(ires)

  return res
  
def _createMessageXmlParsing(exception, strXml):
  """add string xml parsed if not too long"""
  errorMessage = str(exception)
  mess = "\nsome errors occurs parsing xml data"
  if len(strXml) < 300:
    mess += "\n"+strXml
  mess += "\n"+errorMessage
  return mess

def assumeTagXml(tagXmlOrStrXml):
  """if string Xml in entry, create and return an ET Xmltree"""
  
  if issubclass(tagXmlOrStrXml.__class__, ET.Element("").__class__):
    return tagXmlOrStrXml
  
  elif issubclass(tagXmlOrStrXml.__class__, str):
    try:
      tagXml = ET.fromstring(tagXmlOrStrXml) #internal deserialization
      return tagXml
    except Exception as e:
      raise Exception(_createMessageXmlParsing(e, tagXmlOrStrXml))
 
  else: #dont know to do with that... try
    try:
      # print("\nET.fromstring %s" % type(tagXmlOrStrXml)
      tagXml = ET.fromstring(tagXmlOrStrXml.decode("utf-8")) #internal deserialization
      return tagXml
    except Exception as e:
      mess = "unknown type data %s" % type(tagXmlOrStrXml)
      raise Exception(_createMessageXmlParsing(e, mess))
 
def fromFileXml(fileName, **kwargs):
  logger.info("load file Xml '%s'" % fileName)
  root = ET.parse(fileName).getroot()
  if verbose: print("fromFileXml readen '%s'" % fileName)
  try:
    res = fromXml(root, **kwargs)
  except:
    traceback.print_exc() #better explicit verbose problem
    raise Exception("Cannot set from file: %s" % fileName)
  return res

def getEtFromFileXml(fileName):
  logger.info("getEtFromFileXml " + fileName)
  try:
    root = ET.parse(fileName).getroot()
  except:
    traceback.print_exc() #better explicit verbose problem
    raise Exception("Cannot set from file: %s" % fileName)
  return root
  
def fromXml(tagXmlOrStrXml, **kwargs):
  if verbose: print("fromXml ",tagXmlOrStrXml)
  #initial method from root tag have not name attribute
  tagXml = assumeTagXml(tagXmlOrStrXml)
  aBaseXyz, attName = _fromXml(tagXml, **kwargs)
  try:
    aBaseXyz.setNameObject(tagXml.tag)
  except:
    pass #if return elementary aDict or aList
  return aBaseXyz
  
def _fromXml(tagXml, **kwargs):
  """current called method in recursivity in parentAsAttribute settatr"""
  from xyzpy.baseXyz import BaseFreeXyz #, ListOfBaseXyz
  from xyzpy.intFloatListXyz import _XyzImmBase
  import xyzpy.classFactoryXyz as CLFX
  
  if verbose: print("utilsXyz._fromXml:'%s' '%s' '%s'" % ( tagXml.tag, tagXml.text, str(tagXml.attrib) ))
  attribs = tagXml.attrib
  #in future dictOfXyzClass be imported"
  dictOfXyzClass = CLFX.getAllXyzClasses() #something like {"BaseFreeXyz": BaseFreeXyz, "ListOfBaseXyz": ListOfBaseXyz}
  
  try:
    nameClass = attribs["typeClass"]
  except:
    logger.warning("inexisting attribute 'typeClass' in tag %s" % tagXml.tag)
    nameClass = "BaseFreeXyz"
  
  if nameClass == "dict":
    return _dictFromXml(tagXml, **kwargs)
  if nameClass == "list":
    return _listFromXml(tagXml, **kwargs)
    
  try:
    typeClass = dictOfXyzClass[nameClass]
  except:
    logger.warning('unknown class in getAllXyzClasses: %s' % nameClass)
    typeClass = BaseFreeXyz
  
  if issubclass(typeClass, _XyzImmBase):   #tagXml.text != "" #tagXml.text__repr__()
    #print "%%%%%%%%%%%%%tagXml.text:", tagXml.text.__repr__()
    aBaseXyz = typeClass(tagXml.text) #init immutable class with value from .tagXml.text
  else: #BaseFreeXyz, etc...
    aBaseXyz = typeClass()
  for childTag in tagXml:
    try:
      aChildBaseXyz, attName = _fromXml(childTag, **kwargs)
      if issubclass(aBaseXyz.__class__, list): #ListOfBaseXyz
        aBaseXyz.append(aChildBaseXyz)
      else:
        aBaseXyz.__setattr__(attName, aChildBaseXyz)
    except Exception as e:
      logger.critical("Problem importing Xml, skipped at your own risk for now.\n%s" % str(e))
  name = tagXml.tag
  aBaseXyz.setIsSet(True)
  return (aBaseXyz, name)

def _listFromXml(tagXml, **kwargs):
  """current called method in recursive with parentAsAttribute settatr"""
  res = []
  for childTag in tagXml:
    aRes, name = _fromXml(childTag, **kwargs)
    res.append(aRes)
  return (res, "aDict")
  
def _dictFromXml(tagXml, **kwargs):
  """current called method in recursive with parentAsAttribute settatr"""
  res = {}
  for childTag in tagXml:
    aRes, name = _fromXml(childTag, **kwargs)
    res[name] = aRes
  return (res, "aList")

def QMessageBoxWarning(parent, mess):
  from PyQt5 import QtWidgets
  QtWidgets.QMessageBox.warning(parent, "warning", mess)

def getDateTimeNow():
  """returns string yymmdd_HHMMSS (151231_235900) based on datetime.now()"""
  now = datetime.now()
  return now.strftime("%y%m%d_%H%M%S")

def getTimeNow():
  """returns string HHMMSS (235900) based on datetime.now()"""
  now = datetime.now()
  return now.strftime("%H%M%S")

def getDateNow():
  """returns string yymmdd (151231) based on datetime.now()"""
  now = datetime.now()
  return now.strftime("%y%m%d")

def newElement(name, text=None, kwargs=None):
  """kwargs are as xml tag attribute"""
  tag =  ET.Element(name)
  #tag.text = UXYZ.toStrForXml(text) #just a little more human readable
  if text != None: tag.text = text #just almost human readable
  if kwargs != None: 
    strkwargs = {}
    for k, v in list(kwargs.items()):
      strkwargs[str(k)] = str(v)
    tag.attrib = strkwargs
  return tag

##############################################################################
# subprocess utilities, with logger functionalities (trace etc.)
##############################################################################

def getPythonVersion():
  """returns python sys.version_info, something like (2, 5, 2, 'final', 0)"""
  import sys
  return sys.version_info

def getPyQtVersion():
  """returns python QtCore.QT_VERSION_STR , something like (5, 9, 1)"""
  core = None
  try:
    import PyQt4.QtCore as core
  except:
    pass
  if core is None:
    try:
      import PyQt5.QtCore as core
    except:
      pass
  if core == None:
    return None
  else:
    return tuple(int(i) for i in core.QT_VERSION_STR.split('.'))

def toString(value):
  """
  assume conversion/decode if value is byte from unicode python2/3
  see https://docs.python.org/3/howto/unicode.html
  """
  typ = type(value)
  if typ is str:
    return value
  if typ is bytes:
    res = value.decode("utf-8", "ignore") # may be latin-1 better
    return res
  logger.warning("unexpected type %s for python2-3 string coding of ''" % (typ, value))
  return value

def exec2or3(aString, namespace=None):
  """python 2-3 exec() compatibility, return a namespace"""
  if verboseExec: logger.warning("exec2or3 '%s'" % aString) # TODO comment this
  if namespace is None:
    namespace = {}
  try:
    exec(aString, namespace)
  except Exception as e:
    msg = "error in exec(%r)\n%s" % (aString, e)
    raise Exception(msg)
  return namespace

def getenv(value):
  """
  replace all environment variable if inexisting (not set)
  with default to avoid None,
  supposedly salome SOURCES/INSTALL configuration
  """
  strValue = str(value)
  res = os.getenv(strValue)
  if res is not None:
    return res

  ROOT_DIR = os.path.realpath(originDir + "/../INSTALL")
  PACKAGESPY_ROOT_DIR = os.path.realpath(originDir + "/../INSTALL/PACKAGESPY")

  if strValue == "PACKAGESPY_ROOT_DIR":
    res = PACKAGESPY_ROOT_DIR

  if strValue == "AMITEXCODE_ROOT_DIR":
    res = ROOT_DIR + "/../INSTALL/AMITEXCODE"

  if strValue == "NUMODISCODEHOME":
    res = ROOT_DIR + "/../INSTALL/NUMODISCODE/INSTALL/NUMODIS"

  if strValue == "URANIESYS":
    res = ROOT_DIR + "/../INSTALL//uranie"

  if res == None:
    logger.critical("Environment variable '%s' not set, fix it or default value." % strValue)
    return None
  else:
    res = os.path.realpath(res)
    # logger.debug("Environment variable '%s' have value '%s'" % (strValue, res)
    return res


def Popen(command, shell=True, cwd=None, env=None, stdout=SP.PIPE, stderr=SP.PIPE, logger=None):
  """
  make subprocess.Popen(cmd), with
  call logger.trace and logger.error if problem as returncode != 0
  returns RCO.ReturnCode
  communicate seems valid under linux
  for windows Popen_Windows is called instead

  | example of use:
  | >>> import as UXYZ
  | >>> rc = UXYZ.Popen(cmd, logger=logger)
  | >>> if rc.isOk():
  | >>>   res_out = rc.getValue()
  | >>>   ...
  | >>> else:
  | >>>   res_out = rc.getWhy()
  | >>>   ...
  """
  cmd = str(command)
  import platform
  if platform.system() == "Windows":
    res = Popen_Windows(cmd, shell=shell, cwd=cwd, env=env, stdout=stdout, stderr=stderr, logger=logger)
    return res

  if True: #try: #TODO
    proc = SP.Popen(cmd, shell=shell, cwd=cwd, env=env, stdout=stdout, stderr=SP.STDOUT)
    res_out, res_err = proc.communicate() # res_err = None as stderr=SP.STDOUT
    rc = proc.returncode
    if verbose: print("Popen end rc = '%s' cmd=\n%s" % (str(rc), PP.pformat(cmd)))
    res_out = res_out.decode("utf-8") # python3 stdout is b'...'
    # DBG.write("Popen logger returncode", (rc, res_out))

    if rc == 0:
      if logger is not None:
        logger.info("<OK> launch command rc=%s cwd=%s:\n%s" % (rc, cwd, cmd))
        logger.debug("<OK> result command stdout&stderr:\n%s" % res_out)
      return RCO.ReturnCode("OK", "Popen command done", value=res_out)
    else:
      if logger is not None:
        logger.warning("<KO> launch command rc=%s cwd=%s:\n%s" % (rc, cwd, cmd))
        logger.warning("<KO> result command stdout&stderr:\n%s" % res_out)
      return RCO.ReturnCode("KO", "Popen command problem", value=res_out)
  else: #except Exception as e: # TODO
    logger.error("<KO> launch command cwd=%s:\n%s" % (cwd, cmd))
    logger.error("launch command exception:\n%s" % e)
    return RCO.ReturnCode("KO", "Popen command problem")


def Popen_Windows(command, shell=True, cwd=None, env=None, stdout=SP.PIPE, stderr=SP.PIPE, logger=None):
  if len(command.split("\n")) > 1:
    logger.critical("TOFIX: Problem windows Popen with multiples command lines:\n%s" % PP.pformat(command))

  if True: #try
    # cmd = "ECHO COUCOU\nSTART /B " + str(command) + "\nEXIT /B %errorlevel%"
    cmd = "START /B " + str(command)
    if verbose: print("Popen_Windows begin cmd = '%s'" % cmd)
    proc = SP.Popen(cmd, shell=shell, cwd=cwd, env=env, stdout=stdout, stderr=SP.STDOUT)
    proc.wait() # this freeze main windows if not 'START /B'
    rc = proc.returncode
    if verbose: print("Popen_Windows end rc = '%s'" % str(rc))
    res_out = "No stdout under windows"
    if rc == 0:
      if logger is not None:
        logger.info("<OK> launch command rc=%s cwd=%s:\n%s" % (rc, cwd, cmd))
      return RCO.ReturnCode("OK", "Popen_Windows command done", value=res_out)
    else:
      if logger is not None:
        logger.warning("<KO> launch command rc=%s cwd=%s:\n%s" % (rc, cwd, cmd))
      return RCO.ReturnCode("KO", "Popen_Windows command problem", value=res_out)
  else: #except Exception as e:
    logger.error("<KO> launch command cwd=%s:\n%s" % (cwd, cmd))
    logger.error("launch command exception:\n%s" % e)
    return RCO.ReturnCode("KO", "Popen command Exception")


def getFirstCharacterInFile(aFile):
  if not os.path.isfile(aFile):
    logger.warning("inexisting file '%s'" % aFile)
    return None
  with open(aFile, "r") as f:
    return f.read(1)
