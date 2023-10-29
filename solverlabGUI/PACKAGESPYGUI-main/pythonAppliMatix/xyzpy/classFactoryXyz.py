#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""\
| factory XyzClasses of xyzpy:
| 
| - define and store commons inherited classes of baseXyz etc...
"""

import traceback
import pprint as PP
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()

__dictOfXyzClass__ = {}

verbose = False

def appendAllXyzClasses(localsDictOrlistOfclasses, raiseIfProblem=False):
  """
  | factory pattern using xyzpy.utilsXyz.__dictOfXyzClass__ :
  |
  | - tricky way to get access to (future) user defined classes, in others packages.
  |   but it works.
  | - with one more line in end of future (...other packages) user files
  |   userDefinedMyBaseXyz.py and userDefinedMyIntFloatListXyz.py
  |   (please use this sort of name)
  | 
  | usage:
  | >>> #simply add this ended line:
  | >>> import xyzpy.classFactoryXyz as CLFX
  | >>> CLFX.appendAllXyzClasses( locals() )
  | >>> #or
  | >>> CLFX.appendAllXyzClasses( [oneUserXyzClass, anotherOneUserXyzClass,...] )
  | >>> #and so associated with
  | >>> CLFX.getAllXyzClasses()
  """
  
  import xyzpy.helpsFactoryXyz as HLFX
  currentlistOfclasses = None
  if isinstance(localsDictOrlistOfclasses, dict):
    keys = [key for key in list(localsDictOrlistOfclasses.keys()) if key[-3:] == "Xyz"]
    currentlistOfclasses = [ localsDictOrlistOfclasses[key] for key in keys]
  elif isinstance(localsDictOrlistOfclasses, list):
    currentlistOfclasses = localsDictOrlistOfclasses
  else:
    raise Exception("unexpected type for parameter localsDictOrlistOfclasses (have to be dict or list): " + str(type(localsDictOrlistOfclasses)))
  for iclass in currentlistOfclasses:
    name = iclass.__name__
    if name in list(__dictOfXyzClass__.keys()):
      # have to fix it in program
      # could be present if two files from forks as ElementDrt in PACKAGESPY AND ALSO IRADINAGUI
      msg = "appendXyzClass class known yet: '%s'" % name
      if raiseIfProblem:
        raise Exception(msg)
      else:
        # continue, hope it is not for crash
        logger.warning(msg)
    else:
      if verbose: logger.info("appendAllXyzClasses %s " % name)
      __dictOfXyzClass__[name]=iclass
      #print "appendAllXyzClasses",name, iclass
      try:
        HLFX.addInCommonHelp(iclass())
      except:
        traceback.print_exc()
        raise Exception("problem for class " + name)
  return

def getAllXyzClasses():
  """
  | factory pattern using xyzpy.utilsXyz.__dictOfXyzClass__
  | and so we get a set of All current trans-package ClassesXyz
  """
  if __dictOfXyzClass__ == {}:
    logger.warning("no class in getAllXyzClasses try 'import intFloatListXyz' for a minimum")
    try:
      import xyzpy.intFloatListXyz as IFLX #append factory classes
    except:
      raise Exception("no class in getAllXyzClasses even trying 'import intFloatListXyz'")
  return __dictOfXyzClass__

def getXyzClassFromName(nameClass):
  """
  return a class from his name string
  
  | usage:
  | >>> import xyzpy.classFactoryXyz as CLFX
  | >>> aFloatClass = CFLX.getXyzClassFromName("FloatXyz")
  | >>> anInstance = aFloatClass(1.234)
  """
  #something like {"BaseFreeXyz": BaseFreeXyz, "ListOfBaseXyz": ListOfBaseXyz}
  dictOfXyzClass = getAllXyzClasses() #something like {"BaseFreeXyz": BaseFreeXyz, "ListOfBaseXyz": ListOfBaseXyz}
  if type(nameClass) == dict: #as an attrib from xml
    try:
      aNameClass = nameClass["typeClass"]
    except:
      logger.warning('unknown class in nameClass: %s' % str(nameClass))
      return None
  else: #as a name str
    aNameClass = nameClass
  try:
    typeClass = dictOfXyzClass[aNameClass]
  except:
    logger.warning("unknown class in getAllXyzClasses: '%s', known are:\n%s" % \
                    (aNameClass, PP.pformat(sorted(dictOfXyzClass.keys()))))
    typeClass = None
  return typeClass
  
def getXyzInstanceClassFromNameAndValue(nameClass, value):
  import xyzpy.intFloatListXyz as IFLX
  aClass = getXyzClassFromName(nameClass)
  acceptedClass = [IFLX._XyzImmBase,  str,  int,  int,  float]
  for c in acceptedClass:
    if issubclass(value.__class__, c): #immutables only
      try:
        anInstance = aClass(value)
        return anInstance
      except Exception as e:
        logger.error("Uncastable class:\n%s" % e)
        raise Exception("Uncastable class '%s' to class '%s' with immutable value '%s'\n" % \
                        (value.__class__.__name__, nameClass, value))

  #ListOfBaseXyz
  ListOfBaseXyz = getXyzClassFromName("ListOfBaseXyz")
  if verbose: 
    logger.info("issubclass(value) %s %s" % \
                (issubclass(value.__class__, ListOfBaseXyz), issubclass(aClass, ListOfBaseXyz)))
  if issubclass(value.__class__, ListOfBaseXyz) and issubclass(aClass, ListOfBaseXyz):
    if True: # TODO try:
      anInstance = aClass(value)
    else: # TODO except:
      raise Exception("Uncastable ListOf class '%s' to class '%s' with value '%s'\n" % \
                      (value.__class__.__name__, nameClass, str(value)))
  
  #other may be possible, if specified in aClass.__init__(value)
  try:
    anInstance = aClass(value)
    return anInstance
  except Exception as e:
    logger.error("Uncastable class:\n%s" % e)
    raise Exception("Uncastable class '%s' to class '%s' with value '%s'\n" % \
                    (value.__class__.__name__, nameClass, str(value)))

def toClass(aClass):
  """Return a Class for modelXyz if aClass is type string"""
  if type(aClass) is str:
    return getXyzClassFromName(aClass)
  # supposed class
  return aClass

def toClassList(aList):
  """Return a list of Class for modelXyz from list like _allowedClasses"""
  if type(aList) is not list:
    raise Exception("toClassList expect a list instead '%s'" % type(aList))
  return [toClass(i) for i in aList]

def toString(aClass):
  """Return a string for modelXyz if aClass is type class"""
  if type(aClass) is str:
    return aClass
  else:
    return aClass.__name__

def toStringList(aList):
  """Return a list of string for modelXyz from list like _allowedClasses"""
  if type(aList) is not list:
    raise Exception("toStringList expect a list instead '%s'" % type(aList))
  return [toString(i) for i in aList]