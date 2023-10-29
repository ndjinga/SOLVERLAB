#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

"""
factory helps for tooltips
define and store commons help(s)
for multiple inherited classes of baseXyz.
one global factory singleton
goal is to EZ retrive an set tootip and whatis
from everywhere, and from all widgets tree or else views of MVC pattern
"""

import pprint as PP
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()
verbose = False

__commonHelps__ = {}
__doneOnce__ = [False]


def getCommonHelp(nameClass, nameAttribute):
  
  import xyzpy.intFloatListXyz as IFLX
  import xyzpy.baseXyz as BXYZ

  try:
    aDict = __commonHelps__[nameClass]
    anHelp = aDict[nameAttribute]
    res = BXYZ.HelpXyz()
    res.shortHelp = IFLX.StrXyz(anHelp[0])
    res.longHelp = IFLX.StrXyz(anHelp[1])
    return res
  except Exception as e:
    if verbose: 
      logger.warning("getCommonHelp problem %s.%s:\n'%s'" % (nameClass, nameAttribute, e))
    return None

def getCommonToolTip(nameClass, nameAttribute):
  res = getCommonHelp(nameClass, nameAttribute)
  if res is not  None:
    res = res.shortHelp
  if verbose:
    logger.warning("getCommonToolTip %s.%s is '%s'" % (nameClass, nameAttribute, res))
  return res

    
"""
#do not work, no parent information for tagXml
def getCommonToolTipFromTagXml(tagXml):
  #print "CommonHelps:\n", toStrXml()
  nameAttribute = tagXml.tag
  print "tagXml", nameAttribute, tagXml.attrib
  nameClass = "CrescendoCnf" #TODO name of parent typeClass
  res = getCommonToolTip(nameClass, nameAttribute)
  if res != None: print ("tootip", tagXml.tag, res)
  return res
"""

def toXml():
  import xyzpy.utilsXyz as UXYZ
  return UXYZ.toXml(__commonHelps__)

def toStrXml():
  import xyzpy.utilsXyz as UXYZ
  return UXYZ.prettyPrintET(UXYZ.toXml(__commonHelps__, styleXml="withTypeClass,withoutTreePyName"))

def addInCommonHelp(anXyzConstrainBase):
  name = anXyzConstrainBase.__class__.__name__
  
  if name in list(__commonHelps__.keys()):
    logger.warning("helpDict of class '%s' set yet" % name)
    return False
    
  if not hasattr(anXyzConstrainBase, "_helpDict"):
    # is normal for elementaries class
    # cvw TODO logger.debug("helpDict inexistent in class '%s'" % name)
    return False
  
  if verbose:
    logger.warning("addInCommonHelp for class %s:\n%s" % (name, PP.pformat(anXyzConstrainBase._helpDict)))
  __commonHelps__[name] = anXyzConstrainBase._helpDict
  return True
   
