#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


__doc__ = """\
class to assume simulations cases WORKDIR set, and create directories
"""

import os
from datetime import datetime
from xyzpy.baseXyz import _XyzConstrainBase
import xyzpy.classFactoryXyz as CLFX
import xyzpy.utilsXyz as UXYZ
import salomepy.utilsWorkdir as UTW

import xyzpy.loggingXyz as LOG
logger = LOG.getLogger()

verbose = False


########################################################################################
class DataInformationsXyz(_XyzConstrainBase):
  """
  general general informations about cases
  """
  _attributesList = [ #list, not a dict because sequential order list is used in files Cis
    ("name", "StrXyz"),
    #("fileXmlOrigin", "FileXyz"),
    ("directory", "DirectoryXyz"),
    ("release", "ReleaseXyz"),
    #("dateFile", "DateXyz"),
  ]

  _helpDict = {
    "name": ("As name of case directory", ""),
    "directory": ("As parent directory of case directory", ""),
    "release": ("Release of files", ""),
    #"dateFile": (u"Date of last modification from file system", u""),
  }
  _icon = "datainformation"
  _defaultSimName = "Simulation"
  _defaultWorkdirEnvVar = ""

  def __init__(self):
    super(DataInformationsXyz, self).__init__()
    self.setIsCast(True)
    self.setDefaultValues()
    
  def setFromStream(self, f, verbose=False):
    if verbose: print("DataInformationsXyz.setFromStream")
    fileName = f.name
    v = verbose
    try: ##########################
      aDir, aName = os.path.split(os.path.realpath(fileName))
      aDir, aName = os.path.split(aDir)
      self.name = aName
      self.directory = aDir
      #self.release = "1.5"
      #print "last modified: %s" % time.ctime(os.path.getmtime(fileName))
      #self.dateFile = time.ctime(os.path.getmtime(fileName))
      #print "last modified: %s" % time.ctime(os.path.getmtime(file))

    except:
      traceback.print_exc() #better explicit verbose problem
      logger.error(_loc("Cannot read file: ") + fileName)
      raise Exception("Cannot read file: " + fileName)
    self.setIsSet(True)
    return 

  def setFromFile(self, fileName, verbose=False):
    #TODO test name of file
    f = open(fileName, 'r')
    self.setFromStream(f, verbose)

  def toStr(self):
    lf = self._lf
    res = ""
    res += "# name"
    res += UXYZ.toStrForXml(self.name)+lf
    res += "# directory"
    res += UXYZ.toStrForXml(self.directory)+lf
    #res += "# release"
    #res += UXYZ.toStrForXml(self.release)+lf
    #res += "# dateFile"
    #res += UXYZ.toStrForXml(self.dateFile)+lf
    return res

  def setDefaultValues(self):
    aDirw = UTW.getWorkdirDefault(self._defaultWorkdirEnvVar)
    aDirName = os.path.join(aDirw, self._defaultSimName + UXYZ.getDateTimeNow())
    aDir, aName = os.path.split(os.path.realpath(aDirName))
    self.name = aName
    self.directory = aDir
    self.release = "0.0.0"
    #print "last modified: %s" % time.ctime(os.path.getmtime(fileName))
    #self.dateFile = "" #no file saved time.ctime(os.path.getmtime(fileName))

  def getLaunchDir(self):
    res = os.path.join(self.directory, self.name)
    res = os.path.realpath(os.path.expandvars(res))
    return res

  def createLaunchDir(self):
    UTW.makeDir(self.getLaunchDir())


#factory pattern using xyzpy.utilsXyz._dictOfXyzClass
CLFX.appendAllXyzClasses( [DataInformationsXyz] )
