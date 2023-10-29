#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


__doc__ = """\
class to assume view of directories
TODO NOT IMPLEMENTED YET
"""

import os
from xyzpy.baseXyz import _XyzConstrainBase
import xyzpy.classFactoryXyz as CLFX
import xyzpy.utilsXyz as UXYZ
import salomepy.utilsWorkdir as UTW
import filewatcherpy.fileWatcher as FWR

import xyzpy.loggingXyz as LOG
logger = LOG.getLogger()

verbose = False


###############################################################
class DirectoryManagerXyz(_XyzConstrainBase):
  """
  general informations about
  """
  _attributesList = [ #list, not a dict because sequential order list is used in files
    ("directoryName", "StrXyz"),
    ("displayDirectories", "BoolXyz"),
    ("displayFiles", "BoolXyz"),
    ("displayDepth", "IntPosXyz"),
    #("directories", "ListOfDirectories"),
    #("files", "ListOfFiles"),
  ]

  _helpDict = {
    "directoryName": ("Name of directory (, with or without environ variable", ""),
  }

  def __init__(self):
    super(DirectoryManagerXyz, self).__init__()
    self.setIsCast(True)

  def setDefaultValues(self):
    self.directoryName = os.getcwd()
  
  def getRealPath(self, aPathWithEnvVar):
    """
    resolve env variable as $HOME/toto etc... 
    with subprocess shell interpretation of env var
    """
    logger.warning("TODO: DirectoryManagerXyz NOT IMPLEMENTED YET, see iradinaGUI version")
    return FWR.getRealPath(aPathWithEnvVar)



#factory pattern using xyzpy.utilsXyz._dictOfXyzClass
CLFX.appendAllXyzClasses( [DirectoryManagerXyz] )
