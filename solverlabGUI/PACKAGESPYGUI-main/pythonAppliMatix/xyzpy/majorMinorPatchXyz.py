#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


__doc__ = """\
class to assume version X.Y.Z
"""

import os
from xyzpy.baseXyz import _XyzConstrainBase
import xyzpy.classFactoryXyz as CLFX

###############################################################
class MajorMinorPatchXyz(_XyzConstrainBase):
  """
  general informations about version X.Y.Z 
  used only for utility widget dialog from ReleaseXyz
  """
  _attributesList = [ #list, not a dict because sequential order list is used in files
    ("major", "IntPosXyz"),
    ("minor", "IntPosXyz"),
    ("patch", "IntPosXyz"),
  ]

  _helpDict = {
    "major": ("Major release (major.minor.patch)", ""),
    "minor": ("Minor release (major.minor.patch)", ""),
    "patch": ("Patch release (major.minor.patch)", ""),
  }

  _defaultVersion = "0.0.0"

  def __init__(self):
    super(MajorMinorPatchXyz, self).__init__()
    self.setIsCast(True)

  def setDefaultValues(self):
    major, minor, patch = self._defaultVersion.split(".")
    self.major = major
    self.minor = minor
    self.patch = patch

  def getVersion(self):
    return ("%i.%i.%i" %(self.major, self.minor, self.patch))


#factory pattern using xyzpy.utilsXyz._dictOfXyzClass
CLFX.appendAllXyzClasses( [MajorMinorPatchXyz] )
