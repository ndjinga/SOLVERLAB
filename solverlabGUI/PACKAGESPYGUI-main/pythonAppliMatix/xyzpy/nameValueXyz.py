#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


__doc__ = """\
class to assume expression 'name=value'
"""

import os
from xyzpy.baseXyz import _XyzConstrainBase
import xyzpy.classFactoryXyz as CLFX

###############################################################
class NameValueXyz(_XyzConstrainBase):
  """
  define expression 'name=value'
  """
  _attributesList = [ #list, not a dict because sequential order list is used in files
    ("name", "StrXyz"),
    ("value", "StrXyz"),
  ]

  _helpDict = {
    "name": ("left-hand term 'name' of expression as 'name=value'", ""),
    "value": ("right-hand term 'value' of expression as 'name=value'", ""),
  }
  _defaultNameValue = "NoName=NoValue"

  def __init__(self):
    super(NameValueXyz, self).__init__()
    self.setIsCast(True)
    self.setDefaultValues()

  def setDefaultValues(self):
    name, value = self._defaultNameValue.split("=")
    self.name = name
    self.value = value

  def getValue(self):
    return ("%s=%s" % (self.name.strip(), self.value.strip()))

#factory pattern using xyzpy.utilsXyz._dictOfXyzClass
CLFX.appendAllXyzClasses( [NameValueXyz] )
