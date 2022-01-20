#!/usr/bin/env python

from xyzpy.baseXyz import _XyzConstrainBase
import xyzpy.classFactoryXyz as CLFX

import solverlabpy.equationSvl as EQ
from xyzpy.intFloatListXyz import StrInListXyz, XyzQComboBox, StrNoEditionXyz

from PyQt5 import QtWidgets as QTW

class MandatoryTransportEq(_XyzConstrainBase):
  _attributesList = [
    ("transport_velocity", ""),
    ("rod_temp_field", ""),
    ("fluid_temp_field", ""),
    ("", ""),
  ]
