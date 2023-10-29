#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import re
from PyQt5 import QtCore, QtGui, QtWidgets
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()

# Regular expression to find floats. Match groups are the whole string, the
# whole coefficient, the decimal part of the coefficient, and the exponent
# part.
_float_re = re.compile(r'(([+-]?\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)')
_AMinFloat = -1e100
_AMaxFloat = +1e100

verbose = False

def valid_float_string(string):
  match = _float_re.search(string)
  return match.groups()[0] == string if match else False

class FloatValidator(QtGui.QValidator):
  def init(self, parent = None):
    super(FloatValidator, self).__init__(parent = parent)

  def validate(self, string, position):
    if valid_float_string(string):
      return (self.Acceptable, string, position)
    if string == "" or str(string[position-1]) in 'e.-+':
      return (self.Intermediate, string, position)
    return (self.Invalid, string, position)
     
  def fixup(self, text):
    match = _float_re.search(text)
    return match.groups()[0] if match else ""

class SciQDoubleSpinBox(QtWidgets.QDoubleSpinBox):
 
  def __init__(self, *args, **kwargs):
    self.verbose = False
    self._finalValue = 0.
    self._initialValue = 0.
    super(SciQDoubleSpinBox, self).__init__(*args, **kwargs)
    self.setMinimum(_AMinFloat)
    self.setMaximum(_AMaxFloat)
    self.validator = FloatValidator(self)
    self.setDecimals(12)
    self.setValue(self._finalValue)

  """def validate(self, text, position):
    #test pyqt5 validate
    res = super(SciQDoubleSpinBox, self).validate(text, position)
    print "!!!validate Qt5", res #!!!validate Qt5 (2, u'1.2', 0)
    return res"""
    
  def validate(self, text, position):
    #http://doc.qt.io/qt-5/qabstractspinbox.html#validate
    #QValidator::State QAbstractSpinBox::validate(QString &input, int &pos)
    #print("validate text '%s' position '%s'" % (text , position))
    res = self.validator.validate(text, position)
    if self.verbose: 
      logger.info("validate res '%s' text '%s'" % (res , text))
    #QValidator::State 3 items
    #!!!validate Qt5 (2, u'1.2', 0)
    return res
   
  def fixup(self, text):
    res = self.validator.fixup(text)
    if verbose: logger.info("fixup %s %s" % (res , text))
    self.setValue(res)
    return res
   
  def valueFromText(self, text):
    res = float(text)
    if verbose: logger.info("valueFromText %s %s %s" % (res, self.minimum(), self.maximum()))
    if res < self.minimum(): res = self.minimum()
    if res > self.maximum(): res = self.maximum()
    self._finalValue = res
    return res
   
  def textFromValue(self, value):
    if verbose: logger.info("textFromValue %s %s" % (value, self._finalValue))
    #res = format_float(value)
    res = format_float(self._finalValue)
    return str(res)
   
  def stepBy(self, steps):
    text = self.cleanText()
    groups = _float_re.search(text).groups()
    decimal = float(groups[1])
    decimal += steps
    new_string = "{0:g}".format(decimal) + (groups[3] if groups[3] else "")
    if verbose: print("stepBy1", steps , new_string)
    
    res = float(new_string)
    if res < self.minimum(): decimal = self.minimum()
    if res > self.maximum(): decimal = self.maximum()
    new_string = "{0:g}".format(decimal) + (groups[3] if groups[3] else "")
    if verbose: print("stepBy2", steps , new_string)
    
    self.lineEdit().setText(new_string)
    
  def setValue(self,  value):
    if verbose: print("setValue", value, type(value))
    fval = float(value)
    new_string = format_float(fval)
    self._initialValue = float(value)
    self.lineEdit().setText(new_string)
    return super(SciQDoubleSpinBox, self).setValue(float(value))
    
  def value(self):
    """return string, not float"""
    if verbose: print("value", self._finalValue)
    return format_float(self._finalValue)

  def getValue(self):
    return self._finalValue


def format_float(value):
  """Modified form of the 'g' format specifier."""
  string = "{0:.7g}".format(value).replace("e+", "e")
  #string = "{0:.7e}".format(value).replace("e+", "e")
  string = re.sub("e(-?)0*(\d+)", r"e\1\2", string)
  #print "format_float", value , string
  return string
