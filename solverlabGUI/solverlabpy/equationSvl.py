#!/usr/bin/env python

from xyzpy.baseXyz import _XyzConstrainBase, ListOfBaseXyz
import xyzpy.classFactoryXyz as CLFX
from xyzpy.intFloatListXyz import StrInListXyz, IntRangeXyz


class Int20Xyz(IntRangeXyz):
  _allowedRange = [-2, 0]
  _defaultValue = _allowedRange[0]


class FieldOptionSvl(_XyzConstrainBase):
  """
  Generic field selector with advanced settings
  """
  _attributesList = [
    ("field_name", "SelectField"),
    ("time_iteration", "IntPosXyz"),
    ("time_sub_iteration", "IntPosXyz"),
    ("meshLevel", "Int20Xyz"),
  ]
  _helpDict = {
    "field_name": ("Name of the field created by solverlab", ""),
    "time_iteration": ("?", ""),
    "time_sub_iteration": ("?", ""),
    "meshLevel": ("Dimension of the field compared to the dim of the mesh", ""),
  }

  def __init__(self):
    super(FieldOptionSvl, self).__init__()
    self.setIsCast(True)
    self.setDefaultValues()

  def setDefaultValues(self):
    self.time_iteration = 0
    self.time_sub_iteration = 0
    self.meshLevel = 0


class MethodList(StrInListXyz):
  _allowedList = ["Explicit", "Implicit"]
  _defaultValue = _allowedList[0]


class LinearSolverList(StrInListXyz):
  _allowedList = ["GMRES", "BICGSTAB"]
  _defaultValue = _allowedList[0]


class PreconditionerList(StrInListXyz):
  _allowedList = ["LU", "ILU", "None"]
  _defaultValue = _allowedList[0]


class FormatList(StrInListXyz):
  _allowedList = ["VTK", "MED"]
  _defaultValue = _allowedList[0]


class ComputationSvl(_XyzConstrainBase):
  """
  Generic computation parameters needed to launch any solverlab calculus
  """
  _attributesList = [
    ("max_time_step", "IntXyz"),
    ("save_frequency", "IntXyz"),
    ("max_time", "IntXyz"),
    ("cfl", "FloatXyz"),
    ("precision", "FloatXyz"),
    ("file_name", "StrXyz"),
    ("file_format", "FormatList")
  ]
  _helpDict = {
    "max_time_step": ("Number max of time step\ndefault value is 10", ""),
    "save_frequency": ("Number of time step between two save\ndefault value is 1", ""),
    "max_time": ("Maximum time for the solverlab calculus\ndefault value is 10", ""),
    "cfl": ("convergence condition (Courant–Friedrichs–Lewy)\ndefault value is 1", "https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition"),
    "precision": ("default value is 1e-6", ""),
    "file_name": ("Name of the result file", ""),
    "file_format": ("Format of the result file", ""),
  }

  def __init__(self):
    super(ComputationSvl, self).__init__()
    self.setIsCast(True)
    self._setAllAttributesList()
    self.setDefaultValues()

  def setDefaultValues(self):
    self.max_time_step = 10
    self.save_frequency = 1
    self.max_time = 10
    self.cfl = 1
    self.precision = 0.000001
    self.file_name = "result"


class NumericalSvl(_XyzConstrainBase):
  """
  Generic numerical parameters for all solverlab calculus
  """
  _attributesList = [
    ("method", "MethodList"),
    ("linear_solver", "LinearSolverList"),
    ("preconditioner", "PreconditionerList"),
    ("max_iteration", "IntSupEq1Xyz")
  ]
  _helpDict = {
    "method": ("type of method", ""),
    "linear_solver": ("method use to solve linear systems", ""),
    "preconditioner": ("preconditioner for the linear_solver", ""),
  }

  def __init__(self):
    super(NumericalSvl, self).__init__()
    self.setIsCast(True)
    self._setAllAttributesList()


def getListField(file):
  """Return list of all the field name in file (need to be a .med)"""
  import MEDLoader as ML
  data = ML.MEDFileData(file)
  alist = data.getFields()
  rep = [f.getName() for f in alist]
  return rep


class ListOfEquation(ListOfBaseXyz):
  """List of all solverlab model implemented"""
  _allowedClasses = ["DiffusionEq"]

  def getActionsContextMenu(self):
    """append action 'Append equation'"""
    actions = super(ListOfEquation, self).getActionsContextMenu()
    if self.getRoot().getMedfile().endswith(".med"):
      return actions
    return actions[:3]


CLFX.appendAllXyzClasses([Int20Xyz, FieldOptionSvl,
                          MethodList, LinearSolverList, PreconditionerList, FormatList,
                          ComputationSvl, NumericalSvl, ListOfEquation, ])
