#!/usr/bin/env python

from xyzpy.baseXyz import _XyzConstrainBase
import xyzpy.classFactoryXyz as CLFX

import solverlabpy.equationSvl as EQ
from xyzpy.intFloatListXyz import StrInListXyz, XyzQComboBox, StrNoEditionXyz
from xyzpy.baseXyz import ListOfBaseXyz

from PyQt5 import QtWidgets as QTW

class CalculationMethod(StrInListXyz):
  _allowedList = ["Element", "Volume"]
  _defaultValue = _allowedList[0]

class MandatoryDiffusionEq(_XyzConstrainBase):
  """
  specific value mandatory to launch Diffusion Equation
  """
  _attributesList = [
    ("specific_heat", "FloatXyz"),  # specific heat
    ("solid_density", "FloatXyz"),  # solid density
    ("thermal_conductivity", "FloatXyz"),  # thermal conductivity
  ]
  _helpDict = {
    "specific_heat": ("heat capacity (J/Kg/K)", ""),
    "solid_density": ("density (Kg/m^dim)", ""),
    "thermal_conductivity": ("conductivity (W/m/k)", ""),
  }

  def __init__(self):
    super(MandatoryDiffusionEq, self).__init__()
    self.setIsCast(True)
    # self._setAllAttributesList()
    self.setDefaultValues()

  def setDefaultValues(self):
    self.specific_heat = 1
    self.solid_density = 1
    self.thermal_conductivity = 0


class ModeFieldSvl(StrInListXyz):
  _allowedList = ["Scalar", "Field"]
  _defaultValue = _allowedList[0]

class OptionalDiffusionEq(_XyzConstrainBase):
  """
  specific optional value for the diffusion equation
  there are all set at default value
  """
  _attributesList = [
    ("mode_fluid_temp", "ModeFieldSvl"),
    ("fluid_temp_scalar", "FloatXyz"),
    ("fluid_temp_field", "FieldOptionSvl"),
    ("heat_transfer_coeff", "FloatXyz"),
    ("mode_heat_power", "ModeFieldSvl"),
    ("heat_power_scalar", "FloatXyz"),
    ("heat_power_field", "FieldOptionSvl"),
  ]
  _helpDict = {
    "mode_fluid_temp": ("Choice for fluid temperature simple or advanced", ""),
    "fluid_temp_scalar": ("fluid mean temperature (K)\ndefault value is 0", ""),
    "fluid_temp_field": ("Select fluid temperature field from med", ""),
    "heat_transfer_coeff": ("fluid/solid exchange coefficient\ndefault value is 0", ""),
    "mode_heat_power": ("Choice for heat power simple or advanced", ""),
    "heat_power_scalar": ("heat power density (J/m^dim)\ndefault value is 0", ""),
    "heat_power_field": ("Select heat power field from med", "")
  }

  def __init__(self):
    super(OptionalDiffusionEq, self).__init__()
    self.setIsCast(True)
    self._setAllAttributesList()
    self.setDefaultValues()

  def setDefaultValues(self):
    self.fluid_temp_scalar = 0
    self.heat_transfer_coeff = 0
    self.heat_power_scalar = 0
    self.fluid_temp_field.field_name = ""
    self.heat_power_field.field_name = ""

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    # print("Optional diffusion: %s.isHidden" % (nameAttr))
    if nameAttr == "fluid_temp_scalar" and self.mode_fluid_temp == "Field":
        return True
    if nameAttr == "fluid_temp_field" and self.mode_fluid_temp == "Scalar":
        return True
    if nameAttr == "heat_power_scalar" and self.mode_heat_power == "Field":
        return True
    if nameAttr == "heat_power_field" and self.mode_heat_power == "Scalar":
        return True
    return False


class SelectBondary(StrInListXyz):
  _allowedList = ["Neumann", "Dirichlet"]
  _defaultValue = _allowedList[0]


class OneBorder(_XyzConstrainBase):
  _attributesList = [
    ("name", "StrNoEditionXyz"),
    ("type", "SelectBondary"),
  ]

  def __init__(self):
    super(OneBorder, self).__init__()
    self.setIsCast(True)
    self._setAllAttributesList()


class BoundaryDiffusionEq(ListOfBaseXyz):  #TODO read boundary condition in .med
  _allowedList = [OneBorder]

  def __init__(self):
    super(BoundaryDiffusionEq, self).__init__()

  def update(self):
    try:
      import MEDLoader as ML
    except Exception as e:
      QTW.QMessageBox.warning(self.getController().getDesktop(), "warning", "problem import MEDLoader : %s" % e)
      return

    file = self.getRoot().getMedfile()
    data = ML.MEDFileData(file)
    meshes = data.getMeshes()

    for i, m in enumerate(meshes):
      names = m.getGroupsNames()
      if len(names) == 0:
        continue
      else:
        print(names)
        for n in names:
          border = OneBorder()
          border.name = n
          self.append(border)


class SelectField(StrNoEditionXyz):
  def createEditor(self, parent):
    file = self.getRoot().getMedfile()
    fieldList = EQ.getListField(file)
    if len(fieldList) == 0:
      QTW.QMessageBox.warning(self.getController().getDesktop(), "warning", "No field found in\n'%s'" % file)
      return None
    combo = XyzQComboBox(parent)
    combo.addItems(fieldList)
    combo.setCurrentIndex(0)
    return combo

  def getActionsContextMenu(self):
    self.hidden = not self.hidden
    actions = super(SelectField, self).getActionsContextMenu()
    return actions


class DiffusionEq(_XyzConstrainBase):
  """
  Entry for the diffusion equation model
  """
  _attributesList = [
    ("field_name", "SelectField"),
    ("field_option", "FieldOptionSvl"),
    ("space_dim", "Int13Xyz"),
    ("time_iteration", "IntXyz"),  # TODO verifier les incohérences
    ("finite_method", "CalculationMethod"),
    ("mandatory_values", "MandatoryDiffusionEq"),
    ("optional_values", "OptionalDiffusionEq"),
    ("boundary_condition", "SelectBondary"),
    ("numerical_parameters", "NumericalSvl"),
    ("computation_parameters", "ComputationSvl")
  ]
  _helpDict = {
    "field_name": ("The field to use on your .med", ""),
    "space_dim": ("Space dimension", ""),
    "time_iteration": ("Time iteration on field (med functionality)", ""),
    "finite_method": ("Finite method used for the calculus", ""),
    "mandatory_values": ("Information needed to launch the simulation", ""),
    "optional_values": ("Right click to add optional physical values", ""),
    "boundary_condition": ("Boundary condition", ""),
    "numerical_parameters": ("Method parameters", ""),
    "computation_parameters": ("Parameters for the computation", ""),
  }

  def __init__(self):
    super(DiffusionEq, self).__init__()
    self.setIsCast(True)
    self._setAllAttributesList()
    # self.setDefaultValues()

  def setDefaultValues(self):
    # file = self.getRoot().getMedfile()
    # self.fieldName._allowedList = EQ.getListField(file)
    # self.spaceDim._allowedRange = [1,3]
    self.field_name = ""
    self.time_iteration = 0
    #self.mandatory_values = MandatoryDiffusionEq()
    #self.optional_values = OptionalDiffusionEq()
    #self.numerical_parameters = EQ.NumericalSvl()


CLFX.appendAllXyzClasses([MandatoryDiffusionEq, SelectField, ModeFieldSvl, OptionalDiffusionEq, CalculationMethod,
                          SelectBondary, OneBorder, BoundaryDiffusionEq, DiffusionEq, ])
