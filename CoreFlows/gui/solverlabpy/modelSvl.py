#!/usr/bin/env python

import pprint as PP
import os

import xyzpy.classFactoryXyz as CLFX
from xyzpy.baseXyz import _XyzConstrainBase, ListOfBaseXyz
import xyzpy.dataFromFileXyz
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()

import solverlabpy.analysisSvl
import solverlabpy.equationSvl as EQ
import solverlabpy.diffusionEq


class MedSvl(_XyzConstrainBase):
  """
  Class containing the path to med file and list of equation to apply on it
  """
  _attributesList = [
    ("GeometryMed", "DataFromFileMedXyz"),
    ("equation", "ListOfEquation"),
  ]
  _helpDict = {
    "GeometryMed": ("The file .med for initial geometry and fields etc.", ""),
    "equation": ("Right click to choose the equation you want to apply on the .med", ""),
  }

  def __init__(self):
    super(MedSvl, self).__init__()
    self.setIsCast(True)
    self._setAllAttributesList()
    self.setDefaultValues()

  def setDefaultValues(self):
    # self.FileMed = "Right click to browse file"
    # self.equation = EQ.EquationGeneralSvl()
    return


class ListOfMedSvl(ListOfBaseXyz):
  _allowedClasses = [MedSvl]

  def getActionsContextMenu(self):
    """append action 'Append file med'"""
    actions = super(ListOfMedSvl, self).getActionsContextMenu()
    # actions.append( self._createAction('Append file materiaux', None, 'Browse looking for file materiaux.xml', self.browseDialog, "materiauxamt") )
    return actions

  def browseDialog(self):
    for i in self._allowedClasses:
      self.addItem(i)
    f = self[-1].fileMateriaux
    f.browseDialog()
    return


class ModelSvl(_XyzConstrainBase):
  """
  General model for solverlab
  """
  _icon = "solverlabpy.resources.solverlabgui"
  _attributesList = [
    ("GeometryMed", "DataFromFileMedXyz"),
    ("Model", "ListOfEquation"),
    ("Analysis", "AnalysisSvl"),
  ]
  _helpDict = {
    "GeometryMed": ("The file .med for initial geometry and fields etc.", ""),
    "equation": ("Right click to choose the equation you want to apply on the .med", ""),
    "Analysis": ("Solverlab launch case parameters", ""),
  }

  def __init__(self):
    super(ModelSvl, self).__init__()
    self.setIsCast(True)
    self._defautNameAsRoot = "Solverlab"
    self._setAllAttributesList()

  def setDefaultValues(self):
    self.Analysis.setDefaultValues()

  def getMedfile(self):
    res = self.GeometryMed.fileMed
    return os.path.realpath(os.path.expandvars(res))

  def getListEquation(self):
    return self.Model

  def getBackground(self):
    return self.Analysis.caseSolverlab.launchInBackground

  def getSelectedEquation(self):
    try:
      i = int(self.Analysis.caseSolverlab.Equation)
      logger.info("Equation %s selected" % i)
      return self.Model[i]
    except:
      logger.info("No equation")
      return None

  def getEtudeWorkdirExpanded(self):
    res = self.Analysis.dataInformations.getEtudeWorkdir()
    res = os.path.realpath(os.path.expandvars(res))
    if "$" in res:
      logger.error("impossible to expand etude directory: '%s'" % res)
      return None
    return res

  def getEtudeWorkdirBrut(self, expanded=True):
    return self.Analysis.dataInformations.getEtudeWorkdirBrut()

  def getHistoryFile(self):
    return self.Analysis.historyFileManager

  def userExpand(self):
    """
        filename patterns '*,??'
        warning not for '[' ']' as 'alist[*]'
        no found way to quote meta-character
        https://docs.python.org/2/library/fnmatch.html
        """
    expanded = r"""
    Iradina
    Analysis
    dataInformations
    dataManager
    macroManager
    functions
    libraries
    macros
    userFileManager
    """.split()
    expanded = []  # TODO use expanded above
    logger.debug("userExpand\n%s" % PP.pformat(expanded))
    self.getController().UserExpandSignal.emit(expanded)


CLFX.appendAllXyzClasses([MedSvl, ListOfMedSvl, ModelSvl, ])
