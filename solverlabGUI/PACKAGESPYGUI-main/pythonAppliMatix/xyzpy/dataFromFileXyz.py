#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import traceback
import pprint as PP

try:
  import MEDLoader as ML
except:
  print("WARNING : no MEDLoader library dataFromFileXyz.py")
  ML = None



from xyzpy.baseXyz import _XyzConstrainBase, BaseFreeXyz, ListOfBaseXyz
from xyzpy.intFloatListXyz import StrNoEditionXyz, StrXyz, StrNoEditionXyz, FileXyz, BoolNoEditionXyz

import xyzpy.classFactoryXyz as CLFX
from xyzpy.guiXyz.dialogXmlXyz import DialogXyz, QTextEditXyz
import xyzpy.loggingXyz as LOG

from PyQt5 import QtGui, QtWidgets

import xyzpy.stringIO as IOX

logger = LOG.getLogger()
verbose = False
debug = False

if debug: print("TODO: dataFromFileXyz.py: set debug to False")

_example_DataArray = {
  'first name': ['Jason', 'Molly', 'Tina', 'Jake', 'Amy'],
  'last_name': ['Miller', 'Jacobson', "?", 'Milner', 'Cooze'],
  'age': [42, 52, 36, 24, 73],
  'preTestScore': [4, 24, 31, "<", ">"],
  'postTestScore': [25, 94, 57, 62, 71.5]
  }

try:
  import pandas as pd
  _example_DataFrame = pd.DataFrame(_example_DataArray, columns=sorted(_example_DataArray.keys()))
  pandasOk = True
except:
  pandasOk = False
  print("WARNING: no pandas library dataFromFileXyz.py")
  _example_DataFrame = None

###############################################################
# utilities
###############################################################
_loadedFiles = {}

def getMEDFileData(afile):
  """
  load an returns contents file med (supposedly constant)
  as MEDLoader.MEDFileData returns
  load file(s) only one time
  """
  if ML is None:
    return None

  if afile in _loadedFiles:
    return _loadedFiles[afile]
  else:
    c = ML.MEDFileData(afile)
    # print("MEDFileData contents\n", str(c))
    # help(c)
    # print("MEDFileData API %s\" % (type(c))))
    _loadedFiles[afile] = c
  return c

def getResumeMEDFileData(medData):
  """
  returns BaseFree as resume of MEDFileData contents
  """
  res = BaseFreeXyz()
  if ML is None:
    res.problem = "can't import MEDLoader"
    return None

  # getHeader getMeshes getParams
  # serialize -> DataArrayByte
  header = medData.getHeader()
  res.Header = StrNoEditionXyz(" ".join(str(header).split()))
  meshes = medData.getMeshes()
  if len(meshes) == 0:
    res.meshes = StrNoEditionXyz("No meshes")
  else:
    res.meshes = ListOfBaseXyz()
    for i, m in enumerate(meshes):
      # getName
      # getFamilies getFamiliesNames getFamilyInfo
      # getGroups getGroupInfo getGroupsNames
      # getAllGeoTypes
      # getDescription
      # getIteration
      # getJoints
      # getSpaceDimension getMeshDimension
      attr = "mesh_%i" % i
      # if i == 0: help(m)
      spacedim = m.getSpaceDimension()
      meshdim = m.getMeshDimension()
      groupsNames = m.getGroupsNames()
      groupsNames = ", ".join(groupsNames)
      # strm =  str(m).split("(* FAMILIES OF")[0]

      resm = BaseFreeXyz()
      resm.name = StrNoEditionXyz(m.getName())
      resm.spaceDim = StrNoEditionXyz(spacedim)
      resm.meshDim = StrNoEditionXyz(meshdim)
      resm.groupsNames = StrNoEditionXyz(groupsNames)
      # resm.meshInfo = StrNoEditionXyz(strm)
      res.meshes.append(resm)

  fields = medData.getFields()
  if len(fields) == 0:
    res.fields = StrNoEditionXyz("No fields")
  else:
    res.fields = ListOfBaseXyz()
    for i, f in enumerate(fields):  # MEDFileFieldMultiTS
      # field(self, iteration: 'int', order: 'int', mesh: 'MEDFileMesh') -> 'MEDCoupling::MEDCouplingFieldDouble *'
      # field(MEDFileFieldMultiTS self, int iteration, int order, MEDFileMesh mesh) -> MEDCouplingFieldDouble
      # getFieldAtLevel(self, type: 'MEDCoupling::TypeOfField', iteration: 'int', order: 'int', meshDimRelToMax: 'int', renumPol: 'int' = 0) -> 'MEDCoupling::MEDCouplingFieldDouble *'
      # getFieldAtLevel(MEDFileFieldMultiTS self, MEDCoupling::TypeOfField type, int iteration, int order, int meshDimRelToMax, int renumPol=0) -> MEDCouplingFieldDouble
      # getFieldAtTopLevel(self, type: 'MEDCoupling::TypeOfField', iteration: 'int', order: 'int', renumPol: 'int' = 0) -> 'MEDCoupling::MEDCouplingFieldDouble *'
      # getFieldAtTopLevel(MEDFileFieldMultiTS self, MEDCoupling::TypeOfField type, int iteration, int order, int renumPol=0) -> MEDCouplingFieldDouble
      # getFieldOnMeshAtLevel(self, *args) -> 'MEDCoupling::MEDCouplingFieldDouble *'
      # getFieldOnMeshAtLevel(MEDFileFieldMultiTS self, MEDCoupling::TypeOfField type, int iteration, int order, int meshDimRelToMax, MEDFileMesh mesh, int renumPol=0) -> MEDCouplingFieldDouble
      # getFieldOnMeshAtLevel(MEDFileFieldMultiTS self, MEDCoupling::TypeOfField type, int iteration, int order, MEDCouplingMesh mesh, int renumPol=0) -> MEDCouplingFieldDouble
      # LoadSpecificEntities(fileName: 'std::string const &', fieldName: 'std::string const &', entities: 'PyObject *', loadAll: 'bool' = True) -> 'MEDCoupling::MEDFileFieldMultiTS *'
      # LoadSpecificEntities(std::string const & fileName, std::string const & fieldName, PyObject * entities, bool loadAll=True) -> MEDFileFieldMultiTS
      # getInfo(self)
      # getIterations(self) -> 'PyObject *'
      # getIterations(MEDFileAnyTypeFieldMultiTS self) -> PyObject *
      # getMeshName(self) -> 'std::string'
      # getMeshName(MEDFileAnyTypeFieldMultiTS self) -> std::string
      # getName(self) -> 'std::string'
      # getName(MEDFileAnyTypeFieldMultiTS self) -> std::string
      # getNumberOfComponents(self) -> 'int'
      # getNumberOfComponents(MEDFileAnyTypeFieldMultiTS self) -> int
      # loadArrays(self) -> 'void'
      # loadArrays(MEDFileAnyTypeFieldMultiTS self)
      # resetContent(self) -> 'void'
      # resetContent(MEDFileFieldGlobsReal self)
      # serialize(self) -> 'MEDCoupling::DataArrayByte *'
      # serialize(MEDFileWritableStandAlone self) -> DataArrayByte

      # if i == 0: help(f)
      resf = BaseFreeXyz()
      resf.fieldName = StrNoEditionXyz(f.getName())
      resf.meshName = StrNoEditionXyz(f.getMeshName())
      resf.numberOfComponents = StrNoEditionXyz(f.getNumberOfComponents())
      resf.info = StrNoEditionXyz(", ".join(f.getInfo()))
      res.fields.append(resf)

  return res


###############################################################
class FileDataXyz(FileXyz):
  """
  edition file only
  """
  _envvars = "HOME WORDIRDEFAULT WORDIR4CASSIS".split()
  _title = 'Select file'
  _filter = '(*.csv;*.*)'
  
  def getActionsContextMenu(self):
    actions = []
    actions.append( self._createAction("Browse file", None, "Browse", self.browseDialog, "browsefile") )
    actions.append( self._createAction("Quick asci edit", None, "Edit in central widget", self.quickEdit, "editor") )
    return actions


###############################################################
class DataArrayXyz(StrNoEditionXyz):
  """
  designed for Xyz model leaf of pandas.dataframe (and numpy.array?)
  """
  _default_data = {}
  _csv_sep = " "
  _debug = False
  
  def __init__(self, value = None):
    super(DataArrayXyz, self).__init__()
    if value == None:
      if pandasOk:
        self._Data = pd.DataFrame(self._default_data, columns=sorted(self._default_data.keys())) #is not immutable
      else:
        print("ERROR: no pandas library for DataArrayXyz, fix it.")
        self._Data = None
    else:
      self.setData(value)
    if self._debug: print("DataArrayXyz.__init__:\n%s" % self._Data)
      
  def setData(self, data):
    if self._debug: print("DataArrayXyz.setData %s\n%s" % (type(data), data))
    if type(data) == dict:
      self._Data = pd.DataFrame(data)
    elif type(data) == str:
      self._Data = self.fromStrCsv(data)
    else:
      self._Data = data.copy()

  def getData(self):
    return self._Data.copy()
  
  def __repr__(self):
    #__repr__ goal is to be unambiguous
    cl = self.__class__
    return "%s('\n%s')" % (cl.__name__, self._Data.to_string()) #supposed long string
    
  def __str__(self):
    cl = self.__class__
    return "%s(%s)" % (cl.__name__, self._Data.columns.values.tolist()) #avoid self._Data is long string array
    
  def toXml(self, **kwargs):
    # print("DataArrayXyz kwargs:\n%s" % PP.pformat(kwargs))
    res = super(DataArrayXyz, self).toXml(**kwargs)
    res.text = self.toStrCsv() #self._Data.to_string()
    return res

  def toStrCsv(self, sep=None):
    if sep == None:
      theSep = self._csv_sep
    else:
      theSep = sep
    outStrIO = IOX.StringIO()
    # big format used for save float as xml asci....
    # print("DataArrayXyz _Data:\n%s" % PP.pformat(self._Data))
    self._Data.to_csv(outStrIO, sep=theSep, float_format="%.14e")
    res = outStrIO.getvalue()
    if self._debug: print("DataArrayXyz.toStrCsv\n%s" % res)
    return res
    
  def fromStrCsv(self, aStr, sep=None):
    if sep == None:
      theSep = self._csv_sep
    else:
      theSep = sep

    if theSep == " ":
      inStrIO = IOX.StringIO(aStr)
      #index_col=0 #no create Unnamed: 0
      #delim_whitespace=True #fix multiple whitespace
      res = pd.read_csv(inStrIO, index_col=0, delim_whitespace=True)
      if self._debug: print("DataArrayXyz.fromStrCsv sep whitespace:\n",aStr,"->\n",res)
    else:
      #strip whitespace (in case of)
      aStr2 = aStr.strip()
      while ' ' in aStr2: 
        aStr2 = aStr2.replace('  ', '')
      inStrIO = IOX.StringIO(aStr2)
      #index_col=0 #no create Unnamed: 0
      res = pd.read_csv(inStrIO, index_col=0, sep=theSep)
      if self._debug: print("DataArrayXyz.fromStrCsv sep '%s':" % theSep,aStr,"->",res)
    return res

  def getDataAsDataFrame(self):
    #res = self.getData(str(self))
    res = self.getData()
    if self._debug: print("DataArrayXyz.getDataAsDataFrame:",res)
    return res

###############################################################
class DataFromFileXyz(_XyzConstrainBase): #BaseFreeXyz):
  """
  class to evaluate data from file
  designed for data user defined type pandas.dataframe or numpy.array
  X__InitialData__ and X__CurrentData__ are implicitly present in save/restore xml
  but are hidden because big str
  """

  _debug = False
  _show_X__CurrentData__ = False #("wambeke" in os.getenv("USERNAME"))
  _icon = "datafromfile"
  _attributesList = [ #list, not a dict because sequential order list is used in files Xyz
    ("SourceFile", "FileXyz"),
    ("DataModified", "BoolNoEditionXyz"),
    ("X__InitialData__", "DataArrayXyz"),
    ("X__CurrentData__", "DataArrayXyz"),
  ]
  _icon = "TODO"
  _helpDict = {
    "SourceFile": ("Data from file csv or else...", ""),
    "DataModified": ("Current data from file is modified", ""),
  }

  def __init__(self, value=None):
    super(DataFromFileXyz, self).__init__()
    if verbose: print("DataFromFileXyz.__init__")
    self.setIsCast(True)
    self._defautNameAsRoot = "DataFromFile"
    self._setAllAttributesList()
    #self.setIsSet(True)
    #print "DataFromFileXyz type value:",type(value)
    if type(value) == str:
      self.SourceFile = value
      return

    if pandasOk:
      if type(value) == type(pd.DataFrame()):
        logger.debug("DataFromFileXyz type value not filename, suppose array for dataframe")
        self.setData(value)
        return

    if value is None: return

    logger.warning("DataFromFileXyz type value unknown: %s" % type(value))
    return


  '''
  def __setattr__(self, name, value):
    super(DataFromFileXyz, self).__setattr__(name, value)
  '''

  def setDefaultValues(self):
    self.SourceFile = ""
    self.DataModified = False
    self.X__InitialData__ = None
    self.X__CurrentData__ = None
    
  def setData(self, data=None):
    """reset both Initial and Current Data"""
    self.setInitialData(data)
    self.setCurrentData(data)
    self.setIsSet(True)
    
  def setInitialData(self, data=None):
    if verbose: print("DataFromFileXyz.setInitialData", data)
    self.X__InitialData__.setData(data)
    return
    
  def setCurrentData(self, data=None):
    if verbose: print("DataFromFileXyz.setCurrentData", data)
    self.X__CurrentData__.setData(data)
    return
    
  def setCurrentDataFromFile(self, aFile=None):
    if verbose: print("DataFromFileXyz.setCurrentDataFromFile", data)
    aStrFile = str(aFile)
    self.SourceFile = aStrFile
    if not os.path.isfile(aStrFile):
      logger.error("not a file: '%s'" % aStrFile)
      return
    #TODO filter dataframe from file extension JSON, XML, or other...
    #for now accept csv separator " " as string
    with open(aStrFile, 'r') as f: data = f.read()
    self.setCurrentData(data)
    return
    
  def clearInitialData(self):
    self.setInitialData(None)
    
  def clearCurrentData(self):
    self.setCurrentData(None)
    
  def resetCurrentData(self):
    self.setcurrentData(self.getInitialData())
    
  def getInitialData(self):
    return self.X__InitialData__.getData()
    
  def getCurrentData(self):
    return self.X__CurrentData__.getData()

  def getActionsContextMenu(self):
    """no browse"""
    if verbose: print("%s.getContextMenu" % self._className)
    actions = []
    actions.append( self._createAction('Browse file', None, 'Load data from file (csv or else...)', self.browseDialog, 'browsefile') )
    actions.append( self._createAction('Visualize/Edit data', None, 'Visualize/Edit Data', self.createEditorData, 'editor') )
    return actions

  def browseDialog(self, parent=None):
    self.SourceFile.browseDialog()
    aStrFile = str(self.SourceFile)
    if aStrFile =="": 
      return
    if not os.path.isfile(aStrFile):
      logger.error("not a file: '%s'" % aStrFile)
      return
    #TODO filter dataframe from file extension JSON, XML, or other...
    #for now accept csv separator " " as string
    with open(aStrFile, 'r') as f: data = f.read()
    
    controller = self.getController()
    if controller == None:
      self.setData(data)
    else:
      controller.setModelItemValueSignalList.emit( ["%s.setData(args[1])" % self.getTreePyName(), data ] )
    return
  
  def createEditorData(self, parent=None):
    if verbose: print("DataFromFileXyz.createEditorData parent", parent)
    import pandaspy.pandasMainWidgetXyz as PDMW
    controller = self.getController()
    aDialog = PDMW.PandasMainDialogXyz("Pandas Viewer")
    aDialog.setValue(self.X__CurrentData__._Data) #setValue do copy...
    aDialog.setMinimumSize(500, 350)
    aDialog.setWindowTitle(self.getTreePyName())

    aDialog.exec_()
    if aDialog.choice == "Cancel": return True

    newValue = aDialog.getValue()
    if verbose: 
      print("DataFromFileXyz.createEditorData: new value:'\n%s'" % newValue)
      
    if controller == None:
      logger.error("DataFromFileXyz.createEditorData: no controller for action View/Edit")
    else:
      controller.setModelItemValueSignalList.emit( ["%s.setCurrentData(args[1])" % self.getTreePyName(), newValue ] )
      controller.setModelItemValueSignal.emit( "%s.DataModified = True" % self.getTreePyName() )
    return True

  def createEditor(self, parent):
    #logger.warning("DataFromFileXyz.createEditor : NO createEditor")
    return None

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    if nameAttr == "X__InitialData__": return True
    if nameAttr == "X__CurrentData__": 
      if self._show_X__CurrentData__: 
        return False
      else:
        return True
    return False


###############################################################
class DataFromFileMedXyz(_XyzConstrainBase):  # BaseFreeXyz):
  """
  class to evaluate data from med file (as contents branch)
  """
  _attributesList = [  # list, not a dict because sequential order list is used in files Xyz
    ("fileMed", "FileOnlyBrowseXyz"),
    ("contents", "BaseFreeXyz"),
  ]
  _helpDict = {
    "fileMed": ("Med file path", ""),
    "contents": ("Current file resume data", ""),
  }
  _icon = "datafromfilemed"

  def __init__(self, value=None):
    super(DataFromFileMedXyz, self).__init__()
    self._defautNameAsRoot = "DataFromFileMed"
    self.setIsCast(True)
    self._setAllAttributesList()

  def setDefaultValues(self):
    super(DataFromFileMedXyz, self).setDefaultValues()
    self.fileMed = ""

  def __setattr__(self, name, value):
    super(DataFromFileMedXyz, self).__setattr__(name, value)
    if name == "fileMed":
      self.on_attributesChange(False)
    return

  def on_attributesChange(self, verbose=False):
    if verbose: print("%s.on_attributesChange" % self._className)
    self.checkValues(verbose)
    controller = self.getRoot().getController()
    #aStr = self.toStrXml()
    if controller != None:
      if verbose:
        print("refresh views signal of controller %s views %s" % \
              ( controller.objectName(), str([str(v.objectName()) for v in controller.getViews()]) ))
      controller.refreshViewsSignal.emit()
    else:
      if verbose: print("%s new values: controller None, no refresh views" % self._className) #,"\n",aStr

  def checkValues(self, verbose=True):

    def clear_infos(self):
      self.contents = BaseFreeXyz()
      return

    # self.name = ""
    afile = self.fileMed.getNameExpanded()

    if afile == "":
      clear_infos(self)
      return # empty value as no data yet...
    try:
      # if non reaching nfs path (or else) ET.parse(fileElement) have infinite loop
      ok = os.path.isfile(afile)  # test non infinite loop
    except:
      mess = "problem file '%s'" % afile
      if mess not in _messDone:
        logger.error(mess)
        _messDone.append(mess)
      return

    try:
      if not os.path.isfile(afile):
        logger.error("inexisting file '%s'" % afile)
        # self.name = "inexisting file '%s'" % afile
        clear_infos(self)
        return

      # read file
      contents = self.fromFileMed(afile)
      self.contents = contents
      return

    except:
      if verbose:
        traceback.print_exc() # better explicit verbose problem
        logger.error("can't read file med '%s'" % afile)
        clear_infos(self)
      return

  def fromFileMed(self, afile):
    """restore from file med with salome MEDLoader api"""
    c = getMEDFileData(afile) # get
    contents = getResumeMEDFileData(c) # to basefree
    return contents

  def getActionsContextMenu(self):
    actions = self.fileMed.getActionsContextMenu()
    return actions

  def createEditor(self, parent):
    return None

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    return False


#factory pattern using xyzpy.utilsXyz._dictOfXyzClass
CLFX.appendAllXyzClasses( [
  FileDataXyz, DataArrayXyz, DataFromFileXyz,
  DataFromFileMedXyz ] )


