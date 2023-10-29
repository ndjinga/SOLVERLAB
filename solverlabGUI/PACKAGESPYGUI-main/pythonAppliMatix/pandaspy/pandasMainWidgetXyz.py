#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


"""
Visualizing Data in Qt applications
There is experimental support for visualizing DataFrames in PyQt5 and PySide applications.
At the moment you can display and edit the values of the cells in the DataFrame.
Qt will take care of displaying just the portion of the DataFrame that is currently visible
and the edits will be immediately saved to the underlying DataFrame
To demonstrate this we will create a simple PySide application 
that will switch between two editable DataFrames.
For this will use the DataFrameModel class that handles 
the access to the DataFrame, and the DataFrameWidget,
which is just a thin layer around the QTableView.
http://pandas.pydata.org/pandas-docs/stable/faq.html
http://pandas.pydata.org/pandas-docs/stable/io.html#io-read-csv-table
https://matplotlib.org/users/mathtext.html#mathtext-tutorial
"""

import os, sys
import traceback
from collections import OrderedDict

import pprint as PP
import numpy as np
import pandas as pd
from PyQt5 import QtGui, QtCore, QtWidgets

#import xyzpy.actionsFactoryXyz as ACFX
import xyzpy.utilsXyz as UXYZ
import salomepy.iconsUser as IUSR
from matplotlibpy.matplotlibWindowToolbar import MatplotlibWindowToolbar
from xyzpy.guiXyz.dialogXmlXyz import DialogXyz
import widgetpy.messageBoxDialog as MBD
import iradinapy.dateTime as DATT

verbose = False
debug = False

"""
# example from dict: sorted columns keys
_example_DataArray = {
  'first name': ['Jason', 'Molly', 'Tina', 'Jake', 'Amy'],
  'last_name': ['Miller', 'Jacobson', "?", 'Milner', 'Cooze'],
  'age': [42, 52, 36, 24, 73],
  'preTestScore': [4, 24, 31, "<", ">"],
  'postTestScore': [25, 94, 57, 62, 71.5]
  }
 
_example_DataFrame = pd.DataFrame(
     _example_DataArray, 
     columns=sorted(_example_DataArray.keys()))
"""

# example from list: columns keys as it
_example_DataArray = OrderedDict([
  ('first name', ['Jason', 'Molly', 'Tina', 'Jake', 'Amy']),
  ('last_name', ['Miller', 'Jacobson', "?", 'Milner', 'Cooze']),
  ('age', [42, 52, 36, 24, 73]),
  ('preTestScore', [4, 24, 31, "<", ">"]),
  ('postTestScore', [25, 94, 57, 62, 71.5]),
  ])
 
# FutureWarning: from_items is deprecated
_example_DataFrame = pd.DataFrame.from_dict(_example_DataArray)

'''example as expected file .csv, separator " ":
 age "first name" last_name postTestScore preTestScore
0 42 Jason Miller 25.0 4
1 52 Molly Jacobson 94.0 24
2 36 Tina ? 57.0 31
3 24 Jake Milner 62.0 <
4 73 Amy Cooze 71.5 >
'''

########################################################
def mergeDicts(dict1, dict2):
  """
  compatible python2-3 res={**dict1, **dict2}
  overriding duplicates keys from dict2
  """
  res = {}
  for k, v in dict1.items():
    res[k] = v
  for k, v in dict2.items():
    res[k] = v
  return res

########################################################
def getMode():
  """
  get a value for mode as iradinaGUI 'simple' 'advanced' etc.
  as simple menus, or not
  """
  try:
    import iradinapy.configIra as CFGIRA
    mode = CFGIRA.getCurrentMode()
  except:
    mode = "advanced"
    pass
  if verbose:
    print("pandas main widget user mode %s" % mode)
  return mode


########################################################
def getExampleDataFrameEllipse(periods = 1000):
  df = pd.DataFrame({"teta": np.linspace(0, np.pi*2*5, periods)} )
  # vectorized
  r = np.linspace(1, 0.1, periods)
  df["X"] = np.cos(df.teta)*r
  df["Y"] = np.sin(df.teta)*r
  aDict = getDictInfoFromDataFrame(df)
  infos = {"suptitle": "from getExampleDataFrameEllipse", 
           "title": "ellipse", 
           "xtitle": "x", 
           "ytitle": "y=f(x)", 
           "periods": periods, 
           }
  aDict.update(infos)
  df_infos = pd.DataFrame(aDict, index=[0]).transpose()
  """
  aPlot = df.plot(x="X", y="Y", label="ellipse", legend=False)
  lim = 1.5
  aPlot.set_xlim(-lim,lim)
  aPlot.set_ylim(-lim,lim)
  """
  return (df_infos, df)


########################################################
def getExampleDataFrame_1():
  df2 = pd.DataFrame(_example_DataArray, columns=sorted(_example_DataArray.keys()))
  df1 = getInfoFromDataFrame(df2)
  return (df1, df2) 

########################################################
def read_csv_type(aNameFile, typeCsv="Xyz"):
  """
  read a csv file from different expected formats
  return 2 pandas DataFrame(s):
    infos: from header, 2 columns
    contents: 2d array of data
  known types are 'Xyz' 'FirstCol' 'IradinaPlot'
  """
  if typeCsv == "Xyz": return read_csv_type_Xyz(aNameFile)
  if typeCsv == "FirstCol": return read_csv_type_FirstCol(aNameFile)
  if typeCsv == "IradinaPlot":
    try:
      if UXYZ.getFirstCharacterInFile(aNameFile) == "#": # header iradina > 1.2.0
        return read_csv_type_IradinaPlot(aNameFile)
      else:  # no header
        return read_csv_type_IradinaPlot_without_header(aNameFile)
    except Exception as e:
      wid = MBD.getMessageBoxDialog(exception=e)
      wid.exec_()
      raise Exception("problem reading file csv '%s'" % aNameFile)
  raise Exception("type file csv '%s' unknown" % typeCsv)


########################################################
def read_csv_type_Xyz(aNameFile, sep=" "):
  """
  read a csv file from Xyz package
  return 2 pandas DataFrame(s):
    infos: from header, 2 columns
    contents: 2d array of data
  separator = " " for Xyz
  """
  # contents
  df2 = pd.read_csv(aNameFile, index_col=0, sep=sep)
  # infos
  aDict = getDictInfoFromDataFrame(df2)
  aDict["initialFile"] = aNameFile
  df1 = pd.DataFrame(aDict, index=[0]).transpose()
  return (df1, df2)


########################################################
def read_csv_type_FirstCol(aNameFile, sep=" "):
  """
  read a csv file from Xyz package
  return 2 pandas DataFrame(s):
    infos: from header, 2 columns
    contents: 2d array of data
  separator = " " for Xyz
  """
  # contents
  df2 = pd.read_csv(aNameFile, index_col=False, sep=sep)
  # infos
  aDict = getDictInfoFromDataFrame(df2)
  aDict["initialFile"] = aNameFile
  df1 = pd.DataFrame(aDict, index=[0]).transpose()
  return (df1, df2)

########################################################
def read_csv_type_IradinaPlot_3_columns_before_iradina120(aNameFile, sep="\t", logger=None):
  """
  read a iradina output file
  return 2 pandas DataFrame(s):
    infos: from header, 4 columns (indice axe x, indice axe y, indice axe  z, value)
    contents: 2d array of data
  separator = tabulation (from iradina code)
  """
  # contents
  # title of fourth column from filename
  aDir, name = os.path.split(aNameFile)
  tit = name.split(".")[1:] # remove "ira." as first postfix of namefile (ex: "ira.ions.replacements")
  tit = ".".join(tit)
  namesCols = "ix iy iz %s" % tit
  df2 = pd.read_csv(aNameFile, sep=sep, header=None, names=namesCols.split())
  if verbose: print("read_csv_type_IradinaPlot\n%s" % df2)
  # infos
  aDict = getDictInfoFromDataFrame(df2)
  aDict["initialFile"] = aNameFile
  df1 = pd.DataFrame(aDict, index=[0]).transpose()

  # search for info(s) in expected existing output/../Structure.in, if existing
  # Structure file for iradina (example)
  """
  ...
  [Target]
  cell_count_x=250
  cell_count_y=1
  cell_count_z=1
  cell_size_x=1.0
  cell_size_y=1000.0
  cell_size_z=1000.0
  periodic_boundary_x=0
  periodic_boundary_y=1
  periodic_boundary_z=1
  CompositionFileType=0
  CompositionFileName=./Composition.in
  UseDensityMultiplicator=0
  DensityMultiplicatorFileName=none
  special_geometry=0
  ...
  """
  structFile = os.path.realpath(os.path.join(aDir, "..", "Structure.in"))
  if verbose: print("structFile %s", structFile)
  if os.path.exists(structFile):
    import configparserpy.configParserUtils as CPAU
    cfg = CPAU.getConfigFromFile(structFile)
    # cfg = cfg.toCatchAll()
    # print("target file structure contents:\n%s" % cfg.Target)
    cfg = cfg.toDict()
    target = cfg["Target"]
    if verbose: print("target file structure contents:\n%s" % PP.pformat(target))
    aDict = mergeDicts(aDict, target) # merge overriding aDict existing duplicate keys
    df1 = pd.DataFrame(aDict, index=[0]).transpose()
    for i in "x y z".split(): # create x y z columns in (nm)
      ini = "i" + i    # ix iy iz
      end = i + "(nm)" # x(nm) y(nm) z(nm)
      size_i = float(target["cell_size_" + i])
      df2[end] = df2[ini] * size_i
  else:
    raise Exception("inexisting target file structure '%s'" % structFile)
    if logger is None:
      print("WARNING: %s" % msg)
    else:
      logger.warning(msg)
  return (df1, df2)

########################################################
def read_csv_type_IradinaPlot_without_header(aNameFile, sep="\t", logger=None):
  """
  read a iradina output file without heade
  return 2 pandas DataFrame(s):
    infos: from header, number of columns unknown
    contents: 2d array of data
  separator = tabulation (from iradina code)
  """
  # contents
  # title of fourth column from filename
  aDir, name = os.path.split(aNameFile)
  tit = name.split(".")[1:] # remove "ira." as first postfix of namefile (ex: "ira.ions.replacements")
  tit = ".".join(tit)

  # namesCols = "ix iy iz %s" % tit
  df2 = pd.read_csv(aNameFile, sep=sep, header=None)
  # infos
  aDict = getDictInfoFromDataFrame(df2)
  aDict["initialFile"] = aNameFile
  df1 = pd.DataFrame(aDict, index=[0]).transpose()
  nbCols = df2.shape[1]
  names = "ix iy iz".split()
  # print("nbCols", nbCols)
  if nbCols > 3:
    names = "ix iy iz".split()
    renames = {}
    for i in range(nbCols):
      if i < 3:
        renames[i] = names[i]
      else:
        renames[i] = "C%i" %i
    #print("renames", renames)
    df2 = df2.rename(renames, axis=1)  # new method
  # print("read_csv_type_IradinaPlot\n%s" % df2)
  # search for info(s) in expected existing output/../Structure.in, if existing
  # Structure file for iradina (example)
  """
  ...
  [Target]
  cell_count_x=250
  cell_count_y=1
  cell_count_z=1
  cell_size_x=1.0
  cell_size_y=1000.0
  cell_size_z=1000.0
  periodic_boundary_x=0
  periodic_boundary_y=1
  periodic_boundary_z=1
  CompositionFileType=0
  CompositionFileName=./Composition.in
  UseDensityMultiplicator=0
  DensityMultiplicatorFileName=none
  special_geometry=0
  ...
  """
  structFile = os.path.realpath(os.path.join(aDir, "..", "Structure.in"))
  if verbose: print("structFile %s", structFile)
  if os.path.exists(structFile):
    import configparserpy.configParserUtils as CPAU
    cfg = CPAU.getConfigFromFile(structFile)
    # cfg = cfg.toCatchAll()
    # print("target file structure contents:\n%s" % cfg.Target)
    cfg = cfg.toDict()
    target = cfg["Target"]
    if verbose: print("target file structure contents:\n%s" % PP.pformat(target))
    aDict = mergeDicts(aDict, target) # merge overriding aDict existing duplicate keys
    df1 = pd.DataFrame(aDict, index=[0]).transpose()

    for i in "x y z".split(): # create x y z columns in (nm)
      ini = "i" + i    # ix iy iz
      end = i + "(nm)" # x(nm) y(nm) z(nm)
      size_i = float(target["cell_size_" + i])
      df2[end] = df2[ini] * size_i
  else:
    raise Exception("inexisting target file structure '%s'" % structFile)
    if logger is None:
      print("WARNING: %s" % msg)
    else:
      logger.warning(msg)
  return (df1, df2)

########################################################
def read_header_csv(aNameFile):
  """
  if ok, file close is not done, supposedly caller continue read next array...
  """
  f = open(aNameFile, "r")
  header = {}
  for i in range(100): # precaution...
    try:
      line = f.readline()[0:-1] # remove trailing RC
    except: # may be EOF
      print("ERROR read_header_csv file '%s'" % aNameFile)
      f.close()
      return (header, None)
      pass

    # print("read_header_csv line '%s'" % line)
    line = line.strip()
    if line == "":
      return (header, f) # end header caller continue read next array
    else:
      add_tag_csv(header, line)
  print("ERROR read_header_csv file header too long '%s'" % aNameFile)
  f.close()
  return (header, None)

########################################################
def add_tag_csv(header, line):
  if line[0] != "#":
    raise Exception("add_tag_csv needs '#' as first character in line\n'%s'" % line)

  try:
    idx = line.index(':')
  except:
    raise Exception("add_tag_csv needs ':' in '%s'" % line)

  tag = line[1:idx].strip()
  value = line[idx+1:].strip()
  if verbose: print("add_tag_csv found %s" % PP.pformat([tag, value]))
  try:
    idx = line.index('|')
  except:
    idx = None

  if idx is not None: # value is list of COLUMN_NAMES (for example)
    value = [v.strip() for v in value.split('|')]
  header[tag] = value
  return

########################################################
def read_csv_type_IradinaPlot(aNameFile, sep="\t", logger=None):
  """
  read a iradina output file (after iradina code version 1.2.0)
  return 2 pandas DataFrame(s):
    infos: from header, 4 columns (indice axe x, indice axe y, indice axe  z, value)
    contents: 2d array of data
  separator = tabulation (from iradina code)
  """
  # contents
  # title of fourth column from filename
  aDir, name = os.path.split(aNameFile)
  tit = name.split(".")[1:] # remove "ira." as first postfix of namefile (ex: "ira.ions.replacements")
  tit = ".".join(tit)
  header, fileStream = read_header_csv(aNameFile)
  if verbose: print("read_csv_type_IradinaPlot header:\n%s" % PP.pformat(header))
  # namesCols = "ix iy iz %s" % tit
  namesCols = header['COLUMN_TITLES']
  df2 = pd.read_csv(fileStream,  sep=sep, header=None, names=namesCols)
  if verbose: print("read_csv_type_IradinaPlot\n%s" % df2)
  # infos
  aDict = getDictInfoFromDataFrame(df2)
  aDict["initialFile"] = aNameFile
  df1 = pd.DataFrame(aDict, index=[0]).transpose()
  aDict = mergeDicts(aDict, header) # merge overriding aDict existing duplicate keys
  if verbose: print("read_csv_type_IradinaPlot aDict:\n%s" % PP.pformat(aDict))


  # search for info(s) in expected existing output/../Structure.in, if existing
  # Structure file for iradina (example)
  """
  ...
  [Target]
  cell_count_x=250
  cell_count_y=1
  cell_count_z=1
  cell_size_x=1.0
  cell_size_y=1000.0
  cell_size_z=1000.0
  periodic_boundary_x=0
  periodic_boundary_y=1
  periodic_boundary_z=1
  CompositionFileType=0
  CompositionFileName=./Composition.in
  UseDensityMultiplicator=0
  DensityMultiplicatorFileName=none
  special_geometry=0
  ...
  """
  structFile = os.path.realpath(os.path.join(aDir, "..", "Structure.in"))
  if verbose: print("structFile %s", structFile)
  if os.path.exists(structFile):
    import configparserpy.configParserUtils as CPAU
    cfg = CPAU.getConfigFromFile(structFile)
    # cfg = cfg.toCatchAll()
    # print("target file structure contents:\n%s" % cfg.Target)
    cfg = cfg.toDict()
    target = cfg["Target"]
    if verbose:
      print("target file structure contents:\n%s" % PP.pformat(target))
    aDict = mergeDicts(aDict, target) # merge overriding aDict existing duplicate keys
    aDictString = {}
    for k, value in aDict.items():
      if type(value) is list:
        aDictString[k] = PP.pformat(value)
      else:
        aDictString[k] = value
    df1 = pd.DataFrame(aDictString, index=[0]).transpose()
    # for i in "x y z".split(): # create x y z columns in (nm)
    #   ini = "i" + i    # ix iy iz
    #   end = i + "(nm)" # x(nm) y(nm) z(nm)
    #   size_i = float(target["cell_size_" + i])
    #   df2[end] = df2[ini] * size_i
  else:
    msg = "inexisting target file structure '%s'" % structFile
    if logger is None:
      print("WARNING: %s" % msg)
    else:
      logger.warning(msg)
  return (df1, df2)


########################################################
def getInfoFromDataFrame(df):
  aDict = getDictInfoFromDataFrame(df)
  if verbose:
    print("getInfoFromDataFrame %s" % aDict)
  df1 = pd.DataFrame(aDict, index=[0]).transpose()
  return df1


########################################################
def getDictInfoFromDataFrame(df):
  aDict = {}
  nbRows, nbCols = df.shape
  aDict["nbCols"] = str(nbCols)
  aDict["nbRows"] = str(nbRows)
  return aDict


########################################################
class DataFrameModel(QtCore.QAbstractTableModel):
  """
  data model for a DataFrame class
  initial source:
  from pandas.sandbox.qtpandas import DataFrameModel, DataFrameWidget
  """
  
  def __init__(self):
    super(DataFrameModel, self).__init__()
    self._df = pd.DataFrame()

  def setDataFrame(self, dataFrame):
    self._df = dataFrame.copy()

  def signalUpdate(self):
    """
    tell viewers to update their data (this is full update, not efficient)'
    """
    self.layoutChanged.emit()

  def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
    if role != QtCore.Qt.DisplayRole:
        return QtCore.QVariant()

    if orientation == QtCore.Qt.Horizontal:
        try:
            return self._df.columns.tolist()[section]
        except (IndexError, ):
            return QtCore.QVariant()
    elif orientation == QtCore.Qt.Vertical:
        try:
            # return self.df.index.tolist()
            return self._df.index.tolist()[section]
        except (IndexError, ):
            return QtCore.QVariant()

  def data(self, index, role=QtCore.Qt.DisplayRole):
    if role != QtCore.Qt.DisplayRole:
        return QtCore.QVariant()

    if not index.isValid():
        return QtCore.QVariant()

    """
    https://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-integer
    Warning: Starting in 0.20.0, the .ix indexer is deprecated, 
             in favor of the more strict .iloc and .loc indexers. 
    Fast scalar value getting and setting
    Since indexing with [] must handle a lot of cases (single-label access, 
    slicing, boolean indexing, etc.), 
    it has a bit of overhead in order to figure out what you’re asking for. 
    If you only want to access a scalar value, 
    the fastest way is to use the at and iat methods, 
    which are implemented on all of the data structures.

    Similarly to loc, at provides label based scalar lookups, 
    while, iat provides integer based lookups analogously to iloc
    """

    #if debug: print "data DisplayRole", self._df.iat[index.row(), index.column()]
    return QtCore.QVariant(str(self._df.iat[index.row(), index.column()]))

  def flags(self, index):
    flags = super(DataFrameModel, self).flags(index)
    flags |= QtCore.Qt.ItemIsEditable
    return flags

  def setData(self, index, value, role):
    try:
      self._df.at[self._df.index[index.row()], self._df.columns[index.column()]] = value.value()
    except:
      # type value as str in user gui modifiying
      try:
        self._df.at[self._df.index[index.row()], self._df.columns[index.column()]] = value
      except Exception as e:
        print("problem pandas serie set value: %s\n%s" % (value, e))
        return False
    return True

  def rowCount(self, index=QtCore.QModelIndex()):
    return self._df.shape[0]

  def columnCount(self, index=QtCore.QModelIndex()):
    return self._df.shape[1]
    
  def getDataFrame(self):
    return self._df
  

########################################################
class DataFrameWidget(QtWidgets.QWidget):
  """
  a simple widget for using DataFrames in a gui
  from pandas.sandbox.qtpandas import DataFrameModel, DataFrameWidget
  """
  
  def __init__(self, dataFrame, parent=None):
    super(DataFrameWidget, self).__init__(parent)
    self.dataModel = DataFrameModel()
    self.dataTable = QtWidgets.QTableView()
    layout = QtWidgets.QVBoxLayout()
    layout.addWidget(self.dataTable)
    self.setLayout(layout)
    self.setDataFrame(dataFrame)
    self.dataTable.horizontalHeader().setStretchLastSection(True)
    #self.header().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
    #self.dataTable.setColumnWidth(0, 300)

  def resizeColumnsToContents(self):
    self.dataTable.resizeColumnsToContents()
      
  def setDataFrame(self, dataFrame):
    self.dataModel.setDataFrame(dataFrame)
    self.dataTable.setModel(self.dataModel)
    self.dataModel.signalUpdate()
  
  def getSelected(self):
    select = self.dataTable.selectionModel();
    #print "has selection", select.hasSelection()
    #print "rows selection", select.selectedRows()
    #print "columns selection", select.selectedColumns()
    return (select.hasSelection(), select.selectedRows(), select.selectedColumns())
    
  def getDataFrame(self):
    return self.dataModel.getDataFrame()

  def setSelectedColumns(self, columns):
    """colums as [0,2,...] for example"""
    #print("setSelectedColumns", columns)
    model = self.dataModel
    tableView = self.dataTable
    select = tableView.selectionModel()
    #print("type(select)", type(select)) # PyQt5.QtCore.QItemSelectionModel
    indexes = [model.index(0, c) for c in columns]
    #print("indexes", indexes)
    mode = select.Select | select.Columns
    select.clearSelection()
    for index in indexes:
      select.select(index, mode)
    return

########################################################
class MyQMessageBox(QtWidgets.QMessageBox):

  def resizeEvent(self, Event):
    # print("MyQMessageBox.resizeEvent")
    QtWidgets.QMessageBox.resizeEvent(self, Event);
    self.setFixedWidth(700);
    # self.setFixedHeight(200);


########################################################
class PandasTabWidget(QtWidgets.QTabWidget):
  
  TAB_INFO=0
  TAB_CONTENTS=1
  TAB_PLOTS=2
  index = [0]
  
  def __init__(self, plotView=None):
    super(PandasTabWidget, self).__init__()
    self.setObjectName("PandasTabWidget"+str(self.index))
    self.setWindowTitle(self.objectName())
    self.index[0] += 1
    
    self.infos = DataFrameWidget(pd.DataFrame(np.arange(0))) #not empty
    self.contents = DataFrameWidget(pd.DataFrame(np.arange(0))) #not empty
    self.contents.dataTable.horizontalHeader().setStretchLastSection(False)
    if plotView==None:
      self.plots = MatplotlibWindowToolbar()
    else:
      self.plots = plotView #set from a controller for example
    
    self.insertTab(self.TAB_INFO, self.infos, "Infos")
    self.insertTab(self.TAB_CONTENTS, self.contents, "Contents")
    self.insertTab(self.TAB_PLOTS, self.plots, "Plots")
    self.setCurrentIndex(self.TAB_CONTENTS)
    self.resize(600,700)
    self.messDone = False
    
  def read_csv(self, aNameFile, typeCsv="Xyz"):
    try:
      self._df1, self._df2 = read_csv_type(aNameFile, typeCsv=typeCsv)
      self.infos.setDataFrame(self._df1)
      self.contents.setDataFrame(self._df2)
      self.setWindowTitle(os.path.basename(aNameFile))
      #self.infos.resizeColumnsToContents()
      #self.contents.resizeColumnsToContents()
    except Exception as e:
      mess = "Problem reading file\n%s" % (aNameFile)
      wid = MBD.getMessageBoxDialog(self, mess, exception=e)
      ret = wid.exec_()
    return
    
  def setInfos(self, value=None):
    df1 = value
    if df1 == None: df1 = getInfoFromDataFrame(self.contents.getDataFrame())
    self.infos.setDataFrame(df1)

  def getInfosAsDict(self):
    df = self.infos.getDataFrame() # current displayed dataframe in info
    res = OrderedDict()
    for tag, value in df.iterrows():
      # print("type %s %s" % (type(tag), type(value)))
      try:
        expr = value.iloc[0]
      except:
        print("problem get tag '%s'", tag)
        expr = ""
      try:
        cmd = "res['%s'] = %s" % (tag, expr)
        # print("exec command %s" % cmd)
        exec( cmd , {"res": res} ) # python 2-3
      except:
        try:
          cmd = "res['%s'] = '%s'" % (tag, expr)
          # print("exec command %s" % cmd)
          exec( cmd , {"res": res} ) # python 2-3
        except:
          print("problem conversion tag '%s' value '%s'", tag, expr)
          print("from DataFrame %s" % df)

    # print("self.infos.getDataFrame() %s" % df)
    # print("res\n%s" % PP.pformat(res))
    return res
        
  def on_read_click(self):
    """read a csv in DataFrame"""
    aDir = os.getenv("WORKDIR")
    if aDir == None:
      if self.messDone:
        QtWidgets.QMessageBox.warning(self, "warning", "You could set $WORKDIR as your directory of files csv.")
        self.messDone = true
      aDir = os.path.split(__file__)[0]
    aNameFile = str(QtWidgets.QFileDialog.getOpenFileName(self, 'Load file .csv', aDir, "(*.csv *.CSV)")[0])
    if aNameFile == '': return
    self.read_csv(aNameFile)

  def on_example_plot1_click(self):
    """plot from 2 DataFrames columns"""
    df_infos, df_contents = getExampleDataFrameEllipse(periods=1000)
    self.infos.setDataFrame(df_infos)
    infos = df_infos.transpose()
    #print "infos", df_infos
    self.contents.setDataFrame(df_contents)
    #x, y1, suptitle, title, xlabel, ylabel
    self.plots.set_curve(df_contents.X, df_contents.Y, "", infos.title[0], infos.xtitle[0], infos.ytitle[0])
    self.setCurrentIndex(self.TAB_PLOTS)

  def on_plot1_click(self):
    """plot from 2 DataFrames columns"""
    has, rows, columns = self.contents.getSelected()
    if has == False:
      QtWidgets.QMessageBox.warning(self, "warning",
          "Nothing to plot. you have to select 2 columns in contents.")
      return
    if len(columns) != 2:
      QtWidgets.QMessageBox.warning(self, "warning", "You have to select 2 columns in contents.")
      return
    col0, col1 = columns[0].column(), columns[1].column()
    if debug: print("selected columns", col0, col1)
    contents = self.contents.getDataFrame()
    names = contents.columns.values.tolist()
    x = contents.loc[:, names[col0]]
    y = contents.loc[:, names[col1]]
    suptitle = ""
    title = names[col1]
    xtitle = names[col0]
    ytitle =  ""
    self.plots.set_curve(x, y, suptitle, title, xtitle, ytitle)
    self.setCurrentIndex(self.TAB_PLOTS)

  def check_infos_iradina(self, infos):
    """complete as posible dict infos for iradina version < 1.2.0 without header"""
    nbCols = infos["nbCols"]

    try:
      units = infos['COLUMN_UNITS']
    except:
      units = [""] * nbCols
      infos['COLUMN_UNITS'] = units

    try:
      factor = infos['COLUMN_FACTOR']
    except:
      factor = ["1"] * nbCols
      infos['COLUMN_FACTOR'] = factor

    try:
      title = infos['TITLE']
    except:
      title = infos['initialFile']
      title = os.path.basename(title)
      infos['TITLE'] = title

    filename = infos['initialFile']
    filename = os.path.basename(title)

    try:
      title = infos['TITLE']
    except:
      infos['TITLE'] = filename

    try:
      name = infos['NAME']
    except:
      infos['NAME'] = filename

    try:
      datec = infos['DATE']
    except:
      now = DATT.DateTime("now")
      infos['DATE'] = str(now)

    return infos

  def on_plot1_click_iradina(self, infos=None):
    """plot from 2 DataFrames columns"""
    if infos is None:
      cinfos = tabWid.getInfosAsDict()
    else:
      cinfos = infos
    cinfos = self.check_infos_iradina(cinfos)
    self.contents.setSelectedColumns([0, 3]) # TODO something clever
    has, rows, columns = self.contents.getSelected()
    if has == False:
      QtWidgets.QMessageBox.warning(self, "warning",
          "Nothing to plot. you have to select 2 columns in contents.")
      return
    if len(columns) != 2:
      QtWidgets.QMessageBox.warning(self, "warning", "You have to select 2 columns in contents.")
      return
    col0, col1 = columns[0].column(), columns[1].column()
    # print("iradina selected columns", col0, col1)
    contents = self.contents.getDataFrame()
    names = contents.columns.values.tolist()
    units = cinfos['COLUMN_UNITS']
    factor = [float(i) for i in infos['COLUMN_FACTOR']]
    x = contents.loc[:, names[col0]] * factor[col0]
    y = contents.loc[:, names[col1]] * factor[col1]
    suptitle = ""
    title = cinfos["TITLE"]
    comment = "%s - %s" % (os.path.basename(cinfos["NAME"]), cinfos["DATE"])
    xtitle = "%s (%s)" % (names[col0], units[col0])
    ytitle = "%s (%s)" % (names[col1], units[col1])
    self.plots.set_curve(x, y, suptitle, title, xtitle, ytitle, comment=comment)
    self.setCurrentIndex(self.TAB_PLOTS)
    self.create_plot1_iradina_csv(x, y, suptitle, title, xtitle, ytitle, cinfos)

  def get_iradinaKeysCsv(self):
    """get valid iradina keys copy for creating files _plot_csv"""
    return ['NAME', 'DATE', 'TITLE', 'initialFile']

  def create_plot1_iradina_csv(self, x, y, suptitle, title, xtitle, ytitle, infos):
    try:
      # print("create_plot1_iradina_csv\n%s", PP.pformat(infos))
      cinfos = OrderedDict()
      validKeys = self.get_iradinaKeysCsv()
      for k in validKeys:
        try:
          cinfos[k] = infos[k]
        except:
          pass # as possible, no warning in missed
      cinfos["COLUMN_NAMES"] = "x | y"
      cinfos["COLUMN_TITLES"] = "%s | %s" % (xtitle, ytitle)
      aFile = cinfos["initialFile"] + "_plot.csv"
      with open(aFile, "w") as f:
        for k, v in cinfos.items():
          f.write("#%s: %s\n" % (k, v))
        f.write("\n")
        for xx, yy in zip(x,y):
          f.write("%3g\t%1g\n" % (xx, yy))
    except Exception as e:
      wid = MBD.getMessageBoxDialog(exception=e)
      wid.exec_()
    return

  def on_plot2_click(self):
    """plot 2 curves from 3 DataFrames columns (x, y1, y2) for y1=f(x) and y2=f(x)"""
    has, rows, columns = self.contents.getSelected()
    if has == False:
      QtWidgets.QMessageBox.warning(self, "warning", "Nothing to plot. you have to select 3 columns (x,y1,y2) in contents.")
      return
    if len(columns) != 3:
      QtWidgets.QMessageBox.warning(self, "warning", "You have to select 3 columns (x,y1=f(x),y2=f(x)) in contents.")
      return
    col0, col1, col2 = columns[0].column(), columns[1].column(), columns[2].column()
    #print "columns", col0, col1
    contents = self.contents.getDataFrame()
    names = contents.columns.values.tolist()
    x = contents.loc[:, names[col0]]
    y1 = contents.loc[:, names[col1]]
    y2 = contents.loc[:, names[col2]]
    suptitle = ""
    title = ""
    titles_xy = [names[col0], names[col1], names[col2]]
    xy = [x, y1, y2]
    self.plots.set_2_curves(xy, suptitle, title, titles_xy)
    self.setCurrentIndex(self.TAB_PLOTS)

  def on_plot2_click_iradina(self, infos=None):
    """plot 2 curves from 3 DataFrames columns (x, y1, y2) for y1=f(x) and y2=f(x)"""
    if infos is None:
      cinfos = tabWid.getInfosAsDict()
    else:
      cinfos = infos
    cinfos = self.check_infos_iradina(cinfos)
    self.contents.setSelectedColumns([0, 4, 3]) # TODO something clever
    has, rows, columns = self.contents.getSelected()
    if has == False:
      QtWidgets.QMessageBox.warning(self, "warning", "Nothing to plot. you have to select 3 columns (x,y1,y2) in contents.")
      return
    if len(columns) != 3:
      QtWidgets.QMessageBox.warning(self, "warning", "You have to select 3 columns (x,y1=f(x),y2=f(x)) in contents.")
      return
    col0, col1, col2 = columns[0].column(), columns[1].column(), columns[2].column()
    # print("iradina selected columns", col0, col1, col2)
    contents = self.contents.getDataFrame()
    names = contents.columns.values.tolist()
    units = cinfos['COLUMN_UNITS']
    factor = [float(i) for i in cinfos['COLUMN_FACTOR']]
    x = contents.loc[:, names[col0]] * factor[col0]
    y1 = contents.loc[:, names[col1]] * factor[col1]
    y2 = contents.loc[:, names[col2]] * factor[col2]
    suptitle = ""
    title = cinfos["TITLE"]
    comment = "%s - %s" % (os.path.basename(cinfos["NAME"]), cinfos["DATE"])
    titles_xy = ["%s (%s)" % (names[i], units[i]) for i in [col0, col1, col2]]
    xy = [x, y1, y2]
    self.plots.set_2_curves(xy, suptitle, title, titles_xy, comment=comment)
    self.setCurrentIndex(self.TAB_PLOTS)
    self.create_plot2_iradina_csv(xy, suptitle, title, titles_xy, cinfos)

  def create_plot2_iradina_csv(self, xy, suptitle, title, titles_xy, infos):
    try:
      # print("create_plot2_iradina_csv\n%s", PP.pformat(infos))
      cinfos = OrderedDict()
      validKeys = self.get_iradinaKeysCsv()
      for k in validKeys:
        try:
          cinfos[k] = infos[k]
        except:
          pass # as possible, no warning in missed
      cinfos["COLUMN_NAMES"] = "x | y1 | y2"
      cinfos["COLUMN_TITLES"] = "%s | %s | %s" % (titles_xy[0], titles_xy[1], titles_xy[2])
      aFile = cinfos["initialFile"] + "_plot.csv"
      with open(aFile, "w") as f:
        for k, v in cinfos.items():
          f.write("#%s: %s\n" % (k, v))
        f.write("\n")
        for xx, yy1, yy2 in zip(xy[0], xy[1], xy[2]):
          f.write("%3g\t%12g\t%12g\n" % (xx, yy1, yy2))
    except Exception as e:
      wid = MBD.getMessageBoxDialog(exception=e)
      wid.exec_()
    return

  def on_plotn_click(self):
    """plot n curves from n+1 DataFrames columns (x, y1, y2,... yn) for y1=f(x) and y2=f(x)..."""
    has, rows, columns = self.contents.getSelected()
    if has == False:
      QtWidgets.QMessageBox.warning(self, "warning", "Nothing to plot. you have to select n+1 columns (x,y1,y2,yn) in contents.")
      return
    if len(columns) <= 1:
      QtWidgets.QMessageBox.warning(self, "warning", "You have to select n+1 columns (x,y1=f(x),... yn=f(x)) in contents.")
      return
    contents = self.contents.getDataFrame()
    names = contents.columns.values.tolist()
    icols = [col.column() for col in columns]
    # xy = [contents[[icol]] for icol in icols]
    xy = [contents.loc[:, names[coli]] for coli in icols]
    suptitle = ""
    title = ""
    titles_xy = [names[icol] for icol in icols]
    self.plots.set_n_curves(xy, suptitle, title, titles_xy)
    self.setCurrentIndex(self.TAB_PLOTS)

  def on_plotd_click(self):
    """plot distinct n curves from n+1 DataFrames columns (x, y1, y2,... yn) for y1=f(x) and y2=f(x)..."""
    has, rows, columns = self.contents.getSelected()
    if has == False:
      QtWidgets.QMessageBox.warning(self, "warning", "Nothing to plot. you have to select n+1 columns (x,y1,y2,yn) in contents.")
      return
    if len(columns) <= 1:
      QtWidgets.QMessageBox.warning(self, "warning", "You have to select n+1 columns (x,y1=f(x),... yn=f(x)) in contents.")
      return
    contents = self.contents.getDataFrame()
    names = contents.columns.values.tolist()
    icols = [col.column() for col in columns]
    # xy = [contents[[icol]] for icol in icols]
    xy = [contents.loc[:, names[coli]] for coli in icols]
    suptitle = ""
    title = ""
    titles_xy = [names[icol] for icol in icols]
    self.plots.set_n_plots(xy, suptitle, title, titles_xy)
    self.setCurrentIndex(self.TAB_PLOTS)

  def on_ploth_click(self):
    """plot historam 1 DataFrame columns"""
    has, rows, columns = self.contents.getSelected()
    if has == False:
      QtWidgets.QMessageBox.warning(self, "warning", "Nothing to plot. you have to select at least 1 columns in contents.")
      return
    contents = self.contents.getDataFrame()
    names = contents.columns.values.tolist()
    icols = [col.column() for col in columns]
    suptitle = ""
    title = ""
    titles_xy = [names[icol] for icol in icols]
    # xy = [contents[name].values.tolist() for name in titles_xy]
    xy = [contents.loc[:, names[coli]] for coli in icols]
    self.plots.set_n_histograms(xy, suptitle, title, titles_xy)
    self.setCurrentIndex(self.TAB_PLOTS)

  def on_add_column_click(self):
    """add column in dataframe from pandas expression"""
    wid = PandasExpressionDialog("set your pandas expression", self)
    res = wid.exec_()
    if res != wid.Accepted:
      return
    df = self.contents.getDataFrame()
    pdExpr = str(wid.combo.currentText()) # pandas expression of user, no precautions!

    # make alias columns of dataframe, 'C0' as 'C[0]' etc. for user simplicity
    C = [i for i in df.columns]
    namespace = {"C": C, "df": df}
    for i in range(len(C)):
      cmd = "C%i = C[%i]" % (i, i)
      namespace = UXYZ.exec2or3(cmd, namespace)  # TODO python 2-3

    res = None
    try: # expression of user everything could happen !
      namespace = UXYZ.exec2or3(pdExpr, namespace)  # TODO python 2-3
    except:
      trace = traceback.format_exc()
      QtWidgets.QMessageBox.warning(self, "error", "on pandas expression:\n'%s'\n\n%s" % (pdExpr, trace))
      return
    df = namespace["df"]
    if type(df) is not pd.DataFrame:
      QtWidgets.QMessageBox.warning(self, "error", "on pandas expression: df %s is not a DataFrame any more" % type(df))
      return
    else:
      self.contents.setDataFrame(df)
    self.setCurrentIndex(self.TAB_CONTENTS)

"""
# example  pandas 0.23.3
import pandas as pd
import numpy as np
df = pd.DataFrame({
     'X' : ['A', 'A', 'B', np.nan, 'D', 'C'],
     'Y' : [2, 1, 9, 8, 7, 4],
     'Z': [0, 1, 9, 4, 2, 3], })
df
    X    Y    Z
0   A    2    0
1   A    1    1
2   B    9    9
3   NaN  8    4
4   D    7    2
5   C    4    3

Sort by col1

df.sort_values(by=['X'])
    X    Y    Z
0   A    2    0
1   A    1    1
2   B    9    9
5   C    4    3
4   D    7    2
3   NaN  8    4

Sort by multiple columns

df.sort_values(by=['X', 'Y'])
    X    Y    Z
1   A    1    1
0   A    2    0
2   B    9    9
5   C    4    3
4   D    7    2
3   NaN  8    4

# Get list from pandas DataFrame column headers
df.columns.values.tolist() # is the fastest
list(df)
"""

########################################################
class PandasExpressionDialog(QtWidgets.QDialog):
  def __init__(self, title, parent=None):
    super(PandasExpressionDialog, self).__init__(parent)
    self.items = [
      '#pandas expression with "df" as Contents dataFrame...',
      'df["Z1"] = df["X"] * 2        # example add column "Z1"',
      'df["Z2"] = df["X"] * df["Y"]  # example add column "Z2"',
      'df["Z3"] = df["Y"] / df["Z"]  # example add column "Z3 with /0 as missing data"',
      'df = df.sort_values(by="X", ascending=False)     # example sort by column "X"',
      'df = df[["X","Y"]]            # example remove all columns except "X" and "Y"',
      'df = df[df["Z"] > 0]          # example conditional indexing',
      'df = df.dropna()              # example drop rows or columns with missing data',
      'df = df.fillna(aValue)        # example fill missing data with aValue',
      'df = df[(df["Z"] > 0) & (df["Y"] < 2)]      # example conditional indexing',
      'df = df.rename(index=str, columns={"X": "X(cm)"})   # example rename column',
      #'df = pd.DataFrame({"teta": np.linspace(0, np.pi*2, 10)}) #example create ab initio df',
    ]
    
    combo = QtWidgets.QComboBox(self)
    font = QtGui.QFont(self.font()) #copy: self.font() is const
    font.setFamily("Monospace")
    font.setPointSize(9)
    combo.setFont(font)
    combo.addItems(self.items)
    combo.setCurrentIndex(0)
    combo.setEditable(True)
    combo.setDuplicatesEnabled(False)
    self.combo = combo

    dial = QtWidgets.QDialogButtonBox
    buttonBox = dial(dial.Ok | dial.Cancel | dial.Help)
    buttonBox.accepted.connect(self.verify)
    buttonBox.rejected.connect(self.reject)
    buttonBox.helpRequested.connect(self.help)
    
    mainLayout = QtWidgets.QVBoxLayout() #QGridLayout()
    mainLayout.addWidget(combo) #, 0, 0)
    mainLayout.addWidget(buttonBox) #, 3, 0, 1, 3)
    self.setLayout(mainLayout)
    self.setWindowTitle(title)

  def verify(self):
    #http://nullege.com/codes/show/src@p@y@PyQt5-HEAD@examples@richtext@orderform.py/231/PyQt5.QtWidgets.QDialogButtonBox
    self.accept()
    return
      
  def help(self):
    import xyzpy.utilsXyz as UXYZ
    nameBrowser = UXYZ.getBrowser()
    helpPandas1 = "http://pandas.pydata.org/pandas-docs/version/0.23"
    # helpPandas1 = "http://pandas.pydata.org/pandas-docs/version/0.15.2/pandas.pdf"
    helpPandas2 = "http://pandas.pydata.org/index.html" 
    #"http://manishamde.github.io/blog/2013/03/07/pandas-and-python-top-10"
    try:
      os.system("%s %s &" % (nameBrowser, helpPandas1))
      os.system("%s %s &" % (nameBrowser, helpPandas2))
    except:
      pass

          
########################################################
class PandasMainWidgetXyz(QtWidgets.QWidget):
  def __init__(self, parent=None, plotView=None):
    super(PandasMainWidgetXyz, self).__init__(parent)

    self.widget = PandasTabWidget(plotView = plotView)

    self.__createActions()
    self.__addToolBars()
    
    #hbox = QtWidgets.QHBoxLayout()
    #hbox.addWidget(self.button_read)
    hbox = self.toolBars[0]

    # Set the layout
    vbox = QtWidgets.QVBoxLayout()
    #vbox.addLayout(hbox)
    vbox.addWidget(hbox)
    vbox.addWidget(self.widget)

    self.setLayout(vbox)
    self.resize(600,700)
    self.setWindowTitle("PandasMainWidget")
    self._initialValue = None 
    
  def __addToolBars(self):
    self.toolBars = []
    tb = QtWidgets.QToolBar("Edit")   #self.addToolBar("Edit")
    for action in self.actions:
      tb.addAction(action)
    #act = ACFX.getCommonActionByName("GeneralHelp")
    #if act != None: tb.addAction(act)
    self.toolBars.append(tb)

  def __createCommonsActions(self):
    """create actions for self widget AND other widgets through ACFX.addInCommonActions"""
    logger.debug("create actions %s" % (self.objectName()))
    
    action = ACFX.QActionXyz(name="ReadFile", text="Read .csv")
    ok = action.setAction( slot=self.widget.on_read_click, signal=aSignal,
                           shortcut=None, tooltip="Read file csv", icon="open" )
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)
    
    action = ACFX.QActionXyz(name="Plot1", text="Plot y=f(x)")
    ok = action.setAction( slot=self.widget.on_plot1_click, signal=aSignal,
                           shortcut=None, tooltip="plot y=f(x) 2 selected colums", icon="plot" )
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    action = ACFX.QActionXyz(name="Plot2", text="Plot y1,y2=f(x)")
    ok = action.setAction( slot=self.widget.on_plot2_click, signal=aSignal,
                           shortcut=None, tooltip="plot y=f(x) 3 selected colums", icon="plot2" )
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    action = ACFX.QActionXyz(name="Plotn", text="Plot y1,...yn=f(x)")
    ok = action.setAction( slot=self.widget.on_plotn_click, signal=aSignal,
                           shortcut=None, tooltip="plot y1,...yn=f(x) n+1 selected colums", icon="plotn" )
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    action = ACFX.QActionXyz(name="DistinctPlots ", text="Distinct Plot y1,...yn=f(x)")
    ok = action.setAction( slot=self.widget.on_plotd_click, signal=aSignal,
                           shortcut=None, tooltip="distinct plot y1,...yn=f(x) n+1 selected colums", icon="plotd" )
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    action = ACFX.QActionXyz(name="HistogramPlot ", text="Histogram Plot")
    ok = action.setAction( slot=self.widget.on_ploth_click, signal=aSignal,
                           shortcut=None, tooltip="histogram plot 1 or more selected colums", icon="ploth" )
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    mode = getMode()

    if mode == "advanced":
      action = ACFX.QActionXyz(name="ExamplePlot", text="Plot example ellipse")
      ok = action.setAction( slot=self.widget.on_example_plot1_click, signal=aSignal,
                            shortcut=None, tooltip="plot example ellipse", icon="plot" )
      if ok:
        ACFX.addInCommonActions(action)
        self.actions.append(action)

    action = ACFX.QActionXyz(name="AddColumn", text="execute pandas expression")
    ok = action.setAction( slot=self.widget.on_example_plot1_click, signal=aSignal,
                           shortcut=None, tooltip="execute pandas expression", icon="addColumn" )
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

  def __createActions(self):
    """create local widget actions"""
    self.actions = []
    self.actions.append(self.__createAction( 
      "ReadFile", None, "Read file csv", self.widget.on_read_click, "opencsv" ) )
    self.actions.append(self.__createAction( 
      "Plot1", None, "Plot y=f(x) 2 selected columns", self.widget.on_plot1_click, "plot" ) )
    self.actions.append(self.__createAction( 
      "Plot2", None, "Plot yi=f(x) 3 selected columns (x, y1, y2)", self.widget.on_plot2_click, "plot2" ) )
    self.actions.append(self.__createAction( 
      "Plotn", None, "Plot yn=f(x) n+1 selected columns (x, y1,.. yn)", self.widget.on_plotn_click, "plotn" ) )
    self.actions.append(self.__createAction( 
      "Plotd", None, "distinct Plot yn=f(x) n+1 selected columns (x, y1,.. yn)", self.widget.on_plotd_click, "plotd" ) )
    self.actions.append(self.__createAction( 
      "Ploth", None, "histogram Plot 1 or more selected columns", self.widget.on_ploth_click, "ploth" ) )
    self.actions.append(self.__createAction( 
      "AddColumn", None, "execute pandas expression", self.widget.on_add_column_click, "addColumn" ) )

    mode = getMode()

    if mode == "advanced":
      self.actions.append(self.__createAction(
        "PlotEllipse", None, "Plot example ellipse", self.widget.on_example_plot1_click, "plotEllipse" ) )
    #etc...

  def __createAction(self, Name, Shortcut, ToolTip, Call, Icon=None, Enable=True):
    action = QtWidgets.QAction(Name, self)
    if Shortcut!=None: action.setShortcut(self.prefixShortcut+Shortcut)
    action.setToolTip(ToolTip)
    if Icon!=None:
      action.setIcon(IUSR.getIconFromName(Icon))
    action.setEnabled(Enable)
    action.triggered.connect(Call)
    return action
    
  def on_help(self):
    """help slot"""
    import xyzpy.utilsXyz as UXYZ
    nameBrowser = UXYZ.getBrowser()  #
    helpPandas = "http://pandas.pydata.org/pandas-docs/version/0.23"
    #helpPandas = "http://pandas.pydata.org/pandas-docs/version/0.15.2/pandas.pdf"
    #helpPandas = "http://pandas.pydata.org/index.html" 
    try:
      #os.system("%s %s &" % (nameBrowser, helpPandas)) #pdf
      os.system("%s %s &" % (nameBrowser, helpPandas)) #html
    except:
      pass
  
  def getValue(self):
    print("PandasMainWidgetXyz.getValue")
    res = self.widget.contents.getDataFrame()
    return res

  def setValue(self, value):
    if debug: 
      print("PandasMainWidgetXyz.setValue:")
      print(value)
    self.widget.contents.setDataFrame(value)
    self.widget.setInfos()
    if self._initialValue is None: 
      if debug: print("PandasMainWidgetXyz.setValue set initialValue:\n", value)
      self._initialValue = value.copy()
    return

  def resetValue(self):
    if debug: print("PandasMainWidgetXyz.resetValue")
    self.setValue(self._initialValue)
    return
 

########################################################
class PandasMainDialogXyz(DialogXyz):

  def __init__(self, title="pandas plot", parent=None):
    super(PandasMainDialogXyz, self).__init__(parent)
    self.setUpWidgetLayout(PandasMainWidgetXyz())
    self.setWindowTitle(title)
    #self.buttons["Reset"].hide() #if not needed...?
    
  def setValue(self, value):
    #if debug: 
    #  print "PandasMainDialogXyz.setValue:"
    #  print value
    self.getUpWidgetLayout().setValue(value)
    return

  def setValueFromFile(self, aNameFile, typeCsv="FirstCol"):
    self.getUpWidgetLayout().widget.read_csv(aNameFile, typeCsv=typeCsv)
      
  def getValue(self):
    if debug: print("PandasMainDialogXyz.getValue")
    res = self.getUpWidgetLayout().getValue()
    return res
       


if __name__ == '__main__':
  app = QtWidgets.QApplication(sys.argv)
  #mw = PandasMainWidgetXyz()
  mw = PandasMainDialogXyz("Pandas Viewer")
  mw.setValue(_example_DataFrame)
  mw.show()
  print(mw.getValue())
  app.exec_()
  print(mw.choice)

