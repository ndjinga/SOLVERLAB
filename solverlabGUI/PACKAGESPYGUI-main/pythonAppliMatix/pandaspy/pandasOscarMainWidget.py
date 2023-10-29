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
"""

import os, sys
import traceback

import numpy as np
import pandas as pd

#import xyzpy.actionsFactoryXyz as ACFX
import xyzpy.utilsXyz as UXYZ
import salomepy.iconsUser as IUSR
import pandaspy.utilsPandasOscar as UPO
from matplotlibpy.matplotlibWindowToolbar import MatplotlibWindowToolbar

verbose = False

from PyQt5 import QtGui, QtCore, QtWidgets


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
    self._df = dataFrame

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

    #return QtCore.QVariant(str(self._df.loc[index.row(), index.column()])) #bug if conditional indexing
    return QtCore.QVariant(str(self._df.iloc[index.row(), index.column()]))

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
    
  def read_oscar_csv(self, aNameFile):
    try:
      self._df1, self._df2 = UPO.read_oscar_csv(aNameFile)
      self.infos.setDataFrame(self._df1)
      self.contents.setDataFrame(self._df2)
      self.setWindowTitle(os.path.basename(aNameFile))
      #self.infos.resizeColumnsToContents()
      #self.contents.resizeColumnsToContents()
    except Exception as e:
      trace = traceback.format_exc() #better explicit verbose problem
      mess = "Problem reading file\n'%s'\n%s" % (aNameFile, trace)
      QtWidgets.QMessageBox.warning(self, "error", mess)
    return
    
  def on_read_click(self):
    """read a csv in DataFrame"""
    aDir = os.getenv("OSCAR_WORKDIR")
    if aDir == None:
      if self.messDone:
        QtWidgets.QMessageBox.warning(self, "warning", "You could set $WORKDIR4OSCAR as your directory of results OSCAR, files csv.")
        self.messDone = true
      aDir = os.path.split(__file__)[0]
      #aDir = os.getenv("HOME")
    #an example from ResulatsOSCAR_IHMSortie
    #aDir = os.path.join(os.path.abspath(os.path.split(__file__)[0]), "test")
    #aFile = "Terme_Source_Chaine_PF_11.csv"
    #aNameFile = os.path.join(aDir, aFile)
    aNameFile = str(QtWidgets.QFileDialog.getOpenFileName(self, 'Load file .csv', aDir, "(*.csv *.CSV)")[0])
    if aNameFile == '': return
    self.read_oscar_csv(aNameFile)

  def on_example_plot1_click(self):
    """plot from 2 DataFrames columns"""
    df_infos, df_contents = UPO.getExampleDataFrameEllipse(periods=1000)
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
      QtWidgets.QMessageBox.warning(self, "warning", "Nothing to plot. you have to select 2 columns in contents.")
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
    if res != wid.Accepted: return
    df = self.contents.getDataFrame()
    pdExpr = str(wid.combo.currentText())  # pandas expression of user, no precautions!

    # make alias columns of dataframe, 'C0' as 'C[0]' etc. for user simplicity
    C = [i for i in df.columns]
    namespace = {"C": C, "df": df}
    for i in range(len(C)):
      cmd = "C%i = C[%i]" % (i, i)
      namespace = UXYZ.exec2or3(cmd, namespace)  # TODO python 2-3

    res = None
    try:  # expression of user everything could happen !
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


class PandasExpressionDialog(QtWidgets.QDialog):
    def __init__(self, title, parent=None):
        super(PandasExpressionDialog, self).__init__(parent)
        self.items = [
          '#pandas expression with "df" as Contents dataFrame...',
          'df["Z1"] = df.X * 2           #example add column "Z1"',
          'df["Z1"] = df[[C1]] * 2       #same example add column "Z1"',
          'df["Z1"] = df[["X"]] * 2      #same example add column "Z1"',
          'df["Z2"] = df.X * df.Y        #example add column "Z2"', 
          'df["Z3"] = df.X / df.Y        #example add column "Z3 with missing data"', 
          'df = df.sort(columns="X")     #example sort by column "X"', 
          'df = df[["X","Y"]]            #example remove all columns except "X" and "Y"',
          'df = df[[C1,C2]]              #same example remove all columns except "X" and "Y"',
          'df = df[df["X"] > 0]          #example conditional indexing',
          'df = df.dropna()              #example drop rows or columns with missing data',
          'df = df["Z3"].fillna(aValue)  #example fill missing data with aValue',
          'df = df[(df["Z3"] > 0) & (df["Z3"] < 2)]      #example conditional indexing',
          'df = df.rename(columns={"Y": "otherNameY"})   #example rename column',
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
        nameBrowser = UXYZ.getBrowser() #
        helpPandas1 = "http://pandas.pydata.org/pandas-docs/version/0.23/pandas.pdf"
        helpPandas2 = "http://pandas.pydata.org/index.html" 
        #"http://manishamde.github.io/blog/2013/03/07/pandas-and-python-top-10"
        try:
          os.system("%s %s &" % (nameBrowser, helpPandas1))
          os.system("%s %s &" % (nameBrowser, helpPandas2))
        except:
          pass
            
class PandasMainWidget(QtWidgets.QWidget):
  def __init__(self, parent=None, plotView = None):
    super(PandasMainWidget, self).__init__(parent)

    self.widget = PandasTabWidget( plotView = plotView)
    #self.button_read = QtWidgets.QPushButton('Read OSCAR .csv')
    #self.button_read.clicked.connect(self.widget.on_read_click)

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
    
  def __addToolBars(self):
    self.toolBars = []
    tb = QtWidgets.QToolBar("EditOscar")   #self.addToolBar("Edit")
    for action in self.actions:
      tb.addAction(action)
    #act = ACFX.getCommonActionByName("GeneralHelp")
    #if act != None: tb.addAction(act)
    self.toolBars.append(tb)

  def __createCommonsActions(self):
    """create actions for self widget AND other widgets through ACFX.addInCommonActions"""
    logger.debug("createCommonsActions %s" % (self.objectName()))
    
    action = ACFX.QActionXyz(name="ReadOscarFile", text="Read OSCAR .csv")
    ok = action.setAction( slot=self.widget.on_read_click, signal=aSignal,
                           shortcut=None, tooltip="Read OSCAR file csv", icon="open" )
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
      "ReadOscarFile", None, "Read OSCAR file csv", self.widget.on_read_click, "opencsv" ) )
    self.actions.append(self.__createAction( 
      "Plot1", None, "Plot y=f(x) 2 selected colums", self.widget.on_plot1_click, "plot" ) )
    self.actions.append(self.__createAction( 
      "Plot2", None, "Plot yi=f(x) 3 selected colums (x, y1, y2)", self.widget.on_plot2_click, "plot2" ) )
    self.actions.append(self.__createAction( 
      "Plotn", None, "Plot yn=f(x) n+1 selected colums (x, y1,.. yn)", self.widget.on_plotn_click, "plotn" ) )
    self.actions.append(self.__createAction( 
      "Plotd", None, "distinct Plot yn=f(x) n+1 selected colums (x, y1,.. yn)", self.widget.on_plotd_click, "plotd" ) )
    self.actions.append(self.__createAction( 
      "AddColumn", None, "execute pandas expression", self.widget.on_add_column_click, "addColumn" ) )
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


if __name__ == '__main__':
  app = QtWidgets.QApplication(sys.argv)
  mw = PandasMainWidget()
  mw.show()
  app.exec_()

