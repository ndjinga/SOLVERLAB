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
import numpy as np
import pandas as pd
from PyQt5 import QtGui, QtCore, QtWidgets

"""
qtpandas.py:10: UserWarning: The pandas.sandbox.qtpandas module is deprecated and will be removed 
in a future version. We refer users to the external package here: 
https://github.com/datalyze-solutions/pandas-qt

from pandas.sandbox.qtpandas import DataFrameModel #, DataFrameWidget
import pandas.sandbox.qtpandas as qtpd
print qtpd.__file__
"""

verbose = False

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


########################################################
class DataFrameWidget(QtWidgets.QWidget):
    ''' a simple widget for using DataFrames in a gui '''
    #_defaut_dataFrame = pd.DataFrame(np.arange(0))
    _defaut_dataFrame = pd.DataFrame({})
    
    def __init__(self, dataFrame={}, parent=None):
        super(DataFrameWidget, self).__init__(parent)
        self.dataModel = DataFrameModel()
        self.dataTable = QtWidgets.QTableView()
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.dataTable)
        self.setLayout(layout)
        self.setDataFrame(dataFrame)

    def resizeColumnsToContents(self):
        self.dataTable.resizeColumnsToContents()
        
    def setDataFrame(self, dataFrame={}):
        if type(dataFrame) == dict:
          data = pd.DataFrame(dataFrame)
        else:
          data = dataFrame
        self.dataModel.setDataFrame(data)
        self.dataTable.setModel(self.dataModel)
        self.dataModel.signalUpdate()


########################################################
class PandasMainWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(PandasMainWidget, self).__init__(parent)

        #self.widget = DataFrameWidget(pd.DataFrame(np.arange(0)))
        self.widget = DataFrameWidget()

        self.button_read = QtWidgets.QPushButton('Read OSCAR .csv')
        self.button_read.clicked.connect(self.on_read_click)

        # Set the layout
        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(self.widget)
        hbox = QtWidgets.QHBoxLayout()
        hbox.addWidget(self.button_read)
        vbox.addLayout(hbox)
        self.setLayout(vbox)
    
    def getHeader(self, aNameFile, delimiter=";", nbLinesMax=20):
        """
        get length of header and header as dict.
        length is numer of lines to skip before DataFrameModel as csv
        check info in header as first lines (no more nbLinesMax)
        with no more than one or two delimiter
        """
        #print firsts line... check info in header as first lines with on or two delimiter
        aDict = {}
        nbHeader = 0
        if verbose: print("open '%s'" % aNameFile)
        with open(aNameFile) as myFile:
          for x in range(nbLinesMax):
            line = next(myFile)
            if verbose: print("--",line, end=' ')
            lineSplit = line.split(delimiter)
            nbCols = len(lineSplit)
            if nbCols > 3:
              line = next(myFile)
              print("--",line, end=' ')
              return (nbHeader, aDict, nbCols-1)
            nbHeader += 1
            aDict[lineSplit[0]] = lineSplit[1]
          print("Problem reading header: more than %i lines" % nbLinesMax)
          return (0, {}, 0)

    def on_read_click(self):
        """read a csv in DataFrame"""
        #aDir = "./oscar_original/Livraison_DM2S/ResulatsOSCAR_IHMSortie/Terme_Source/Terme_Source_PF"
        #aFile = "Terme_Source_Chaine_PF_11.csv" #an example from ResulatsOSCAR_IHMSortie
        #aNameFile = os.path.join(aDir, aFile)
        aDir = os.path.join(os.path.abspath(os.path.split(__file__)[0]), "test")
        aNameFile = str(QtWidgets.QFileDialog.getOpenFileName(self, 'Load file .csv', aDir, "(*.csv *.CSV)")[0])
        if aNameFile == '': return
        delimiter = ";"
        nbHeader, aDict, nbCols = self.getHeader(aNameFile, delimiter=delimiter)
        if verbose: print("info: %s" % aDict) #nbHeader, nbCols, aDict
        
        self._df = pd.read_csv(aNameFile, delimiter=delimiter, header=nbHeader, usecols=list(range(nbCols)))
        self.widget.setDataFrame(self._df)
        self.setWindowTitle(os.path.basename(aNameFile))
        self.widget.resizeColumnsToContents()


########################################################
if __name__ == '__main__':
    import sys
    verbose = True
    app = QtWidgets.QApplication(sys.argv)
    mw = PandasMainWidget()
    mw.show()
    app.exec_()

