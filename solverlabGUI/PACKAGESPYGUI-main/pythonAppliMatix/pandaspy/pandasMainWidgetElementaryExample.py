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
"""

import numpy as np
import pandas as pd
from pandas.sandbox.qtpandas import DataFrameModel, DataFrameWidget

#from PySide import QtGui, QtCore
# Or if you use PyQt5:
from PyQt5 import QtGui, QtCore, QtWidgets

class PandasMainWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(PandasMainWidget, self).__init__(parent)

        # Create two DataFrames
        self.df1 = pd.DataFrame(np.arange(9).reshape(3, 3),
                                columns=['foo', 'bar', 'baz'])
        self.df2 = pd.DataFrame({
                'int': [1, 2, 3],
                'float': [1.5, 2.5, 3.5],
                'string': ['a', 'b', 'c'],
                'nan': [np.nan, np.nan, np.nan]
            }, index=['AAA', 'BBB', 'CCC'],
            columns=['int', 'float', 'string', 'nan'])

        # Create the widget and set the first DataFrame
        self.widget = DataFrameWidget(self.df1)

        # Create the buttons for changing DataFrames
        self.button_first = QtWidgets.QPushButton('First')
        self.button_first.clicked.connect(self.on_first_click)
        self.button_second = QtWidgets.QPushButton('Second')
        self.button_second.clicked.connect(self.on_second_click)

        # Set the layout
        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(self.widget)
        hbox = QtWidgets.QHBoxLayout()
        hbox.addWidget(self.button_first)
        hbox.addWidget(self.button_second)
        vbox.addLayout(hbox)
        self.setLayout(vbox)

    def on_first_click(self):
        '''Sets the first DataFrame'''
        self.widget.setDataFrame(self.df1)

    def on_second_click(self):
        '''Sets the second DataFrame'''
        self.widget.setDataFrame(self.df2)

if __name__ == '__main__':
    import sys

    # Initialize the application
    app = QtWidgets.QApplication(sys.argv)
    mw = PandasMainWidget()
    mw.show()
    app.exec_()

