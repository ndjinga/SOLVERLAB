#!/usr/bin/env python
# -*- coding: utf-8 -*-

#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

import sys
from PyQt5 import QtCore, QtGui, QtWidgets
from salomepy.onceQApplication import OnceQApplication


class StyleWidget(QtWidgets.QDialog):
  def __init__(self, parent=None):
    super(StyleWidget, self).__init__(parent)
    self._styles = list(QtWidgets.QStyleFactory.keys())
    layout = QtWidgets.QVBoxLayout()
    self.styleLabel = QtWidgets.QLabel("Choose Style")
    self.styleComboBox = QtWidgets.QComboBox()
    self.styleComboBox.addItems(self._styles) # adding the styles list
    # find current style
    index = self.styleComboBox.findText(
              QtWidgets.qApp.style().objectName(),
              QtCore.Qt.MatchFixedString )
    # set current style
    self.styleComboBox.setCurrentIndex(index)  
    layout.addWidget(self.styleLabel)
    layout.addWidget(self.styleComboBox)
    self.setLayout(layout)
    self.styleComboBox.activated[str].connect(self.on_StyleChanged)
    self.setFixedSize(250,100)

  def getCurrentIndexStyle(self): 
    """get the current index from combobox"""
    return self.styleComboBox.currentIndex()
  
  def getCurrentStyle(self):
    """get the current name from combobox"""
    return self._styles[self.getCurrentIndexStyle()]

  def on_StyleChanged(self, style):
    """slot when changing style"""
    app = OnceQApplication()
    app.setStyle(style)

def run():
  #app = QtWidgets.QApplication(sys.argv)
  app = OnceQApplication()
  wid = StyleWidget()
  wid.show()
  app.exec_()
  print("current style: %s" % wid.getCurrentStyle())

if __name__ == "__main__":
  run()
