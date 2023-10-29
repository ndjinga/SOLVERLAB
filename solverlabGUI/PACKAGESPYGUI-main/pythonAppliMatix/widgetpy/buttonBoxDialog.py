#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

import os
import sys

from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *

import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()

verbose = False

###############################################################
class MyButton(QPushButton):

    def __init__(self, text):
        QPushButton.__init__(self, text)
        self.setMouseTracking(True)
        self.initialStyleSheet = "background-color: rgb(200, 200, 200);"

    def SetStyleSheet(self, s):
        self.setStyleSheet(s)
        self.initialStyleSheet = s

    def mouseMoveEvent(self, event):
        self.setFocus()

    def focusInEvent(self, event):
        self.setStyleSheet("background-color: yellow;")

    def focusOutEvent(self, event):
        self.setStyleSheet(self.initialStyleSheet)


###############################################################
class ButtonBoxDialog(QDialog):
    """
    quick QComboBox with buttons, if items are files with path, display basename only
    if no choice and close dialog, ButtonBoxDialog.choice is None
    """
    def __init__(self, title, items, parent=None):
        super(ButtonBoxDialog, self).__init__(parent)
        self.title = title
        self.setWindowTitle(title)
        self.choice = None
        self.indice = None
        self.setWidgets(items)
        self.closeOnChoice = True
        if verbose: logger.info("ButtonBoxDialog.init: %s" % items)
        
    def setWidgets(self, items):
        self.mainLayout = QVBoxLayout(self)
        self._buttons = []
        self.items = {}
        for indice, i in enumerate(items):
          path, name = os.path.split(i)
          if name in list(self.items.keys()):
            name = i #if same basename display 2nd with path
          aButton = MyButton(name) #text of button
          self.mainLayout.addWidget(aButton)
          self.items[name] = i
          aButton.clicked.connect( lambda status=None, ele=aButton, ind=indice: self.setChoiceButton(status, ele, ind) )

        self.show()
        w = max(len(self.title)*8+50, 40)
        self.setMinimumWidth(w)

    def setChoice(self, text):
        self.choice = text

    def setIndice(self, text):
        self.indice = text

    def setChoiceButton(self, status, aButton, indice):
        choice = self.items[str(aButton.text())]
        if verbose: logger.info("ButtonBoxDialog choice '%s' indice %i" % (choice, indice))
        self.setChoice(choice)
        self.setIndice(indice)
        if self.closeOnChoice:
          self.close()


def test():
    app = QApplication(sys.argv)
    ex = ButtonBoxDialog("Select your choice", ["/tmp/Choice_alonglonglongname1","/tmp1/Choice_2","/tmp/Choice_2","/tmp/Choice_3"])
    ex.exec_()
    print("choice", ex.choice)

if __name__ == '__main__':
    test()

