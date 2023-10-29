#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END

import os
import sys
from PyQt5.QtWidgets import QComboBox, QDialog, QApplication

verbose = False

class ComboBoxDialog(QDialog):
    """
    quick QComboBox, if items are files with path, display basename only
    if no choice and close dialog, ComboBoxDialog.choice is None
    """
    def __init__(self, title, items, parent=None):
        super(ComboBoxDialog, self).__init__(parent)
        self.title = title
        self.setWindowTitle( title )
        self.choice = None
        self.setWidgets(items)
        if verbose: print("ComboBoxDialog.init: ",items)
        self.combo.activated[str].connect(self.onActivated)
        
    def setWidgets(self, items):
        self.combo = QComboBox(self)
        self.items = {}
        for i in items:
          path , name = os.path.split(i)
          if name in list(self.items.keys()):
            name = i #if same basename display 2nd with path
          self.combo.addItem(name)
          self.items[name] = i
        self.combo.setSizeAdjustPolicy( QComboBox.AdjustToContents )
        self.combo.move(10, 10)
        self.show()
        w = max(len(self.title)*8+50, self.combo.width()+20)
        self.setMinimumWidth(w)

    def onActivated(self, text):
        if verbose: "ComboBoxDialog.choice", text
        self.choice=self.items[str(text)]
        self.close()

def test():
    app = QApplication(sys.argv)
    ex = ComboBoxDialog("Select your choice", ["/tmp/Choice_alonglonglongname1","/tmp1/Choice_2","/tmp/Choice_2","/tmp/Choice_3"])
    ex.exec_()
    print("choice",ex.choice)

if __name__ == '__main__':
    test()

