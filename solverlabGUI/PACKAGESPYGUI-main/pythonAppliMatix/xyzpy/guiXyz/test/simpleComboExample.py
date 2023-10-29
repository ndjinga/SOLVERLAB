#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""https://stackoverflow.com/questions/36614635/pyqt-right-click-menu-for-qcombobox"""

import sys
from PyQt5 import QtGui, QtWidgets
from PyQt5.QtCore import Qt, QVariant
from PyQt5.QtGui import QMenu

class Example(QtWidgets.QWidget):

    def __init__(self):
        super(Example, self).__init__()
        self.initUI()

    def initUI(self):
        self.lbl = QtWidgets.QLabel("Ubuntu", self)
        self.combo = QtWidgets.QComboBox(self)
        self.combo.setContextMenuPolicy(Qt.CustomContextMenu)
        self.combo.customContextMenuRequested.connect(self.showMenu)
        self.combo.addItem("Ubuntu")
        self.combo.addItem("Mandriva")
        self.combo.addItem("Fedora")
        self.combo.addItem("Red Hat")
        self.combo.addItem("Gentoooooooooooooooooo")
        self.combo.move(10, 30)
        self.lbl.move(10, 10)
        self.combo.activated[str].connect(self.onActivated)
        self.setGeometry(300, 300, 300, 200)
        self.setWindowTitle('QtWidgets.QComboBox')
        self.show()

    def showMenu(self,pos):
        menu = QMenu()
        clear_action = menu.addAction("Clear Selection")
        action = menu.exec_(self.mapToGlobal(pos))
        if action == clear_action:
            self.combo.setCurrentIndex(0)

    def onActivated(self, text):
        self.lbl.setText(text)
        self.lbl.adjustSize()

def main1():
    app = QtWidgets.QApplication(sys.argv)
    ex = Example()
    ex.show()
    app.exec_()

def main2():
    """https://wiki.python.org/moin/PyQt/Using%20a%20different%20view%20with%20QComboBox"""
    app = QtWidgets.QApplication(sys.argv)
    model = QtGui.QStandardItemModel()
    items = [("ABC", True),("DEF", False),("GHI", False)]
    for text, checked in items:
      text_item = QtGui.QStandardItem(text)
      checked_item = QtGui.QStandardItem()
      checked_item.setData(QVariant(checked), Qt.CheckStateRole)
      model.appendRow([text_item, checked_item])
    view = QtWidgets.QTreeView()
    view.header().hide()
    view.setRootIsDecorated(False)
    combo = QtWidgets.QComboBox()
    combo.setView(view)
    combo.setModel(model)
    combo.setGeometry(300, 300, 300, 200)
    combo.show()
    app.exec_()

if __name__ == '__main__':
    main1()    
    main2()
