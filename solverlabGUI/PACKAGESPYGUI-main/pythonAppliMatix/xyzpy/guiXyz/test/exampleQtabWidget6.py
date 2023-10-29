#!/usr/bin/env python
# -*- coding: utf-8 -*-


from PyQt5 import QtGui, QtWidgets
from PyQt5 import QtCore
import sys


mystruct = {'A': ['aa', 'ab', 'ac'],
            'B': ['ba', 'bb', 'bc'],
            'C': ['ca', 'cb', 'cc']}


##########################################################
class MyTable(QtWidgets.QTableWidget):
  def __init__(self, theStruct, *args):
    super(MyTable, self).__init__(*args)
    self.data = theStruct
    self.setmydata()

  def setmydata(self):
    n = 0
    for key in sorted(self.data.keys()):
      m = 0
      for item in self.data[key]:
        newitem = QtWidgets.QTableWidgetItem(item)
        self.setItem(m, n, newitem)
        m += 1
      n += 1


##########################################################
class ExampleWid(QtWidgets.QWidget):
  no = 0
  def __init__(self, *args, **kwargs):
    super(ExampleWid, self).__init__(*args, **kwargs)
    noc = ExampleWid.no
    self.setObjectName("tab_%i" % noc)
    ExampleWid.no += 1
    #print self.objectName()
    vBoxlayout  = QtWidgets.QVBoxLayout()
    self.setLayout(vBoxlayout)
    
    ii = 0
    for i in range(1):
      groupBox = QtWidgets.QGroupBox("groupBox_%i_%i" % (noc, ii))
      vBoxlayout.addWidget(groupBox)
      vbox =  QtWidgets.QVBoxLayout()
      checks = [QtWidgets.QRadioButton(name) for name in ["bonjour","hello"]]
      for c in checks: vbox.addWidget(c)
      groupBox.setLayout(vbox)
      ii += 1
    
    for i in range(2):
      groupBox = QtWidgets.QGroupBox("groupBox_%i_%i" % (noc, ii))
      vBoxlayout.addWidget(groupBox)
      vbox =  QtWidgets.QVBoxLayout()
      tables = [MyTable(mystruct,5,3) for i in range(1)]
      for t in tables: vbox.addWidget(t)
      groupBox.setLayout(vbox)
      ii += 1


##########################################################
class ExampleTab(QtWidgets.QWidget):
  def __init__(self, tabs=[]):
    super(ExampleTab, self).__init__()
    
    tabsWidget = QtWidgets.QTabWidget()
    for t in tabs:
      tabsWidget.addTab(t, t.objectName())
      
    vBoxlayout  = QtWidgets.QVBoxLayout()
    hBoxlayout  = QtWidgets.QHBoxLayout()
    pushButtons = [QtWidgets.QPushButton(name) for name in ["Apply", "Reset", "Cancel"]]
    vBoxlayout.addWidget(tabsWidget)
    for pb in pushButtons: hBoxlayout.addWidget(pb)
    vBoxlayout.addLayout(hBoxlayout)
    
    self.setLayout(vBoxlayout)
    self.setWindowTitle('aWidgetWithTabs')

    self.resize(400, 450)
    self.move(300, 300) #Move QTabWidget to x:300,y:300
    tabsWidget.currentChanged.connect(self.tabChangedSlot)
    
  def tabChangedSlot(self, argTabIndex):
    QtWidgets.QMessageBox.information(self,"Tab Index Changed!",
          "Current Tab Index: "+str.number(argTabIndex));


##########################################################
if __name__ == '__main__':
  app = QtWidgets.QApplication(sys.argv)
  print([str(i) for i in list(QtWidgets.QStyleFactory.keys())])
  app.setStyle("Cleanlooks") #['Oxygen', 'Windows', 'Motif', 'CDE', 'Plastique', 'GTK+', 'Cleanlooks']
  tabs = [ExampleWid() for i in range(4)]
  wid = ExampleTab(tabs)
  wid.show()
  sys.exit(app.exec_())



