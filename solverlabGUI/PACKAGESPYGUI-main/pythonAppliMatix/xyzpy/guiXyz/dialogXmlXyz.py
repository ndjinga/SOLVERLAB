#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import traceback
import pprint as PP

from PyQt5 import QtCore, QtGui, QtWidgets
import xml.etree.ElementTree as ET
import xyzpy.utilsXyz as UXYZ
import xyzpy.baseXyz as BXYZ
import xyzpy.intFloatListXyz as IFLX
import xyzpy.classFactoryXyz as CLFX
import xyzpy.helpsFactoryXyz as HLFX
import salomepy.iconsUser as IUSR
from xyzpy.controllerXyz import ControllerXyz #local model through local controller
import xyzpy.loggingXyz as LOG

from salomepy.strEvent import *
from time import sleep

logger = LOG.getLogger()

verbose = False
verboseEvent = verbose
verboseHelp = verbose
debug = verbose


##########################################################
class ExampleWid(QtWidgets.QWidget):
  
  _no = 0
  
  def __init__(self, *args, **kwargs):
    super(ExampleWid,self).__init__(*args, **kwargs)
    noc = self._no
    self.setObjectName("tab_%i" % noc)
    self._no += 1
    if verbose: logger.info("ExampleWid.__init__ %s" % self.objectName())
    vBoxlayout  = QtWidgets.QVBoxLayout()
    self.setLayout(vBoxlayout)

    ii = 0
    for i in range(1):
      groupBox = QtWidgets.QGroupBox("groupBox_%i_%i" % (noc, ii))
      vBoxlayout.addWidget(groupBox)
      vbox =  QtWidgets.QVBoxLayout()
      checks = [QtWidgets.QRadioButton(name) for name in ["example_bonjour","example_hello"]]
      for c in checks: 
        vbox.addWidget(c)
      groupBox.setLayout(vbox)
      ii += 1

##########################################################
class LayoutWidgetXyz(QtWidgets.QWidget):
  
  _no = 0
  
  def __init__(self, *args, **kwargs):
    super(LayoutWidgetXyz,self).__init__(*args, **kwargs)
    self.noc = self._no
    self.setObjectName("LayoutWidgetXyz_%i" % self.noc)
    self._no += 1
    if verbose: logger.info("LayoutWidgetXyz.__init__ %" % self.objectName())
    vBoxlayout  = QtWidgets.QVBoxLayout()
    self.setLayout(vBoxlayout)

  def setGroupBox(self, key, value):
    groupBox = QtWidgets.QGroupBox(key)
    self.layout().addWidget(groupBox)
    vbox =  QtWidgets.QVBoxLayout()
    #checks = [QtWidgets.QRadioButton(name) for name in ["bonjour","hello"]]
    #for c in checks: vbox.addWidget(c)
    groupBox.setLayout(vbox)


##########################################################
class TabWidgetXyz(QtWidgets.QWidget):
  
  _no = 0
  
  def __init__(self, tabs=[]):
    super(TabWidgetXyz,self).__init__()
    self.noc = self._no
    tabsWidget = QtWidgets.QTabWidget()
    for t in tabs:
      if verbose: logger.info("TabWidgetXyz.addTab %s" % t.objectName())
      tabsWidget.addTab(t, t.objectName())
    vBoxlayout  = QtWidgets.QVBoxLayout()
    vBoxlayout.addWidget(tabsWidget)
    self.setLayout(vBoxlayout)
    self._margin = 2
    m = self._margin
    self.layout().setContentsMargins(m,m,m,m) #left, top, right, bottom
    self.setWindowTitle("TabWidgetXyz_%i" % self.noc)
    self._no += 1
    tabsWidget.currentChanged.connect(self.tabChangedSlot)

  def tabChangedSlot(self, argTabIndex):
    if debug: logger.info("Current Tab Index: %s" % argTabIndex)
    #QtWidgets.QMessageBox.information(self,"Tab Index Changed","Tab Index: %s" % argTabIndex)


###############################################################
class QTextEditXyz(QtWidgets.QTextEdit):

  def __init__(self, *args, **kwargs):
    super(QtWidgets.QTextEdit, self).__init__(*args, **kwargs)
    self.setFont(QtGui.QFont("Monospace", 9))
    self.initialValue = str(self.toPlainText())

  def setValue(self, value):
    if verbose: logger.info("QTextEditXyz.setValue %s" % value)
    self.setText(value)
    self.initialValue = str(value)

  def getValue(self):
    """
    method implemented to avoid have to use equivalent methods with other names:
    
    - QLineEdit.text()
    - QComboBox.currentText()
    - QSpinBox/QDoubleSpinBox.value()
    - etc...
    """
    return str(self.toPlainText())

  def resetValue(self):
    self.setText(self.initialValue)
    

##########################################################
class DialogXyz(QtWidgets.QDialog):
  """
  a vertical layout with a upWidget with a QVBoxLayout
  and at bottom: Apply, Reset, Cancel, Help buttons
  """

  _no = 0

  def __init__(self, *args, **kwargs):
    super(DialogXyz,self).__init__(*args, **kwargs)
    noc = self._no
    self.setObjectName("DialogXyz_%i" % noc)
    self._className = str(self.__class__.__name__)
    self.setWindowTitle(self.objectName())
    self._no += 1

    self._controller = None #caller controller general model
    
    # for upwidget, and buttons
    vBoxlayout = QtWidgets.QVBoxLayout()
    self.setLayout(vBoxlayout)

    m = 2
    self._margin = m
    vBoxlayout.setContentsMargins(m,m,m,m*2) #left, top, right, bottom

    # buttons, do dot use QtWidgets.QDialogButtonBox(self),
    self.buttonNames = ["Apply", "Reset", "Cancel", "Quit", "Help"]
    self.buttonToolTips = ["Apply modifications and close dialog",
                           "Reload all initial values",
                           "Close dialog with all initial values",
                           "Close dialog",
                           "Help for widget"]
    self.buttons = {}
    for name, toolTip in zip(self.buttonNames, self.buttonToolTips):
      b = QtWidgets.QPushButton(name)
      b.setToolTip(toolTip)
      self.buttons[name] = b
    hBoxlayout  = QtWidgets.QHBoxLayout()
    for name in self.buttonNames:
      hBoxlayout.addWidget(self.buttons[name])

    # upWidget to contents future user widget with setUpWidgetLayout(userWidget)
    self.upWidget = QtWidgets.QWidget()
    self.upWidget.setLayout(QtWidgets.QVBoxLayout())
    self.upWidget.layout().setContentsMargins(m,m,m,m*2) # left, top, right, bottom

    vBoxlayout.addWidget(self.upWidget)    
    vBoxlayout.addLayout(hBoxlayout)

    self.choice = None
    self.buttons["Apply"].clicked.connect(self.on_apply)
    self.buttons["Reset"].clicked.connect(self.on_reset)
    self.buttons["Cancel"].clicked.connect(self.on_cancel)
    self.buttons["Quit"].clicked.connect(self.on_cancel)
    self.buttons["Help"].clicked.connect(self.on_help)

    self.hideButtons("Quit".split())
    self._UpWidget = None # as user appended widget (later), empty if not setUpWidgetLayout

  def hideButtons(self, namesButtons):
    for name in namesButtons:
      self.buttons[name].hide()
    return

  def showButtons(self, namesButtons):
    for name in namesButtons:
      self.buttons[name].show()
    return

  def setController(self, controller):
    if self._controller != None:
      mess = "\n"
      for i in traceback.format_stack():
        mess += (" | %s\n" % i)
      mess += "controller '%s'<->'%s' set yet for dialog '%s'\n" % \
       (self._controller.objectName(), controller.objectName(), self.objectName())
      logger.error(mess)
      return
    self._controller = controller
    
  def getController(self):
    return self._controller
    
  def on_apply(self):
    if verboseEvent: logger.info("DialogXyz.Apply")
    self.choice = "Apply"
    self.close()
   
  def on_cancel(self):
    if verboseEvent: logger.info("DialogXmlXyz.Cancel")
    self.choice = "Cancel"
    controller =  self.getController()
    if controller != None:
      controller.cancelDialogBox(self)
    else:
      self.close()

  def on_reset(self):
    if verboseEvent: logger.info("DialogXyz.Reset")
    self.choice = "Reset"
    self.resetUpWidget()
    
  def on_help(self):
    if verboseEvent: logger.info("DialogXyz.Help")
    self.helpUpWidget()
    
  def resetUpWidget(self):
    try:
      self._UpWidget.resetValue()
    except:
      logger.warning("no method resetValue() for the up widget '%s'" % self._UpWidget.__class__.__name__)

  def helpUpWidget(self):
    try:
      self._UpWidget.on_help()
    except:
      logger.warning("no method on_help() for the up widget '%s'" % self._UpWidget.__class__.__name__)

  def setUpWidgetLayout(self, oneOrSomeWidgets):
    if verboseEvent: logger.info("DialogXyz.setUpWidgetLayout")
    if type(oneOrSomeWidgets) == list:
      for w in oneOrSomeWidgets:  
        self.upWidget.layout().addWidget(w)
      self._UpWidget = oneOrSomeWidgets[0]
    else:
      self.upWidget.layout().addWidget(oneOrSomeWidgets)
      self._UpWidget = oneOrSomeWidgets
    return
    
  def getUpWidgetLayout(self):
    return self._UpWidget
    
  def xx_event(self, event):
    if verbose: logger.info('there is event: %s' % strEvent(event))
    return super(DialogXyz, self).event(event)

  def receiveRequestToView(self, strXmlRequest):
    if verboseEvent: 
      logger.info("%s %s receiveRequest virtual method" % (self.__class.__name__, self.objectName()))
    return True


##########################################################
class DialogXmlXyz(QtWidgets.QDialog):
  """
  a vertical layout with a upWidget with a QVBoxLayout
  and at bottom: Apply, Reset, Cancel buttons
  """
  _no = 0

  def __init__(self, *args, **kwargs):
    super(DialogXmlXyz,self).__init__(*args, **kwargs)
    noc = self._no
    self.setObjectName("DialogXmlXyz_%i" % noc)
    self._className = str(self.__class__.__name__)
    self._no += 1
    self._controller = None #caller controller general model
    self._localController = None #self controller local partial model
    vBoxlayout = QtWidgets.QVBoxLayout()
    self.setLayout(vBoxlayout)
    self.setWindowTitle(self.objectName())
    self.buttonNames = ["Apply", "Reset", "Cancel"]
    self.buttonToolTips = ["Apply modifications and close dialog",
                           "Reload all initial values",
                           "Close dialog with all initial values"]
    self.buttons = {}
    for name, toolTip in zip(self.buttonNames, self.buttonToolTips):
      b = QtWidgets.QPushButton(name)
      b.setToolTip(toolTip)
      self.buttons[name] = b
      
    hBoxlayout  = QtWidgets.QHBoxLayout()
    self.upWidget = QtWidgets.QWidget()
    self.upWidget.setLayout(QtWidgets.QVBoxLayout())
    self._margin = 2
    m = self._margin
    self.upWidget.layout().setContentsMargins(m,m,m,m) #left, top, right, bottom
    vBoxlayout.addWidget(self.upWidget)
    for name in self.buttonNames:
      hBoxlayout.addWidget(self.buttons[name])
    
    vBoxlayout.addLayout(hBoxlayout)
    self.choice = None
    self.dictForResetModel = None
    self.dictForResetEditors = None
    self.dictForHideEditors = None
    self.buttons["Apply"].clicked.connect(self.on_apply)
    self.buttons["Reset"].clicked.connect(self.on_reset)
    self.buttons["Cancel"].clicked.connect(self.on_cancel)
    self.closeOnApply = False #TODO maybe could be always True as DialogXyz... mefiance for ReleaseIra

  def setController(self, controller):
    if self._controller != None:
      mess = "\n"
      for i in traceback.format_stack():
        mess += (" | %s\n" % i)
      mess += "controller '%s'<->'%s' set yet for dialog '%s'\n" % \
       (self._controller.objectName(), controller.objectName(), self.objectName())
      logger.error(mess)
      return
    self._controller = controller
    
  def getController(self):
    return self._controller
    
  def getLocalModel(self):
    """reference, and direct modifications allowed, no total MVC pattern yet"""
    return self._localController._model
    
  def on_apply(self):
    if verboseEvent: logger.info("DialogXmlXyz.Apply")
    if verbose: logger.info("on_apply localModel:\n%s" % self.getLocalModel().toStrXml())
    self.choice = "Apply"
    controller =  self.getController()
    if controller != None:
      controller.applyDialogBox(self)
    if self.closeOnApply: 
      self.close()
    
  def on_cancel(self):
    if verboseEvent: logger.info("DialogXmlXyz.Cancel")
    self.choice = "Cancel"
    controller =  self.getController()
    if controller != None:
      controller.cancelDialogBox(self)
    else:
      self.close()

  def strDict(self, aDict):
    """str(aDict) with sorted keys"""
    #something like: "".join("\n"+[str((k, self.dictForResetEditors[k])) for k in sorted(self.dictForResetEditors.keys())])
    res = "{"
    sl = ""
    for k in sorted(aDict.keys()):
      res += "\n '%s': '%s'" % (k, aDict[k])
      sl = "\n"
    res += sl+"}"
    return res

  def on_reset(self):
    if verboseEvent: logger.info("DialogXmlXyz.Reset")
    self.choice = "Reset"
    self.resetEditors()
    self.resetLocalModel()
    #reset controller model? no! modif are only on local
    return
    
    self.clearUpWidgetLayout()
    aWidget = self.createDialogXyz_0("hello", self.getLocalModel(), [self.createDialogXyz_1, self.createDialogXyz_2])
    self.upWidget.layout().addWidget(aWidget)
    self.choice = "Reset"
    
  def resetLocalModel(self):
    if verboseEvent: logger.info("DialogXmlXyz.resetLocalModel")
    self.getLocalModel().setValueFromDictOfTreePyNameStrValues(self.dictForResetModel)
  
  def resetEditors(self):
    if verboseEvent: logger.info("DialogXmlXyz.resetEditors %s", self.strDict(self.dictForResetEditors))
    self._resetEditorsFromDict(self.dictForResetEditors)
    
  def _resetEditorsFromDict(self, aDict):
    #sorted as root(s) firsts, leaf last
    for k, (qeditor, strvalue) in sorted(aDict.items()):
      #print "resetEditor",k, strvalue
      qeditor.setValue(strvalue)
    return
  
  def _getDictForSetEditors(self):
    res = {}
    model = self.getLocalModel()
    for k, (qeditor, strvalue) in sorted(self.dictForResetEditors.items()):
      cmd = "newValue = model%s" % k
      namespace = {"model": model}
      namespace = UXYZ.exec2or3(cmd, namespace)  # TODO python 2-3
      newValue = namespace["newValue"]
      if verboseEvent:
        logger.info("setEditor %s %s %s" % (k, strvalue, newValue))
      res[k] = (qeditor, newValue)
    return res
  
  def refreshModel(self):
    """warning tree structure is not supposed to change, only leaf"""
    if verboseEvent: logger.info("DialogXmlXyz.refreshModel")
    aDict = self._getDictForSetEditors()
    self._resetEditorsFromDict(aDict)

  def resetUpWidgetLayout(self):
    if verboseEvent: logger.info("DialogXmlXyz.resetUpWidgetLayout")
    self.clearUpWidgetLayout()
    aWidget = self.createDialogXyz_0("hello", self.getLocalModel(), [self.createDialogXyz_1, self.createDialogXyz_2])
    self.upWidget.layout().addWidget(aWidget)
    return
    
  def clearUpWidgetLayout(self):
    """do not clear buttons"""
    """
    http://stackoverflow.com/questions/10416582/replacing-layout-on-a-qwidget-with-another-layout
    That will reparent all the child widgets to that temporary object, 
    and that object is deleted immediately along with its new children 
    because we don't keep a reference to it.
    
    If layout is the layout on a different widget, 
    setLayout() will reparent the layout and
    make it the layout manager for this widget
    """
    
    QtWidgets.QWidget().setLayout(self.upWidget.layout())
    self.upWidget.setLayout(QtWidgets.QVBoxLayout())
    m = self._margin
    self.upWidget.layout().setContentsMargins(m,m,m,m) #left, top, right, bottom
    return
  
  def clearLayout(self):
    """clear buttons also"""
    QtWidgets.QWidget().setLayout(self.layout())
    self.setLayout(QtWidgets.QVBoxLayout())
    return

  def activeStretch(self):
    """valid only once, no desactive stretch yet.."""
    self.layout().addStretch(1)

  def setFromXml(self, data=None, modeView="standard", treePyNameRoot=""):
    """
    modeView == "standard": #first list in tabview, other in groups
    modeView == "allInTab": #all lists in tabview, no groups
    """
    logger.debug("dialogXmlXyz.setFromXml type(data) %s" % type(data))
    if data == None: # an example...
      # print "dialogXmlXyz.setFromXml None data as an example"
      self.upWidget.layout().addWidget(TabWidgetXyz([ExampleWid() for i in range(4)]))
      return

    try:
      dataXml = ET.fromstring(data) #from str Xml or ElementTree
    except:
      dataXml = data # supposedly ET... yet

    # logger.info("dialogXmlXyz.setFromXml type(dataXml) is %s" % type(dataXml))
    
    self.clearUpWidgetLayout()
    try:
      localModel = UXYZ.fromXml(dataXml)
    except:
      logger.critical("dialogXmlXyz.setFromXml unexpected type(dataXml) %s" % type(dataXml))
      return # do nothing
    
    #local controller and garbage collecting
    self.treePyNameRoot = treePyNameRoot
    self._localController = ControllerXyz() #TODO kwarg desktop = desktop?
    self._localController.setModel(localModel)
    if debug: logger.info("self.getLocalModel() inital %s %s" % (self.treePyNameRoot,type(self.getLocalModel())))

    #self.localModel.setController(ctrl) #done in setModel
    self._localController.setViewLocal(self)
    
    self.dictForResetModel = localModel.getDictOfTreePyNameStrValues()
    self.dictForResetEditors = {}
    self.dictForHideEditors = {}

    key = dataXml.tag
    self.setObjectName(key)
    self.setWindowTitle(key)
    
    if debug:
      logger.info("create dialog widget from data:\n%s" % localModel.toStrXml())
      logger.info("controller of local model %s" % localModel.getController())
    
    if verbose:
      logger.info("name %s type %s" % (key, type(localModel)))

    #localRootModel = localModel
    localRootModel = localModel.getValueByTreePyName(self.treePyNameRoot)
    if modeView == "standard":
      #first list in tabview, other in groups
      aWidget = self.createDialogXyz_0(key, localRootModel, [self.createDialogXyz_1, self.createDialogXyz_2, self.createDialogXyz_0])
    elif modeView == "allInTab":
      #all lists in tabview, no groups
      logger.debug("DialogXmlXyz allInTab")
      aWidget = self.createDialogXyz_0(key, localRootModel, [self.createDialogXyz_0])
    self.upWidget.layout().addWidget(aWidget)
    #self.upWidget.layout().addStretch(1)

  def getEnd(self, aList):
    """get end of list without first item: alist[1:], or aList if len(aList)==1"""
    if len(aList) == 1: #same as last
      return aList
    else:
      return aList[1:]
  
  def createDialogXyz_debug(self, name, model, listOfRecursiveCalls):
    """in fine, for debug, we could create QRadioButton"""
    logger.info("createDialogXyz_debug %s" % name)
    res = QtWidgets.QRadioButton(name)
    res.setObjectName(name)
    return res
  
  def createDialogXyz_elementary(self, name, model, listOfRecursiveCalls): 
    """
    create generic dialog widget for elements with labels and editors
    leaf ot tree as immutables and other stuff
    """
    if verbose: logger.info("createDialogXyz_elementary %s %s" % (name,model.__class__.__name__))
    """
    res = QtWidgets.QWidget()
    res.setLayout(QtWidgets.QFormLayout())
    editor = model.createEditor(self)
    res.layout().addRow(name, editor)
    """
    res = QtWidgets.QWidget()
    res.setObjectName(name)
    res.setLayout(QtWidgets.QHBoxLayout())
    qeditor = model.createEditor(self)
    
    """
    if qeditor == None: #invalid qeditor
      qname = QtWidgets.QLabel(name.ljust(15))
      res.layout().addWidget(qname)
      strmodel = str(model)
      if len(strmodel) > 25:
        short = strmodel[0:25] + "..."
      else:
        short = strmodel
      label = QtWidgets.QLabel(short)
      #TODO for test to get contextmenu of leaf
      label = MyQLabel(short) ; label.setModelItem(model)
      label.setToolTip("invalid editor:\n" + strmodel)
      res.layout().addWidget(label)"""
    
    """
    if qeditor == None: #invalid qeditor
      qname = QtWidgets.QLabel(name.ljust(15))
      res.layout().addWidget(qname)
      strmodel = str(model)
      label = IFLX.XyzQLineEditAsQLabel(strmodel)
      label.setModelItem(model)
      label.setToolTip("invalid editor: try contexMenu\n" + strmodel)
      res.layout().addWidget(label)
      self.dictForResetEditors[model.getTreePyName()] = (qeditor, str(model)) #remind initial value in 
    """
    
    if qeditor == None:
      qeditor = IFLX.XyzQLineEditAsQLabel() #no modification
      #print "setModelItem(model) controller", model.getController(), model
      #qeditor.setToolTip("invalid editor: try contexMenu")
      
    if qeditor == None:
      logger.error("%s.createEditor None unexpected on value:\n%s" % (model._className, str(model).ljust(15)))
      return res
      
    if type(qeditor) == str:
      logger.warning("%s.createEditor: no edition allowed on value:\n%s" % (model._className, str(model)))
      #QtWidgets.QMessageBox.warning(None, "warning", qeditor)
      qeditor = QtWidgets.QLabel(str(model).ljust(15))
    
    else: #valid qeditor
      qeditor.setValue(model)
      #only for context menu powered editors (ex: XyzQLineEditBrowseFile)
      if hasattr(qeditor, "setModelItem"): qeditor.setModelItem(model)
      self.dictForResetEditors[model.getTreePyName()] = (qeditor, str(model)) #remind initial value in 

      # python 2 pyqt5
      try:
        qeditor.editingFinished.connect(lambda status=None, name=name, modelItem=model.getTreePyName(): self.somethingClicked(status, name, modelItem))
      except:
        try: 
          qeditor.currentIndexChanged.connect(lambda status=None, name=name, treePyName=model.getTreePyName(): self.somethingClicked(status, name, treePyName))
        except:
          vals = [i for i in dir(qeditor) if 'ed' == i[-2:]] #for debug info only ending "ED": finishED, editED modifiED etc
          logger.error("no existing event editingFinished or currentIndexChanged %s %s :\n%s" % (name, qeditor, PP.pformat(vals)))

      # python 2 pyqt4
      """
      try: 
        qeditor.editingFinished.connect(lambda name=name, modelItem=model.getTreePyName(): self.somethingClicked(name, modelItem))
      except:
        try: 
          qeditor.currentIndexChanged.connect(lambda name=name, treePyName=model.getTreePyName(): self.somethingClicked(name, treePyName))
        ...
      """

        
    qname = QtWidgets.QLabel(name.ljust(15))
    res.layout().addWidget(qname)
    res.layout().addWidget(qeditor)
    self.dictForHideEditors[model.getTreePyName()] = (qeditor, qname)
    hidden = model.isHidden()
    qeditor.setHidden(hidden)
    qname.setHidden(hidden)
      
    m = self._margin
    res.layout().setContentsMargins(m,m,m,0)
    return res
    
  def issubclassInList(self, model, aListOfClasses):
    for aClass in aListOfClasses:
      if issubclass(model.__class__, aClass): return True
    return False
    
  def createDialogXyz_0(self, name, model, listOfRecursiveCalls):
    """in first we create tabs in one tabwidget"""
    if verbose: logger.info("createDialogXyz_0 %s %s" % (name, model.__class__.__name__))
    
    if not self.issubclassInList(model, [BXYZ.BaseFreeXyz, BXYZ.ListOfBaseXyz]):
      #is a simple leaf value
      res = self.createDialogXyz_elementary(name, model, listOfRecursiveCalls)
      return res
    
    else:
      tabs = []
      others = []
      res = QtWidgets.QWidget()
      res.setObjectName(name)
      res.setLayout(QtWidgets.QVBoxLayout())
      m = self._margin
      res.layout().setContentsMargins(m,m,m,0) #left, top, right, bottom
      loops = [i for i in model.itemsByNumero()]
      logger.debug("model.itemsByNumero()\n%s" % PP.pformat(loops))
      for key, value in loops:
        if verbose: logger.info("createDialogXyz_0 key %s" % key)
        help = HLFX.getCommonHelp(model._className, key)
        aWidget = listOfRecursiveCalls[0](key, value, self.getEnd(listOfRecursiveCalls))
        if help is not None: 
          if verboseHelp: logger.info('help_0 %s' % help.shortHelp)
          aWidget.setToolTip(help.shortHelp)
        if self.issubclassInList(value, [BXYZ.BaseFreeXyz, BXYZ.ListOfBaseXyz]):
          if verbose: logger.info("createDialogXyz_0 set in tabs: %s" % key)
          tabs.append(aWidget)
        else:
          if verbose: logger.info("createDialogXyz_0 set in others: %s" % key)
          others.append(aWidget)
      #others top of widget
      for o in others: res.layout().addWidget(o)
      #TabWidget bottom of widget
      if len(tabs) > 0: res.layout().addWidget(TabWidgetXyz(tabs))
      res.layout().setSpacing(0)
      res.layout().addStretch(1)
      return res

  def createDialogXyz_1(self, name, model, listOfRecursiveCalls):
    """in first we create tabs in one tabwidget"""
    if verbose: logger.info("createDialogXyz_1 %s %s" % (name,model.__class__.__name__))
    
    if not self.issubclassInList(model, [BXYZ.BaseFreeXyz, BXYZ.ListOfBaseXyz]):
      #is a simple leaf value
      res = self.createDialogXyz_elementary(name, model, listOfRecursiveCalls)
      return res
    
    else:
      widgets = []
      res = QtWidgets.QWidget()
      res.setObjectName(name)
      res.setLayout(QtWidgets.QVBoxLayout())
      for key, value in model.itemsByNumero():
        #print "  createDialogXyz_1 key", key
        help = HLFX.getCommonHelp(model.__class__.__name__, key)
        aWidget = listOfRecursiveCalls[0](key, value, self.getEnd(listOfRecursiveCalls))
        if help is not None: 
          if verboseHelp: logger.info('help_1 %s' % help.shortHelp)
          aWidget.setToolTip(help.shortHelp)
        widgets.append(aWidget)
      for w in widgets: res.layout().addWidget(w)
      res.layout().addStretch(1)
      return res
   
  def createDialogXyz_2(self, name, model, listOfRecursiveCalls):
    """in second we create QGroupBoxes in one vlayout"""
    if verbose: logger.info("createDialogXyz_2 %s %s" % (name,model.__class__.__name__))
    
    if not self.issubclassInList(model, [BXYZ.BaseFreeXyz, BXYZ.ListOfBaseXyz]):
      #is a simple leaf value
      res = self.createDialogXyz_elementary(name, model, listOfRecursiveCalls)
      return res
    
    else: #create a groupbox
      #print "    model.__class__",model.__class__
      widgets = []
      res = QtWidgets.QGroupBox(name)
      res.setObjectName(name)
      res.setLayout(QtWidgets.QVBoxLayout())
      for key, value in list(model.items()):
        help = HLFX.getCommonHelp(model.__class__.__name__, key)
        aWidget = listOfRecursiveCalls[0](key, value, self.getEnd(listOfRecursiveCalls))
        #print "++createDialogXyz_2 widgets.append",key,type(aWidget)
        if help is not None: 
          if verboseHelp: logger.info('help_2 %s' % help.shortHelp)
          aWidget.setToolTip(help.shortHelp)
        widgets.append(aWidget)
      #widgets in groupbox
      for w in widgets: res.layout().addWidget(w)
      res.layout().addStretch(1)
      return res
      
  # python 2 pyqt5
  def somethingClicked(self, status=None,  name=None, treePyName=None):
    if verboseEvent: 
      logger.info("DialogXmlXyz.somethingClicked '%s' with %s" %  (treePyName, self.sender().__class__.__name__))
    try:
      sender = self.sender()
      newValue = self.sender().getValue()
    except:
      mess = "problem reaching result of editor: %s.getValue() inexisting method" % sender.__class__.__name__
      logger.error(mess)
      # logger.debug("dir(sender):\n%s" % PP.pformat(dir(sender)))

      QtWidgets.QMessageBox.warning(self, "error", mess)
      return
    
    try:
      oldValue = self.getLocalModel().getValueByTreePyName(treePyName)
      if verbose: 
        logger.info("DialogXmlXyz.somethingClicked on %s %s '%s' '%s' -> '%s'" % (name, str(sender.__class__.__name__), treePyName, oldValue, newValue))
    except:
      mess = "problem reaching local model item:\n '%s'" % treePyName
      QtWidgets.QMessageBox.warning(self, "error", mess)
      return
    
    try:
      self.getLocalModel().setValueByTreePyName(treePyName, newValue)
      if debug: logger.info("self.getLocalModel() after setValue %s" % type(self.getLocalModel()))
    except Exception as e:
      mess = "problem setting local model item:\n '%s' '%s' -> '%s'\n%s" % (treePyName, oldValue, newValue, e)
      QtWidgets.QMessageBox.warning(self, "error", mess)
      return
    #set hide or not...
    localModel = self.getLocalModel()
    for treePyName, (qeditor, qname) in list(self.dictForHideEditors.items()):
      #print "somethingClicked: isHidden",treePyName
      modelItem = localModel.getValueByTreePyName(treePyName)
      hidden = modelItem.isHidden()
      qeditor.setHidden(hidden)
      qname.setHidden(hidden)
    return

  def xx_event(self, event):
    if verbose: logger.info("DialogXmlXyz.event: there is event: %s" % strEvent(event))
    return super(DialogXmlXyz, self).event(event)

  def receiveRequestToView(self, strXmlRequest):
    if verboseEvent: 
      logger.info("%s %s receiveRequest method virtual" % (self.__class.__name__, self.objectName()))
    return True

