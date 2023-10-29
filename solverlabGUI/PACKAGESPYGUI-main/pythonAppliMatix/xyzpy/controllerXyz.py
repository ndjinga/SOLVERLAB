#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# see PACKAGESPY/LICENCE file
# %% LICENSE_END


import os
import traceback
import fnmatch
from datetime import datetime
import pprint as PP

from PyQt5 import QtGui, QtCore, QtWidgets
import xml.etree.ElementTree as ET
import xyzpy.loggingXyz as LOG
import xyzpy.intFloatListXyz as IFLX #append factory classes
import xyzpy.utilsXyz as UXYZ
import xyzpy.baseXyz as BXYZ
import salomepy.iconsUser as IUSR


logger = LOG.getLogger()
verbose = False
verboseEvent = False
debug = False

########################################################################################
class ControllerXyz(QtCore.QObject):
  """
  goal is a base class for manage request and action to/from views and model
  as MVC pattern
  request data through QtCore.pyqtSignal signals events, could be strings Xml, not useful (2017)
  one controller for one model and some views
  """
  index = [0]
  requestToModelSignal = QtCore.pyqtSignal(str) #when something have to do in model
  requestToViewsSignal = QtCore.pyqtSignal(str) #when something have to do in views (all views from one model)
  PrintApiIpcSignal = QtCore.pyqtSignal(object) #IPC interprocess communication
  refreshViewsSignal = QtCore.pyqtSignal()
  modelChangeSignal = QtCore.pyqtSignal()
  setModelItemValueSignal = QtCore.pyqtSignal(str)
  setModelItemValueSignalList = QtCore.pyqtSignal(list)
  ExpandSignal = QtCore.pyqtSignal(list)
  CollapseSignal = QtCore.pyqtSignal(list)

  def __init__(self, *args, **kwargs):
    """
    ControllerXyz(..., desktop=None, ...) #from salome desktop or standalone
    """
    super( ControllerXyz, self).__init__(*args)
    self._className = self.__class__.__name__ #shortcut
    self._MODULEDesktop = self.parent()
    self._views = {}
    self.treeViews = []
    self._dialogBoxes = {}
    self._model = None
    self._desktop = kwargs.get('desktop', None) #salome or qmainwindow?
    objectName = self.__class__.__name__+str(self.index)
    self.index[0] += 1 #unambigous objectName
    self.setObjectName(objectName)
    self.prefixShortcut = "Ctrl+"

    #IPC interprocess communication
    self._IpcController = None
    self.PrintApiIpcSignal.connect(self.printApiIpc) #IPC interprocess communication

    self.refreshViewsSignal.connect(self.refreshModelViewsFromSignal)
    self.setModelItemValueSignal.connect(self.setModelItemValueFromSignal)
    self.setModelItemValueSignalList.connect(self.setModelItemValueFromSignalList)
    self.ExpandSignal.connect(self.setExpand)
    self.CollapseSignal.connect(self.setCollapse)

    #contains links to environment variables and used directories and default editor command etc...
    self._settings = None
    #from namespace from parser
    self._args = None

  def ExpandAllAction(self):
    if verboseEvent: print("ExpandAllAction")
    for treeView in self.treeViews:
      treeView.expandAll()
    return True

  def setExpand(self, aList):
    logger.debug("Expand\n%s" % PP.pformat(aList))
    for treeView in self.treeViews:
      treeView.ExpandItems(aList)
    return True

  def setCollapse(self, aList):
    logger.debug("Collapse\n%s" % PP.pformat(aList))
    for treeView in self.treeViews:
      treeView.CollapseItems(aList) # userExpandModel
    return True

  def getIpcController(self):
    return self._IpcController

  def setIpcController(self, controllerIpc):
    if self._IpcController != None:
      logger.error("Ipc controller set yet")
      return
    self._IpcController = controllerIpc
    controllerIpc.setControllerOfModel(self)

  def printApiIpc(self):
    controllerIpc = self.getIpcController()
    if controllerIpc == None:
      logger.warning("No Ipc controller: No API")
    else:
      controllerIpc.IpcPrintApi()

  def setSettings(self, settings):
    self._settings = settings

  def getSettings(self):
    return self._settings

  def setArgs(self, args):
    self._args = args

  def getArgs(self):
    return self._args

  def setSettings(self, settings):
    self._settings = settings

  def getSettings(self):
    return self._settings

  def getDesktop(self):
    if self._desktop == None:
      import salomepy.xsalomesession as XSS
      self._desktop = XSS.getDesktop()
      if self._desktop != None: # if not work (in test ?)
        logger.warning( "set default desktop named '%s' from xsalomesession" % (self._desktop.objectName()) )
      else:
        logger.warning( "set default desktop as None" )
    return self._desktop

  def getMODULEDesktop(self):
    return self._MODULEDesktop
    # res = self._MODULEDesktop  # case inside SALOME
    # if res is None:  # case outside SALOME as standalone
    #  return self._desktop

  def setModel(self, aModel):
    if self._model != None:
      raise Exception('ControllerXyz.setModel done yet in %s' % self.objectName())
       #TODO self.requestToModelSignal.disconnect(aModel.receiveRequest)
    logger.debug("ControllerXyz.setModel")
    self._model = aModel
    aModel.setController(self)
    #link asynchronous requests to model
    self.requestToModelSignal.connect(aModel.receiveRequest)
    #ask to views to refresh
    aRequest = self.getRequest("newModel")
    self.sendRequestToViews(aRequest)
    self.modelChangeSignal.emit()

  def replaceModel(self, aModel):
    logger.debug("ControllerXyz.replaceModel %s" % type(aModel))
    self._model = aModel
    aModel.setController(self)
    #link asynchronous requests to model
    self.requestToModelSignal.connect(aModel.receiveRequest)
    #ask to views to refresh
    aRequest = self.getRequest("newModel")
    self.sendRequestToViews(aRequest)
    self.modelChangeSignal.emit()

  def clearModel(self):
    logger.debug("ControllerXyz.clearModel")
    if self._model != None:
      self._model = None
    #action systematic... for views
    aRequest = self.getRequest("clearModel")
    self.sendRequestToViews(aRequest)
    self.modelChangeSignal.emit()

  def setView(self, aView):
    viewName = aView.objectName() #unambigous name as an ident
    logger.debug("controller '%s' setView '%s'" %  (self.objectName(), viewName))
    if viewName in list(self._views.keys()):
      raise Exception("'%s' '%s' setView done yet for view '%s'" % (self.__class__.__name__, self.objectName(), viewName) )
    try:
      logger.debug("'%s' '%s' setView '%s' type '%s'" % (self.__class__.__name__, self.objectName(), viewName, aView.__class__.__name__))
      aView.setController(self)
      try:
        self.requestToViewsSignal.connect(aView.receiveRequestToView)
      except:
        #problems for controller request signals etc...
        logger.error("controller.setView of view type '%s' without receiveRequestToView method" % aView.__class__.__name__)
    except:
      #problems for controller etc...
      logger.error("controller.setView of view type '%s' without setController method" % aView.__class__.__name__)
    self._views[viewName] = aView

  def setViewLocal(self, aView):
    """for dialogXmlXyz, TODO better MVC"""
    viewName = aView.objectName() #unambigous name as an ident
    logger.debug("controller '%s' setView '%s'" %  (self.objectName(), viewName))
    if viewName in list(self._views.keys()):
      raise Exception("'%s' '%s' setView done yet for view '%s'" % (self._className, self.objectName(), viewName) )
    try:
      logger.debug("'%s' '%s' setView '%s' type '%s'" % (self._className, self.objectName(), viewName, aView.__class__.__name__))
      #aView.setController(self) #set to general controller
      try:
        self.requestToViewsSignal.connect(aView.receiveRequestToView)
      except:
        #problems for controller request signals etc...
        logger.error("controller.setView of view type '%s' without receiveRequestToView method" % aView.__class__.__name__)
    except:
      #problems for controller etc...
      logger.error("controller.setView of view type '%s' without setController method" % aView.__class__.__name__)
    self._views[viewName] = aView

  def getModel(self):
    """
    model is only owned by controller. MVC pattern.
    so return model is a copy as Elementtree
    """
    if self._model == None:
      logger.debug("controller.getModel: no model")
      return None
    return self._model.toXml()

  def getViews(self):
    views = []
    for k in sorted(self._views.keys()):
      views.append(self._views[k])
    return views

  def getExploreDir(self):
    try:
      return self.centralLogView.getTabByName("Explore Dir")
    except:
      return None

  def showTab(self, index):
    self.centralLogView.showTab(index)

  def showTabByName(self, name):
    self.centralLogView.showTabByName(name)

  def showExploreDir(self):
    return self.centralLogView.showTabByName("Explore Dir")

  def getIncrement(self):
    """get an unique once shot str(something), based on datetime.now()"""
    return self.getDateTimeNow()

  def getDateTimeNow(self):
    """get an str(), based on datetime.now()"""
    now = datetime.now()
    #return now.strftime("%Y-%m-%d %H:%M")
    #return now.strftime("%m%d%H%M%S")
    return now.strftime("%y%m%d_%H%M%S")

  def makedirs(self, namedir):
    if os.path.exists(namedir):
      dirbak=namedir+self.getIncrement()
      if os.path.exists(dirbak): shutil.rmtree(dirbak)
      os.rename(namedir,dirbak)
      os.listdir(dirbak) #sert seulement a mettre a jour le systeme de fichier sur certaines machines
    os.makedirs(namedir)

  def isSameFile(self, file1, file2):
    """follow relative paths and symbolic links"""
    rfile1=os.path.realpath(file1)
    rfile2=os.path.realpath(file2)
    return (rfile1==rfile2)

  def getRequest(self, nameRequest):
    """
    goal is view want to request via controller something to model.
    view get a request from controller, and customize it:

    usage:

    >>> #get a XyzControllerRequest instance from factory of controller
    >>> aRequest = controller.getRequest("modification")
    >>> aRequest.xmlModel = self.xmlData
    >>> controller.sendRequestToController(aRequest)   #aRequest have to be const
    >>> #return of this method is typed as simple BaseFreeXyz
    >>> #definition example follow in source
    """

    aRequest = BXYZ.BaseFreeXyz()
    aRequest.typeRequest = IFLX.StrXyz(nameRequest)
    aDatetime = datetime.now()
    aRequest.strdatetime = IFLX.StrXyz(aDatetime.strftime("%Y_%m_%d-%H_%M_%S"))

    """
    #caller have to add his own model as xmlModel and may be other data...
    #aRequest.xmlModel = IFLX.NoneXyz() #TODO could be not set at all if typed aRequest (future)
    #aRequest.xmlModel = IFLX.StrXyz("an unknown xmlModel to modify") #TODO could be not set at all if typed aRequest (future)
    """

    logger.debug("ControllerXyz.getRequest aRequest:\n%s" % aRequest.toStrXml())
    return aRequest

  def sendRequestToController(self, aRequest, verbose=False):
    """
    sender owned aRequest.xmlModel have to be const
    """
    logger.debug("ControllerXyz.sendRequestToController %s at %s" % (aRequest.typeRequest, aRequest.strdatetime))
    #TODO better aRequest.toXml with aRequest.xmlModel
    xmlRequest = UXYZ.toXml(aRequest)
    xmlRequest.tag = "RequestXyz"
    #xmlModel = ET.Element('xmlModel')
    #xmlModel.append(aRequest.xmlModel)
    #xmlRequest.append(xmlModel)
    #set xmlRequest as string Xml
    strXmlRequest = UXYZ.prettyPrintET(xmlRequest)
    #send signal to model... if necessary
    if self._model != None:
      self.requestToModelSignal.emit(strXmlRequest)
      return True
    if verbose: logger.warning("sendRequestToController : No model for controller %s" % (self.objectName()))
    return False

  def sendRequestToViews(self, aRequest):
    """
    sender owned aRequest.xmlModel have to be const
    """
    logger.debug("ControllerXyz.sendRequestToViews %s" % aRequest)
    #TODO better aRequest.toXml with aRequest.xmlModel
    xmlRequest = UXYZ.toXml(aRequest)
    xmlRequest.tag = "RequestXyz"
    #set xmlRequest as string Xml
    strXmlRequest = UXYZ.prettyPrintET(xmlRequest)
    #send signal to views... TODO if necessary ??
    self.requestToViewsSignal.emit(strXmlRequest)
    self.refreshViewsSignal.emit()
    return True

  def getSalomeStudy(self):
    """store result in self.salomeStudy"""
    if self._salomeStudy!=None:
      #if salome crash or GUI user exit
      try:
        logger.info("GetSalomeStudy: (done yet) got study: %s Id: %s",self._salomeStudy._get_Name(),self._salomeStudy._get_StudyId())
        return
      except:
        self._salomeStudy==None
    self._omniOrbConfigFile = os.getenv('OMNIORB_CONFIG')
    if self._omniOrbConfigFile==None:
      logger.error("I can't get a salome-cassis cession because OMNIORB_CONFIG environment variable is None.")
      self._salomeStudy = None
      return
    if not os.path.isfile(self._omniOrbConfigFile):
      logger.warning("I can't see a salome-cassis cession: no file " + self._omniOrbConfigFile)
      logger.warning("I try to launch a salome-cassis cession... take a time...")
      self.launchSalome()
    import salome
    for i in range(0,5): #try twice
      sleep(5)
      self._salomeStudy = None
      logger.info("try getSalomeStudy %s" ,i)
      try:
        salome.salome_init(theStudyId=1) #get the number Id 1, or create it
        self._salomeStudy = salome.myStudy
        logger.info("GetSalomeStudy: got study: %s Id: %s",self._salomeStudy._get_Name(),self._salomeStudy._get_StudyId())
        return
      except:
        pass
    logger.error("I can't get a salome cession. Launch one please: '.../salome7 -k'")

  def launchSalome(self):
    cmd = 'python -u -c "from salomeRunner import SalomeRunner ; runner = SalomeRunner(configFileNames=None) ; runner.go([])" &'
    self.launchSalomeIntoPopen(cmd)
    sleep(1) #with a background "&" the first sleep is killed ...mystery...
    self._salomeStudy = None

  def chmodarwx(self, nameFile):
    """chmod a+rwx"""
    import stat
    st = os.stat(nameFile)
    os.chmod(nameFile, st.st_mode | stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
    return

  def getDefaultContextMenuZeroSelection(self):
    """virtual method"""
    logger.warning("controllerXyz.getDefaulContextMenuZeroSelection standard method")
    return None

  def getDefaultContextMenuMultipleSelection(self, treePaths):
    """virtual method"""
    logger.warning("controllerXyz.getDefaulContextMenuMultipleSelection standard no actions for '%s'" % str(treePaths))
    return None

  def getDefaultContextMenuForItem(self, treePath):
    """controllerXyz.getDefaultMenuForItem method could be overriden in subclass '%s'" % self.__class__.__name__"""
    desktop = self.getDesktop()

    if verbose: print("controllerXyz %s.getDefaultContextMenuForItem for '%s'" % (self.__class__.__name__, str(treePath)))
    value = self._model.getValueByTreePyName(treePath)
    #value._controller = self #TODO better as root model than item

    if not hasattr(value, "getActionsContextMenu"):
      mess = "No context menu for '%s'\ntry direct edition (mouse double left click)" % treePath
      QtWidgets.QMessageBox.warning(desktop, "warning", mess)
      return None
    try:
      actions = value.getActionsContextMenu()
    except:
      traceback.print_exc() #better explicit verbose problem
      mess = "problem getActionsContextMenu for '%s'" % treePath
      QtWidgets.QMessageBox.warning(desktop, "warning", mess)
      return None

    # print("controllerXyz actions:\n" + "\n".join([str((a.ClassName, a.Name)) for a in actions]))
    if len(actions)==0: return None
    menu = QtWidgets.QMenu("ContextMenu%s" % value.__class__.__name__, desktop)
    for a in actions:
      action = self._createAction(a.Name, a.Shortcut, a.ToolTip, a.Call, a.Icon, a.Enable)
      menu.addAction(action)
    return menu

  def _createAction(self, Name, Shortcut, ToolTip, Call, Icon=None, Enable=True):
    #logger.debug("cvw _createAction 2 %s %s" % (self.__class__.__name__, Name, Icon))
    desktop = self.getDesktop()
    action = QtWidgets.QAction(Name, desktop)
    action.setIconVisibleInMenu(True)
    if Shortcut != None: action.setShortcut(self.prefixShortcut + Shortcut)
    action.setToolTip(ToolTip)
    if Icon != None:
      action.setIcon(IUSR.getIconFromName(Icon))
    action.setEnabled(Enable)
    action.triggered.connect(Call)
    return action

  def refreshModelViews(self):
    logger.debug("refreshModelViews")
    #aDataXml =  self._model.toXml() #as a copy, only for debug
    for view in self.getViews():
      if hasattr(view, "refreshModel"):
        logger.debug("refreshModelViews %s" % view.objectName())
        view.refreshModel()
    return True

  def refreshModelViewsFromSignal(self):
    #only for verboseEvent
    logger.debug("refreshModelViewsFromSignal")
    return self.refreshModelViews()

  def setModelItemValueFromSignal(self, aSet):
    model = self._model
    if model == None:
      logger.warning("controller.setModelItemValueFromSignal: no model to change")
      return
    cmd = "model%s" % aSet
    namespace = {"model": model}
    res = UXYZ.exec2or3(cmd, namespace)  # TODO python 2-3
    self.refreshModelViews()
    self.modelChangeSignal.emit()

  def setModelItemValueFromSignalList(self, args):
    """
    args from signal is list, args[0] is TreePyName expression for set, or else user trick.

    | example:
    | >>> controller.setModelItemValueSignalList.emit( ["%s.setSource(args[1])" % self.getTreePyName(), newSource ] )
    | where newSource is args[1] for exec() (yes it is weird..;)
    """
    model = self._model
    if model == None:
      logger.warning("controller.setModelItemValueFromSignal: no model to change")
      return
    cmd = "model%s" % args[0]
    logger.debug("controller.setModelItemValueFromSignalList cmd='%s' args=\n%s" % (cmd, PP.pformat(args)))
    namespace = {"model": model, "args": args}
    res = UXYZ.exec2or3(cmd, namespace)  # TODO python 2-3
    # logger.debug("result setModelItemValueFromSignalList %s" % model.Simulations[0].Variables)
    self.refreshModelViews()
    self.modelChangeSignal.emit()

  def assumeDialogBox(self, widDialog, treePyName):
    """to know an active DialogBox exists on itemModel"""
    logger.debug("controller.assumeDialogBox treePyName '%s'" % treePyName)
    itemModel = self._model.getValueByTreePyName(treePyName)
    if itemModel == None:
      logger.error("controller model do not contains item '%s'" % treePyName)
      return
    widDialog.setController(self)
    self._dialogBoxes[widDialog] = treePyName #persistence of widDialog
    #self._dialogBoxes[treePyName] = widDialog #persistence of widDialog
    widDialog.show() #TODO set in dockview if desktop salome ... or not

  def applyDialogBox(self, widDialog):
    """
    controllerXyz answer to a DialogXmlXyz apply event
    set widDialog.getLocalModel() to controller model with treePyName
    treePyName set from assumeDialogBox
    """
    try:
      treePyName = self._dialogBoxes[widDialog]
    except:
      logger.error("problem controller do not know widDialog")
      return
    try:
      itemModel = self._model.getValueByTreePyName(treePyName)
    except:
      logger.error("problem controller model with unknown item '%s'" % str(treePyName))
      return
    logger.debug("controller.applyDialogBox '%s'" % treePyName)
    if True:
      #replace all model as root
      #replace model tree
      newModel = UXYZ.fromXml(widDialog.getLocalModel().toXml())
      logger.debug("on_apply newModel:\n%s" % widDialog.getLocalModel().toStrXml())

    widDialog.close()
    del (self._dialogBoxes[widDialog])
    self.replaceModel(newModel)

  def obsolete_applyDialogBox(self, widDialog):
    """
    | controllerXyz answer to a DialogXmlXyz apply event
    | set widDialog.getLocalModel() to controller model with treePyName
    | treePyName set from assumeDialogBox
    """
    try:
      treePyName = self._dialogBoxes[widDialog]
    except:
      logger.error("problem controller do not know widDialog")
      return
    try:
      itemModel = self._model.getValueByTreePyName(treePyName)
    except:
      logger.error("problem controller model with unknown item '%s'" % str(treePyName))
      return
    logger.debug("controller.applyDialogBox '%s'" % treePyName)
    #print "old itemModel\n",itemModel
    #print "new itemModel\n",widDialog.getLocalModel()
    if treePyName != "":
      #change part of model tree
      newModel = UXYZ.fromXml(self._model.toXml())
      newValue = UXYZ.fromXml(widDialog.getLocalModel().toXml())
      #print "%%%%%% localModel",newValue.toStrXml()
      try:
        newModel.setValueByTreePyName(treePyName, newValue)
      except:
        logger.error("problem set new value '%s'" % str(treePyName))
        return #no close widDialog ???...
    else:
      #replace all model as root
      #replace model tree
      newModel = UXYZ.fromXml(widDialog.getLocalModel().toXml())

    widDialog.close()
    del (self._dialogBoxes[widDialog])
    self.replaceModel(newModel)

  def cancelDialogBox(self, widDialog):
    widDialog.close()
    del (self._dialogBoxes[widDialog])

  def getLogCMDWidget(self):
    """to get a trace output device as expected widget textedit"""
    try:
      return self.centralLogView.centralWidget().logCMDWidget
    except:
      return None

  def activateLogCMDWidget(self):
    """to show trace output device as widget textedit"""
    try:
      #show in ModuleLog viewer widget the 'Log Run Code'"""
      self.centralLogView.centralWidget().showLogCMDWidget()
      #show in salome desktop the 'ModuleLog viewer'
      self.getMODULEDesktop().activateLogView()
    except:
      logger.error("Problem to activate unreachable 'Log Run Code' widget")

  def showObjectBrowser(self):
    try:
      modDesk = self.getMODULEDesktop()
      modDesk.tryShow(modDesk.findDockByName("Object Browser"))
    except:
      logger.error("Problem to show unreachable 'Object Browser' widget")


  def quickEdit(self, fileName):
    """
    designed for quick display/edit a file
    with existing salome centralWidget WindowForLog
    and MODULEDesktop (as DARTDesktop or else...)
    """
    aFile = str(fileName)
    try: #with salome centralWidget WindowForLog
      self.centralLogView.centralWidget().quickEditFiles([aFile])
      self.getMODULEDesktop().activateLogView() #set front
    except:
      logger.error("%s: problem for quick edit on "% (self._className, fileName))

  def getDefaultView(self):
    viewNames = list(self._views.keys())
    if len(viewNames) == 1:
      return self._views[viewNames[0]]
    else:
      logger.error("%s: problem for getDefaultView on unknown view %s"% (self._className, str(viewNames)))
      return None

  def getDefaultTreeView(self):
    if len(self.treeViews) == 1:
      return self.treeViews[0]
    else:
      logger.error("%s: problem for getDefaultTreeView on unknown view %s"% (self._className, str(self.treeViews)))
      return None

  def getSelectedItems(self, onView=None):
    view = onView
    if view == None:
      view = self.getDefaultTreeView()
    return view.selectedItems()

  def getTreePathsSelected(self, onView=None):
    """good for pattern MVC because do not return elements of model"""
    view = onView
    if view == None:
      view = self.getDefaultTreeView()
    items = self.getSelectedItems(view)
    treePathsSelected = view._getTreePathsOfListItems(items)
    return treePathsSelected

  def getTreePathsSelectedOnPattern(self, pattern="*", onView=None):
    treepaths = self.getTreePathsSelected(onView)
    return fnmatch.filter(treepaths, pattern)

  def _getModelItemsFromTreePath(self, treePaths):
    """
    have to be private
    risky to return pointer to elements of model for MVC pattern,
    so use treepath instead (for external use)
    treePaths is a list
    """
    res = []
    for i in treePaths:
      namespace = {"myself": self}
      cmd = "ele = myself._model%s" % i
      namespace = UXYZ.exec2or3(cmd, namespace)  # TODO python 2-3
      res.append(namespace["ele"])
    return res
